##' @import arrow
NULL

## Shared helpers

.check_q <- function(cs) {
    stopifnot("cs must be a data.frame containing numeric columns chr and pos"=
                  is.data.frame(cs) &&
                  all(c("chr","pos") %in% names(cs)) &&
                  is.numeric(cs$chr) &&
                  is.numeric(cs$pos))
}

.check_window <- function(window, ot = FALSE) {
  stopifnot("window should be a single non-negative number" = is.numeric(window) && length(window) == 1 && window >= 0)
  if (ot && window > 5e+5)
    warning("You have set window very large. This is probably unnecessary and will cause downloads to be slower.")
}

## Resolve a data directory: argument, then option, then a per-session temp folder
.cache_dir <- function(dir, subdir, option="coloc.cache_dir") {
    if(is.null(dir))
        dir <- getOption(option)
  if(is.null(dir)) {
    warning("Downloading files to a temporary directory. Set `dir` or options(",option," = \"path/to/dir\") to keep downloads between R sessions. This will save you time by preventing unnecessary downloads in each new R session.")
    dir <- file.path(tempdir(), subdir)
  }
    if(!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
  normalizePath(dir)
}

## Download url to dir/file unless already there; returns the local path
## An existing file is reused, unless size is given and it doesn't match. A new download is checked
## against size and md5 (either may be NULL/NA); if either is wrong it is deleted and we stop.
.download_once <- function(url, dir, file = basename(url), md5 = NULL, size = NULL) {
  path <- file.path(dir, file)
  if (file.exists(path) && (is.null(size) || is.na(size) || file.size(path) == size))
    return(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  old <- options(timeout = max(3600, getOption("timeout"))) # default 60s is far too short for large files
  on.exit(options(old))
  message("Downloading ", file, " to ", dir)
  ## per-process temp file: an interrupted download never looks like a complete file, and parallel
  ## workers fetching the same file don't write into each other's downloads. rename is atomic, so the
  ## final file only ever appears complete; if two workers download it, the second replaces it with an identical copy
  tmp <- paste0(path, ".part.", Sys.getpid())
  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  bad <- if (!is.null(size) && !is.na(size) && file.size(tmp) != size) "Size" else
         if (!is.null(md5) && !is.na(md5) && unname(tools::md5sum(tmp)) != md5) "Checksum" else NULL
  if (!is.null(bad)) {
    unlink(tmp)
    stop(bad, " mismatch for ", file, "; the download may be corrupt or interrupted. Rerun to retry.", call. = FALSE)
  }
  file.rename(tmp, path)
  path
}

################################################################################
## Open Targets queries

.ot_release <- "26.06"
.ot_s3_region <- "eu-west-1" # bucket region; check with: curl -sI https://open-targets-public-data-releases.s3.amazonaws.com | grep -i x-amz-bucket-region
.ot_s3_bucket <- "open-targets-public-data-releases"
.ot_s3_url <- sprintf("https://%s.s3.%s.amazonaws.com", .ot_s3_bucket, .ot_s3_region) # plain https access to the public bucket
.ot_credset_prefix <- sprintf("platform/%s/output/credible_set/", .ot_release)
.ot_manifest_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/output/study/part-00000-1c0a4198-4d65-425e-82a1-d9d8b54d8549-c000.snappy.parquet", .ot_release)

## Study columns added to the query_ot() annotation, where present in the study file
## See https://platform.opentargets.org/downloads/study/schema
.ot_study_cols <- c("studyId", "projectId", "traitFromSource", "geneId", "biosampleFromSourceId", "pubmedId")

## Session cache for small lookup tables (file index, study manifest)
.ot_env <- new.env(parent = emptyenv())

.ot_dir <- function(dir = NULL) .cache_dir(dir, "opentargets")
.ot_local_dir <- function(dir = NULL) file.path(.ot_dir(dir), paste0("credible_set_", .ot_release))

## Index of the credible set parquet files: name, size, first/last (chromosome, position), and md5 (NA where unknown).
## The files are sorted by chromosome (as text: 1, 10, 11, ..., 2, 20, ...) then lead variant position,
## so each covers one contiguous block. Regenerate with data-raw/ot_credset_index.R for a new release.
.ot_index <- function() {
  if (is.null(.ot_env$index)) {
    f <- system.file("extdata", sprintf("ot_credset_index_%s.tsv", .ot_release), package = "coloc")
    if (!nzchar(f)) stop("Open Targets file index for release ", .ot_release, " not found in inst/extdata")
    .ot_env$index <- utils::read.delim(f, colClasses = c(chr_first = "character", chr_last = "character", md5 = "character"))
  }
  .ot_env$index
}

## Files whose (chromosome, position) range overlaps chr:lower-upper
.ot_files_for <- function(chr, lower, upper) {
  idx <- .ot_index()
  chrs <- unique(c(rbind(idx$chr_first, idx$chr_last))) # chromosome order used in the files
  if (!as.character(chr) %in% chrs) return(idx[0, ])
  key <- function(c, p) match(c, chrs) * 1e10 + p # sortable (chromosome rank, position)
  idx[key(idx$chr_last, idx$pos_last) >= key(as.character(chr), lower) & key(idx$chr_first, idx$pos_first) <= key(as.character(chr), upper), ]
}

## Stop if a large download would go to a temporary directory, where it would be repeated every session
.check_download_size <- function(bytes, dir, what) {
  limit <- getOption("coloc.max_temp_download", 1e9)
  if (is.null(dir) && is.null(getOption("coloc.cache_dir")) && bytes > limit)
    stop(what, " needs a ", round(bytes / 1e9, 1), " GB download, which is too large for a temporary directory. ",
         "Please supply a persistent cache directory with `dir`, or set options(coloc.cache_dir = \"path/to/dir\"), so it is only downloaded once.",
         call. = FALSE)
}

## Make sure the given credible set files are in the local cache, downloading any that are missing; returns local paths
.ot_fetch_files <- function(files, dir = NULL, what = "This Open Targets query") {
  local <- .ot_local_dir(dir)
  paths <- file.path(local, files$file)
  todo <- which(!file.exists(paths) | file.size(paths) != files$size)
  if (length(todo)) {
    .check_download_size(sum(files$size[todo]), dir, what)
    if (!dir.exists(local)) dir.create(local, recursive = TRUE)
    message("Downloading ", length(todo), " Open Targets credible set file(s) (", round(sum(files$size[todo]) / 1e6), " MB) to ", local)
    for (i in todo)
      .download_once(paste0(.ot_s3_url, "/", .ot_credset_prefix, files$file[i]), local, file = files$file[i], md5 = files$md5[i], size = files$size[i])
  }
  paths
}

##' Download a complete local copy of the Open Targets credible sets
##'
##' query_ot() downloads only the credible set files it needs for each
##' query (typically one or two files, ~15 MB each) and caches them. Use
##' ot_fetch() to download all of them (about 3.2 GB for release 26.06)
##' in one go, for example on an HPC login node before running
##' query_ot() on compute nodes without internet access, using the same
##' dir.
##'
##' Because of its size, a persistent dir (or
##' `options(coloc.cache_dir=...)`) is required. Files already present
##' with the expected size are not downloaded again, so an interrupted
##' download can be resumed by calling ot_fetch() again.
##'
##' @inheritParams query_ot
##' @return the path to the local credible set directory, invisibly
##' @export
##' @author Chris Wallace
##' @seealso [query_ot()]
ot_fetch <- function(dir = NULL) {
  .ot_fetch_files(.ot_index(), dir, what = "Downloading all Open Targets credible sets")
  invisible(.ot_local_dir(dir))
}

## Credible set dataset for chr:lower-upper, opened fresh on each call (no cached, long-lived network connections)
## remote = FALSE: download the needed files to the cache dir (if not already there) and read them locally
## remote = TRUE: read the needed files directly from S3, without downloading
.ot_dataset <- function(chr, lower, upper, dir = NULL, remote = FALSE) {
  files <- .ot_files_for(chr, lower, upper)
  if (!nrow(files)) return(NULL)
  if (isTRUE(remote)) {
    fs <- arrow::S3FileSystem$create(anonymous = TRUE, region = .ot_s3_region) # region set explicitly: automatic lookup fails on some HPC systems
    return(arrow::open_dataset(paste0(.ot_s3_bucket, "/", .ot_credset_prefix, files$file), filesystem = fs, format = "parquet")) # file paths within the bucket
  }
  arrow::open_dataset(.ot_fetch_files(files, dir), format = "parquet")
}

## Scan the credible sets; a dropped S3 connection can leave R unstable, so stop with advice rather than retry
.ot_scan <- function(ds, ...) {
  tryCatch(arrow::Scanner$create(ds, ...)$ToTable(), error = function(e) {
    if (grepl("SIGPIPE", conditionMessage(e)))
      stop("The Open Targets S3 connection was dropped (SIGPIPE). Please restart R before continuing, as the session may be unstable. Using remote = FALSE (the default) avoids S3.", call. = FALSE)
    stop(e)
  })
}

## OT manifest
.ot_manifest <- function(dir = NULL) {
    path <- .download_once(.ot_manifest_url, .ot_dir(dir))
    if (is.null(.ot_env$study))
        .ot_env$study <- arrow::open_dataset(path, format = "parquet")
  .ot_env$study
}

## Study annotations for a set of studyIds
.ot_study <- function(study_ids, dir = NULL, cols = .ot_study_cols) {
  ds <- .ot_manifest(dir)
  cols <- intersect(unique(c("studyId", cols)), names(ds))
  ids <- arrow::Array$create(unique(study_ids))
  filt <- arrow::Expression$create("is_in", arrow::Expression$field_ref("studyId"), options = list(value_set = ids, skip_nulls = TRUE))
  as.data.frame(arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable())
}

## Convert one Open Targets locus (struct of variants) to a credible set
.otlocus_to_cs <- function(locus, ...) {
  parts <- strsplit(locus$variantId, "_", fixed = TRUE) # split once, reuse
  data.table(...,
             ## chr = as.integer(vapply(parts, `[[`, character(1), 1)),
             chr = vapply(parts, `[[`, character(1), 1),
             pos = as.integer(vapply(parts, `[[`, character(1), 2)),
             snp = locus$variantId,
             PP = locus$posteriorProbability,
             beta = locus$beta)
}

##' query Open Targets credible sets
##'
##' This function takes an input credible set, uses that to define a
##' region of interest, and searches for Open Targets credible sets
##' that may overlap - ie contain at least one snp within the region
##' +/- window. The resulting list of Open Targets credible sets may
##' then be colocalised with the index credible set using the
##' coloc_query() function.
##'
##' The Open Targets credible sets are stored in 200 files, each
##' covering a contiguous block of the genome. By default, query_ot()
##' downloads only the file(s) overlapping the query region (typically
##' one or two, ~15 MB each) to the cache dir, and reuses them in later
##' queries. Alternatively, remote = TRUE reads those files directly
##' from the Open Targets S3 bucket without downloading. To download
##' everything up front (eg before working on HPC compute nodes without
##' internet access), use [ot_fetch()].
##'
##' S3 access (remote = TRUE) opens a fresh connection for each call.
##' Avoid using it inside parallel::mclapply or other forked workers. If
##' a call fails with a SIGPIPE error, restart R before continuing.
##'
##' @param cs an index data.frame containing columns "chr" "pos". NB
##'     chr should be a numeric: 1 for chromosome 1, not "chr1".
##' @param window numeric value saying how far outside the index
##'     credible set we should search. Given that only overlapping
##'     credible sets can show colocalisation there is no reason to
##'     set this particularly large. The default, 1e+5, has been
##'     adequate in all our tests.
##' @param type all, gwas or qtl. Depending on your use call, you may
##'     wish to colocalise against all study types in the Open Targets
##'     dataset, just gwas datasets (you want to colocalise with human
##'     level traits eg to look for pleiotropy) or qtl datasets (you
##'     want to colocalise with molecular datasets eg to nominate
##'     causal genes).
##' @param dir Open Targets credible set files and the study manifest
##'     (`r .ot_manifest_url`) are downloaded as needed and saved in a
##'     cache dir. If dir is NULL, a temporary directory will be
##'     used. It is strongly recommended to use dir to supply a local
##'     directory where files will persist between R sessions. This
##'     will save time and prevent downloading the same files in each
##'     R session. Downloads larger than
##'     `getOption("coloc.max_temp_download", 1e9)` bytes require a
##'     persistent dir.
##'
##'     You can also set this variable persistently using
##'     `options(coloc.cache_dir="path/to/my/directory")`. See
##'     [coloc-options].
##' @param extra extra columns from the Open Targets manifest
##' @param remote FALSE (default): download the credible set files
##'     needed for this query to dir and read them locally. TRUE: read
##'     them directly from the Open Targets S3 bucket without
##'     downloading. Can be set persistently with
##'     `options(coloc.ot_remote = TRUE)`. See [coloc-options].
##' @return a list of credible sets with an attribute "annot"
##'     containing a data.frame where one row per element in the list
##'     describes the annotation of the list (most importantly,
##'     studyId).
##' @export
##' @author Chris Wallace, Jeffrey Pullin
##' @seealso [coloc_query()], [ot_fetch()]
##' @examples
##' \dontrun{
##' q <- data.table(chr = 6,
##'    pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
##'    snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T",
##'            "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
##'    PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))
##'
##' ot <- query_ot(q)
##' out <- coloc_query(q, ot)
##'
##' out[out$PP.H4 > 0 & out$studyType != "gwas",]
##' }
query_ot <- function(cs, window = 1e+5, type = c("all", "gwas", "qtl"), dir = NULL,
    extra=c("traitFromSource","pubmedId"), remote = getOption("coloc.ot_remote", FALSE)) {
  .check_q(cs)
  .check_window(window, ot = TRUE)
  type <- match.arg(type)
  chr <- unique(cs$chr)
  if (length(chr) != 1) stop("`chr` column must only have one value")
  lower_pos <- as.integer(max(0, min(cs$pos) - window))
  upper_pos <- as.integer(max(cs$pos) + window)

  ds <- .ot_dataset(chr, lower_pos, upper_pos, dir, remote)
  if (is.null(ds)) {
    warning("No Open Targets credible set files cover chromosome ", chr)
    return(structure(list(), annot = data.frame()))
  }

  ## Build the filter as an Arrow expression so it's pushed down to the parquet scan
  pos_ref <- arrow::Expression$field_ref("position")
  type_ref <- arrow::Expression$field_ref("studyType")
  filt <- arrow::Expression$field_ref("chromosome") == as.character(chr) & pos_ref >= lower_pos & pos_ref <= upper_pos
  if (type == "gwas") filt <- filt & type_ref == "gwas"
  if (type == "qtl") filt <- filt & type_ref != "gwas"
  cols <- c("studyLocusId", "chromosome", "studyId", "variantId", "position", "region", "studyType", "locus","beta")

  tbl <- .ot_scan(ds, projection = cols, filter = filt)

  credsets <- lapply(tbl$locus$as_vector(), .otlocus_to_cs)
  ok <- vapply(credsets, function(x) is_prob(sum(x$PP)), logical(1))
  if(any(!ok)) {
      warning(sum(!ok), " / ", length(ok), " credsets have sum(PP) not in [0,1]. Removing: ", paste(tbl$studyLocusId$as_vector()[!ok], collapse = ", "))
      credsets <- credsets[ok]
  }
  annot <- data.frame(leadVariant = tbl$variantId$as_vector(),
                      studyLocusId = tbl$studyLocusId$as_vector(),
                      studyId = tbl$studyId$as_vector(),
                      chromosome = tbl$chromosome$as_vector(),
                      position = tbl$position$as_vector(),
                      studyType = tbl$studyType$as_vector())[ok, , drop = FALSE]
  if (nrow(annot) && length(extra)) {
      study <- .ot_study(annot$studyId, dir, cols = extra) # only the rows and columns we need, not the whole manifest
      extra <- intersect(extra, names(study))
      if(length(extra))
          annot <- cbind(annot, study[match(annot$studyId, study$studyId), extra, drop = FALSE])
  }
  rownames(annot) <- NULL

  names(credsets) <- annot$studyLocusId
  attr(credsets, "annot") <- annot
  credsets
}

################################################################################
## eqtl catalogue queries

.eqtl_cat_url <- "https://zenodo.org/records/22846842"
.eqtl_cat_files <- data.frame(quant_method = c("exon", "ge", "leafcutter", "majiq", "tx", "txrev", "metadata"),
                              file = c(sprintf("eQTL_Catalogue_r8-beta_cs_%s_190926.parquet",
                                               c("exon", "ge", "leafcutter", "majiq", "tx", "txrev")),
                                       "dataset_metadata.tsv"),
                              md5 = c("5085e0b6c61b57ff8acd3aad9f7fc297",
                                      "1c54f8dfab4ced90fe0da3154b9b03f7",
                                      "e83fbe0b283c63c28f843db4e155060b",
                                      "3963dbb8ffc347b46d25f13983ef9d0a",
                                      "c11a60ba590e6cf24417de8935cfb699",
                                      "277bf7c2f592a2e48dbfab59618a4e27",
                                      "52b0305f56bd046eb5eeb2b6f2efa874"),
                              size = c(293.8e6, 779.0e6, 228.5e6, 607.3e6, 284.2e6, 325.1e6, 64.6e3)) # approximate bytes, from Zenodo: used for the download size limit, not to check files

eqtl_cat_dir <- function(dir = NULL) .cache_dir(dir, "eqtl_catalogue")

##' Download eQTL Catalogue credible sets
##'
##' Downloads the eQTL Catalogue credible set files for the chosen
##' quantification methods, plus the dataset metadata, from
##' [Zenodo](`r .eqtl_cat_url`) to the cache directory. Files already
##' present are not downloaded again. query_eqtl_cat() calls this
##' itself; call it directly to download in advance, for example on an
##' HPC login node before running query_eqtl_cat() on compute nodes
##' without internet access, using the same dir.
##'
##' @inheritParams query_eqtl_cat
##' @return named character vector of local file paths, named by
##'     quant_method, with the dataset metadata file named "metadata"
##' @export
##' @author Chris Wallace
##' @seealso [query_eqtl_cat()]
eqtl_cat_fetch <- function(quant_method = "ge", dir = NULL) {
  quant_method <- match.arg(quant_method, setdiff(.eqtl_cat_files$quant_method, "metadata"), several.ok = TRUE)
  files <- .eqtl_cat_files[.eqtl_cat_files$quant_method %in% c(quant_method, "metadata"), ]
  local <- eqtl_cat_dir(dir)
  todo <- !file.exists(file.path(local, files$file))
  .check_download_size(sum(files$size[todo]), dir, "This eQTL Catalogue query")
  getone <- function(f, m) .download_once(paste0(.eqtl_cat_url, "/files/", f, "?download=1"), local, file = f, md5 = m)
  stats::setNames(mapply(getone, files$file, files$md5, USE.NAMES = FALSE), files$quant_method)
}

##' query eQTL catalogue credible sets
##'
##' This function takes an input credible set, uses that to define a
##' region of interest, and searches for eQTL Catalogue credible sets
##' that may overlap - ie contain at least one snp within the region
##' +/- window. The resulting list of eQTL Catalogue credible sets may
##' then be colocalised with the index credible set using the
##' coloc_query() function.
##'
##' @inheritParams query_ot
##' @param quant_method the eQTL Catalogue contains results for up to
##'     six different ways of quantifying gene expression: exon, ge,
##'     leafcutter, majiq, tx, txrev. See
##'     <https://www.ebi.ac.uk/eqtl/Methods/>.
##'     [majiq](https://www.biorxiv.org/content/10.1101/2024.07.02.601792v1)
##'     is a recent splicing quantification that was introduced to
##'     eQTL Catalogue in r8. The default is ge, but you can choose
##'     more than one if you wish.
##' @param dir eQTL Catalogue credible sets are available at
##'     [Zenodo](`r .eqtl_cat_url`) subdivided by quant_method. On
##'     first use, this function will download copies of the
##'     appropriate files to a local temporary directory, and query
##'     those files without download on subsequent use.
##'
##'     It is strongly recommended to use dir to supply a local
##'     directory where files will persist between R sessions. This
##'     will save time and prevent downloading the same files in each
##'     R session.
##'
##'     You can also set this variable persistently using
##'     `options(coloc.cache_dir="path/to/my/directory")`. See
##'     [coloc-options].
##'
##'     Downloads larger than `getOption("coloc.max_temp_download",
##'     1e9)` bytes (eg several quant_methods at once) require a
##'     persistent dir.
##' @return a list of credible sets with an attribute "annot"
##'     containing a data.frame where one row per element in the list
##'     describes the annotation of the list (most importantly,
##'     dataset_id).
##' @export
##' @author Chris Wallace
##' @seealso [coloc_query()]
##' @examples
##' \dontrun{
##' q <- data.frame(chr = 6,
##'    pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
##'    snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T",
##'            "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
##'    PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))
##'
##' ecat <- query_eqtl_cat(q)
##' out2 <- coloc_query(q, ecat)
##' out2[out2$PP.H4 > 0, ]
##' }
query_eqtl_cat <- function(cs, quant_method = "ge", dir = NULL, window = 1e+5) {
  .check_q(cs)
  .check_window(window)
  chr <- unique(cs$chr)
  if (length(chr) != 1) stop("`chr` column must only have one value")
  lower_pos <- as.integer(min(cs$pos) - window)
  upper_pos <- as.integer(max(cs$pos) + window)
  paths <- eqtl_cat_fetch(quant_method, dir)
  path_manifest <- paths[["metadata"]]
  ds <- arrow::open_dataset(unname(paths[names(paths) != "metadata"]), format = "parquet")

  pos_ref <- arrow::Expression$field_ref("position")
  ## pos_type <- ds$schema$position$type
  ## lower <- arrow::Scalar$create(lower_pos, type = pos_type)
  ## upper <- arrow::Scalar$create(upper_pos, type = pos_type)
  filt <- arrow::Expression$field_ref("chromosome") == as.character(chr) & pos_ref >= lower_pos & pos_ref <= upper_pos 
  cols <- c("dataset_id", "molecular_trait_id", "molecular_trait_object_id", "cs_id", "region", "variant", "pip","beta")
  ## Create bounds with the column's own type, so arrow doesn't convert positions for the comparison
  
  tbl <- arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable()  |> as.data.frame()
  names(tbl)[names(tbl) == "variant"] <- "snp"
  names(tbl)[names(tbl) == "pip"] <- "PP"
  tbl$snp <- sub("^chr", "", tbl$snp) # match snp format eg 6_90279406_G_C
  tbl <- unique(tbl)
  credsets <- split(tbl[, c("snp", "PP", "beta")], tbl$cs_id) 
  ok <- sapply(credsets, function(x) is_prob(sum(x$PP)))
  if(any(!ok)) {
      warning(sum(!ok), " / ", length(ok), " credsets have sum(PP) not in [0,1]. Removing")
      print(names(credsets)[!ok])
      credsets <- credsets[ok]
  }
  annot <- unique(tbl[, setdiff(names(tbl), c("snp", "PP"))])
  manifest <- data.table::fread(path_manifest, data.table = FALSE)
  annot <- merge(annot, manifest, by = "dataset_id", all.x = TRUE)
  annot <- annot[match(names(credsets), annot$cs_id), ] 
  rownames(annot) <- NULL
  attr(credsets, "annot") <- annot
  credsets
}

################################################################################
## colocalise
.check_snp_ids <- function(snps) {
  bad <- !grepl("^([1-9]|1[0-9]|2[0-2]|X|Y|MT)_[0-9]+_[ACGT]+_[ACGT]+$", snps)
  if (any(bad)) stop(sum(bad), " snp identifiers are not in the format chrom_position_ref_alt (eg \"6_90279406_G_C\"). First few: ", paste(head(snps[bad], 3), collapse = ", "))
}

##' Colocalisation of an index credible set and a list of credible
##' sets generated by query_ot() or query_eqtl_cat()
##'
##' Note that the two credible sets are matched by column "snp".
##' Therefore it is important that the formatting of the snp column in
##' cs matches that in the object returned by query_ot(). It can be
##' helpful to examine that object as shown in the example below.
##'
##' This function uses the coloc_pip method for colocalisation.
##'
##' @param cs index credible set: a data.frame containing columns
##'     "snp" (character) and "PP" (numeric). "snp" will be used to
##'     match against credible sets in querylist and should be of the
##'     format "chr_pos_ref_alt", eg "6_90279406_G_C" (not
##'     "chr6_90279406_G_C").
##' @param querylist output of query_ot() or query_eqtl_cat(). If
##'     querylist is of length 0 (no credible sets), this function
##'     warns and returns NULL.
##' @return a data.frame containing coloc results, one row per
##'     credible set in querylist, with the annotation from querylist
##' @export
##' @author Chris Wallace
##' @seealso [query_ot()], [query_eqtl_cat()], [coloc_pip()]
##' @examples
##' \dontrun{
##' q <- data.frame(chr = 6,
##'    pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
##'    snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T",
##'            "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
##'    PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))
##'
##' ot <- query_ot(q)
##' head(ot[[1]]) # check the snp format matches that in q
##' head(q)
##' out <- coloc_query(q, ot)
##' out[out$PP.H4 > 0 & out$studyType != "gwas", ]
##' }
coloc_query <- function(cs, querylist) {
    stopifnot("cs must be a data.frame containing character column snp and numeric column PP"=
                  is.data.frame(cs) &&
                  all(c("snp","PP") %in% colnames(cs)) &&
                  is.character(cs$snp) &&
                  is.numeric(cs$PP))
    ## check expected pattern for SNP identifier - this is a common gotcha
    .check_snp_ids(cs$snp)
    if(!length(querylist)) {
        warning("querylist contains no credible sets.")
        return(NULL)
    }
    stopifnot("querylist must be an object returned by query_ot()"=
                  is.list(querylist) &&
                  is.data.frame(querylist[[1]]) &&
                  "annot" %in% names(attributes(querylist)))
    coloc_res <- lapply(querylist, function(x) coloc_pip(cs, x)) |>
        do.call("rbind", args=_) |>
        as.data.frame()
    cbind(attr(querylist, "annot"), coloc_res)
}

################################################################################

## ensembl to symbol
.hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
.gene_env <- new.env(parent = emptyenv())

## Ensembl id -> symbol table from HGNC, downloaded once and read once per session
.hgnc_genes <- function(dir = NULL) {
  if (is.null(.gene_env$hgnc)) {
    path <- .download_once(.hgnc_url, eqtl_cat_dir(dir))
    g <- data.table::fread(path, select = c("symbol", "ensembl_gene_id"), data.table = FALSE)
    .gene_env$hgnc <- g[g$ensembl_gene_id != "", ]
  }
  .gene_env$hgnc
}

##' Translate Ensembl gene ids to HGNC gene symbols
##'
##' Uses the HGNC complete gene set, downloaded on first use to the same
##' directory as the eQTL Catalogue files. Version suffixes (eg
##' ENSG00000160712.13) are ignored.
##'
##' @param ids character vector of Ensembl gene ids
##' @inheritParams query_eqtl_cat
##' @return character vector of gene symbols, the same length as ids, with NA where no symbol is found
##' @export
ensembl_to_symbol <- function(ids, dir = NULL) {
    g <- .hgnc_genes(dir)
    ids <- sub("\\..*","",ids)
  g$symbol[match(sub("\\.[0-9]+$", "", ids), g$ensembl_gene_id)]
}
