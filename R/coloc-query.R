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
.download_once <- function(url, dir, file = basename(url), md5 = NULL) {
  path <- file.path(dir, file)
  if(file.exists(path))
      return(path)
  old <- options(timeout = max(3600, getOption("timeout"))) # default 60s is far too short for large files
  on.exit(options(old))
  message("Downloading ", file, " to ", dir)
  tmp <- paste0(path, ".part") # an interrupted download never looks like a complete file
  utils::download.file(url, tmp, mode = "wb")
  if(!is.null(md5) && unname(tools::md5sum(tmp)) != md5)
    warning("Checksum mismatch for ", file, "; download may be corrupt. Delete ", path, " to force a fresh download.")
  file.rename(tmp, path)
  path
}

################################################################################
## Open Targets queries

.ot_release <- "26.06"
.ot_credset_url <- sprintf("s3://open-targets-public-data-releases/platform/%s/output/credible_set", .ot_release)
.ot_manifest_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/output/study/part-00000-1c0a4198-4d65-425e-82a1-d9d8b54d8549-c000.snappy.parquet", .ot_release)

## Study columns added to the query_ot() annotation, where present in the study file
## See https://platform.opentargets.org/downloads/study/schema
.ot_study_cols <- c("studyId", "projectId", "traitFromSource", "geneId", "biosampleFromSourceId", "pubmedId")

## Session cache for opened datasets
.ot_env <- new.env(parent = emptyenv())

.ot_dir <- function(dir = NULL) .cache_dir(dir, "opentargets")

## Remote credible set dataset, opened once per session
.ot_dataset <- function() {
  if (is.null(.ot_env$ds)) .ot_env$ds <- arrow::open_dataset(.ot_credset_url, format = "parquet")
  .ot_env$ds
}

## OT manifest
.ot_manifest <- function(dir = NULL) {
    path <- .download_once(.ot_manifest_url, .ot_dir(dir))
    if (is.null(.ot_env$study))
        .ot_env$study <- arrow::open_dataset(path, format = "parquet")
  .ot_env$study
}

## Study annotations for a set of studyIds
.ot_study <- function(study_ids, dir = NULL) {
  ds <- .ot_manifest(dir)
  cols <- intersect(.ot_study_cols, names(ds))
  ids <- arrow::Array$create(unique(study_ids))
  filt <- arrow::Expression$create("is_in", arrow::Expression$field_ref("studyId"), options = list(value_set = ids, skip_nulls = TRUE))
  as.data.frame(arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable())
}

## Convert one Open Targets locus (struct of variants) to a credible set
.otlocus_to_cs <- function(locus, ...) {
  parts <- strsplit(locus$variantId, "_", fixed = TRUE) # split once, reuse
  data.table(...,
             chr = as.integer(vapply(parts, `[[`, character(1), 1)),
             pos = as.integer(vapply(parts, `[[`, character(1), 2)),
             snp = locus$variantId,
             PP = locus$posteriorProbability)
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
##' The Open Targets dataset is on an s3 connection. To make query_ot()
##' time efficient, we open that dataset on the first use, and cache
##' it. Subsequent calls should be quicker because they reuse this
##' cached arrow::open_dataset() call. However, this may fail if used
##' inside parallel::mclapply or other forked workers. Use
##' socket-based workers (parLapply, future's multisession) and let
##' each worker open its own dataset.
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
##' @param dir The manifest for Open Target files is downloaded from
##'     `r .ot_manifest_url` and saved in a cache dir. If dir is NULL,
##'     a temporary directory will be used. It is strongly recommended
##'     to use dir to supply a local directory where files will
##'     persist between R sessions. This will save time and prevent
##'     downloading the same files in each R session.
##'
##'     You can also set this variable persistently using
##'     `options(coloc.cache_dir="path/to/my/directory")`
##' @param extra extra columns from the Open Targets manifest
##' @return a list of credible sets with an attribute "annot"
##'     containing a data.frame where one row per element in the list
##'     describes the annotation of the list (most importantly,
##'     studyId).
##' @export
##' @author Chris Wallace, Jeffrey Pullin
##' @seealso [coloc_query()]
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
    extra=c("traitFromSource","pubmedID")) {
  .check_q(cs)
  .check_window(window, ot = TRUE)
  type <- match.arg(type)
  chr <- unique(cs$chr)
  if (length(chr) != 1) stop("`chr` column must only have one value")
  lower_pos <- min(cs$pos) - window
  upper_pos <- max(cs$pos) + window

  ## Build the filter as an Arrow expression so it's pushed down to the parquet scan
  pos_ref <- arrow::Expression$field_ref("position")
  type_ref <- arrow::Expression$field_ref("studyType")
  filt <- arrow::Expression$field_ref("chromosome") == as.character(chr) & pos_ref >= lower_pos & pos_ref <= upper_pos
  if (type == "gwas") filt <- filt & type_ref == "gwas"
  if (type == "qtl") filt <- filt & type_ref != "gwas"
  cols <- c("studyLocusId", "chromosome", "studyId", "variantId", "position", "region", "studyType", "locus")

  tbl <- arrow::Scanner$create(.ot_dataset(), projection = cols, filter = filt)$ToTable()

  credsets <- lapply(tbl$locus$as_vector(), .otlocus_to_cs)
  annot <- data.frame(leadVariant = tbl$variantId$as_vector(),
                      studyLocusId = tbl$studyLocusId$as_vector(),
                      studyId = tbl$studyId$as_vector(),
                      chromosome = tbl$chromosome$as_vector(),
                      position = tbl$position$as_vector(),
                      studyType = tbl$studyType$as_vector())

  if (nrow(annot) & length(extra)) {
      manifest <- .ot_manifest(dir) |> as.data.frame() |> as.data.frame() # need to do this twice to avoid tibbles
      extra <- intersect(extra, names(manifest))
      if(length(extra))
          annot <- cbind(annot, manifest[match(annot$studyId, manifest$studyId), extra, drop = FALSE])
      rownames(annot) <- NULL
  }

  names(credsets) <- annot$studyLocusId
  attr(credsets, "annot") <- annot
  credsets
}

################################################################################
## eqtl catalogue queries

.eqtl_cat_url <- "https://zenodo.org/records/22846842"
.eqtl_cat_files <- data.frame(quant_method = c("exon", "ge", "leafcutter", "majiq", "tx", "txrev",NA),
                              file = c(sprintf("eQTL_Catalogue_r8-beta_cs_%s_190926.parquet",
                                               c("exon", "ge", "leafcutter", "majiq", "tx", "txrev")),
                                       "dataset_metadata.tsv"),
                              md5 = c("5085e0b6c61b57ff8acd3aad9f7fc297",
                                      "1c54f8dfab4ced90fe0da3154b9b03f7",
                                      "e83fbe0b283c63c28f843db4e155060b",
                                      "3963dbb8ffc347b46d25f13983ef9d0a",
                                      "c11a60ba590e6cf24417de8935cfb699",
                                      "277bf7c2f592a2e48dbfab59618a4e27",
                                      "52b0305f56bd046eb5eeb2b6f2efa874"))

eqtl_cat_dir <- function(dir = NULL) .cache_dir(dir, "eqtl_catalogue")

eqtl_cat_fetch <- function(quant_method = "ge", dir = NULL) {
  quant_method <- match.arg(quant_method, setdiff(.eqtl_cat_files$quant_method, "metadata"), several.ok = TRUE)
  dir <- eqtl_cat_dir(dir)
  files <- .eqtl_cat_files[.eqtl_cat_files$quant_method %in% c(quant_method) | is.na(.eqtl_cat_files$quant_method), ]
  getone <- function(f, m) .download_once(paste0(.eqtl_cat_url, "/files/", f, "?download=1"), dir, file = f, md5 = m)
  mapply(getone, files$file, files$md5)
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
##'     `options(coloc.cache_dir="path/to/my/directory")`
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
  lower_pos <- min(cs$pos) - window
  upper_pos <- max(cs$pos) + window
  paths <- eqtl_cat_fetch(quant_method, dir) # last path is for the manifest
  path_manifest <- paths[ length(paths) ]
  paths <- paths[ -length(paths) ]
  ds <- arrow::open_dataset(paths, format = "parquet")

  pos_ref <- arrow::Expression$field_ref("position")
  filt <- arrow::Expression$field_ref("chromosome") == as.character(chr) & pos_ref >= lower_pos & pos_ref <= upper_pos
  cols <- c("dataset_id", "molecular_trait_id", "molecular_trait_object_id", "cs_id", "region", "variant", "pip")

  tbl <- arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable()  |> as.data.frame()
  names(tbl)[names(tbl) == "variant"] <- "snp"
  names(tbl)[names(tbl) == "pip"] <- "PP"
  tbl$snp <- sub("^chr", "", tbl$snp) # match snp format eg 6_90279406_G_C
  credsets <- split(tbl[, c("snp", "PP")], tbl$cs_id)
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
    stopifnot("querylist must be an object returned by query_ot()"=
                  is.list(querylist) &&
                  is.data.frame(querylist[[1]]) &&
                  "annot" %in% names(attributes(querylist)))
    if(!length(querylist)) {
        warning("querylist contains no credible sets.")
        return(NULL)
    }
    coloc_res <- lapply(querylist, function(x) coloc_pip(cs, x)) |>
        do.call("rbind", args=_) |>
        as.data.frame()
    cbind(attr(querylist, "annot"), coloc_res)
}
