##' @import arrow
NULL

## Shared helpers

## Check an index credible set for the ot_query and query_ecat functions.
## return its single chromosome as used in the public data: "6", "X".
## chr may be numeric or character
.check_q <- function(cs) {
    stopifnot("cs must be a data.frame containing columns chr and pos, with pos numeric"=
                  is.data.frame(cs) &&
                  all(c("chr","pos") %in% names(cs)) &&
                  is.numeric(cs$pos))
    stopifnot("cs must have at least one non-missing pos" = any(!is.na(cs$pos)))
    stopifnot("chr must be numeric or character, eg 6, \"6\" or \"X\""=
                  is.numeric(cs$chr) || is.character(cs$chr) || is.factor(cs$chr))
    chr <- unique(sub("^chr", "", as.character(cs$chr), ignore.case = TRUE)) # tolerate "chr6", though not documented
    if (length(chr) != 1) stop("`chr` column must only have one value")
    chr
}

.check_window <- function(window, ot = FALSE) {
    stopifnot("window should be a single non-negative number" = is.numeric(window) && length(window) == 1 && window >= 0)
    if (ot && window > 5e+5)
        warning("You have set window very large. This is probably unnecessary and will cause downloads to be slower.")
}

## Resolve a data directory: the coloc.cache_dir option, else a per-session temp folder
.cache_dir <- function(subdir) {
    dir <- getOption("coloc.cache_dir")
    if(is.null(dir)) {
        warning("Downloading files to a temporary directory. Set options(coloc.cache_dir = \"path/to/dir\") to keep downloads between R sessions. This will save you time by preventing unnecessary downloads in each new R session.")
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
.ot_s3_region <- "eu-west-1" # bucket region, for the https url; check with: curl -sI https://open-targets-public-data-releases.s3.amazonaws.com | grep -i x-amz-bucket-region
.ot_s3_bucket <- "open-targets-public-data-releases"
.ot_s3_url <- sprintf("https://%s.s3.%s.amazonaws.com", .ot_s3_bucket, .ot_s3_region) # plain https access to the public bucket
.ot_credset_prefix <- sprintf("platform/%s/output/credible_set/", .ot_release)
.ot_manifest_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/%s/output/study/part-00000-1c0a4198-4d65-425e-82a1-d9d8b54d8549-c000.snappy.parquet", .ot_release)

## Study columns added to the query_ot() annotation, where present in the study file
## See https://platform.opentargets.org/downloads/study/schema
.ot_study_cols <- c("studyId", "projectId", "traitFromSource", "geneId", "biosampleFromSourceId", "pubmedId")

## Session cache for small lookup tables (file index, study manifest)
.ot_env <- new.env(parent = emptyenv())

.ot_dir <- function() .cache_dir("opentargets")
.ot_local_dir <- function() file.path(.ot_dir(), paste0("credible_set_", .ot_release))

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
.check_download_size <- function(bytes, what) {
    limit <- getOption("coloc.max_temp_download", 1e9)
    if (is.null(getOption("coloc.cache_dir")) && bytes > limit)
        stop(what, " needs a ", round(bytes / 1e9, 1), " GB download, which is too large for a temporary directory. ",
             "Please set a persistent cache directory with options(coloc.cache_dir = \"path/to/dir\"), so it is only downloaded once.",
             call. = FALSE)
}

## Make sure the given credible set files are in the local cache, downloading any that are missing; returns local paths
.ot_fetch_files <- function(files, what = "This Open Targets query") {
    local <- .ot_local_dir()
    paths <- file.path(local, files$file)
    todo <- which(!file.exists(paths) | file.size(paths) != files$size)
    if (length(todo)) {
        .check_download_size(sum(files$size[todo]), what)
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
##' cache directory.
##'
##' Because of its size, a persistent cache directory
##' (`options(coloc.cache_dir=...)`, see [coloc-options]) is required.
##' Files already present
##' with the expected size are not downloaded again, so an interrupted
##' download can be resumed by calling ot_fetch() again.
##'
##' @return the path to the local credible set directory, invisibly
##' @export
##' @author Chris Wallace
##' @seealso [query_ot()]
ot_fetch <- function() {
    .ot_fetch_files(.ot_index(), what = "Downloading all Open Targets credible sets")
    invisible(.ot_local_dir())
}

## Credible set dataset for chr:lower-upper: the files overlapping it, downloaded to the cache dir if not already there
.ot_dataset <- function(chr, lower, upper) {
    files <- .ot_files_for(chr, lower, upper)
    if (!nrow(files)) return(NULL)
    arrow::open_dataset(.ot_fetch_files(files), format = "parquet")
}

## OT manifest
.ot_manifest <- function() {
    path <- .download_once(.ot_manifest_url, .ot_dir())
    if (is.null(.ot_env$study))
        .ot_env$study <- arrow::open_dataset(path, format = "parquet")
    .ot_env$study
}

## Study annotations for a set of studyIds
.ot_study <- function(study_ids, cols = .ot_study_cols) {
    ds <- .ot_manifest()
    cols <- intersect(unique(c("studyId", cols)), names(ds))
    ids <- arrow::Array$create(unique(study_ids))
    filt <- arrow::Expression$create("is_in", arrow::Expression$field_ref("studyId"), options = list(value_set = ids, skip_nulls = TRUE))
    as.data.frame(arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable())
}

## Convert one Open Targets locus (struct of variants) to a credible set
.otlocus_to_cs <- function(locus, ...) {
    parts <- strsplit(locus$variantId, "_", fixed = TRUE) # split once, reuse
    ## ids are chr_pos_ref_alt, except variants on alt/patch contigs,
    ## eg 17_KI270857v1_alt_811563_C_G. Their position is on the
    ## contig, not the chromosome, so it is set to NA
    primary <- lengths(parts) == 4
    pos <- rep(NA_integer_, length(parts))
    pos[primary] <- as.integer(vapply(parts[primary], `[[`, character(1), 2))
    data.table(...,
               chr = vapply(parts, `[[`, character(1), 1),
               pos = pos,
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
##' queries. To download everything up front (eg before working on
##' HPC compute nodes without internet access), use [ot_fetch()].
##'
##' Files are kept in the directory given by
##' `options(coloc.cache_dir = "path/to/dir")`, so they persist between
##' R sessions. Setting this is strongly recommended: otherwise a
##' temporary directory is used, and files are downloaded again in
##' every session. See [coloc-options].
##'
##' A few Open Targets credible sets include variants on alt or patch
##' contigs (eg 17_KI270857v1_alt_811563_C_G). These are returned as
##' Open Targets provides them, with pos set to NA. Because they cannot
##' match variants on the primary chromosomes, their posterior
##' probability counts against colocalisation, so query_ot() warns
##' when it finds them: assess any lack of colocalisation with those
##' credible sets carefully.
##'
##' @param cs an index data.frame containing columns "chr" and
##'     "pos", for a single chromosome. chr may be numeric or
##'     character, eg 6, "6" or "X". Credible sets returned by
##'     query_ot() can be used directly.
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
##' @param extra extra columns from the Open Targets manifest
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
query_ot <- function(cs, window = 1e+5, type = c("all", "gwas", "qtl"),
                     extra=c("traitFromSource","pubmedId")) {
    chr <- .check_q(cs)
    .check_window(window, ot = TRUE)
    type <- match.arg(type)
    lower_pos <- as.integer(max(0, min(cs$pos, na.rm = TRUE) - window)) # NA pos: eg variants on alt contigs
    upper_pos <- as.integer(max(cs$pos, na.rm = TRUE) + window)

    ds <- .ot_dataset(chr, lower_pos, upper_pos)
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

    tbl <- arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable()

    credsets <- lapply(tbl$locus$as_vector(), .otlocus_to_cs)
    ok <- vapply(credsets, function(x) is_prob(sum(x$PP)), logical(1))
    if(any(!ok)) {
        warning(sum(!ok), " / ", length(ok), " credsets have sum(PP) not in [0,1]. Removing: ", paste(head(tbl$studyLocusId$as_vector()[!ok]), collapse = ", "))
        credsets <- credsets[ok]
    }
    annot <- data.frame(leadVariant = tbl$variantId$as_vector(),
                        studyLocusId = tbl$studyLocusId$as_vector(),
                        studyId = tbl$studyId$as_vector(),
                        chromosome = tbl$chromosome$as_vector(),
                        position = tbl$position$as_vector(),
                        studyType = tbl$studyType$as_vector())[ok, , drop = FALSE]
    if (nrow(annot) && length(extra)) {
        study <- .ot_study(annot$studyId, cols = extra) # only the rows and columns we need, not the whole manifest
        extra <- intersect(extra, names(study))
        if(length(extra))
            annot <- cbind(annot, study[match(annot$studyId, study$studyId), extra, drop = FALSE])
    }
    rownames(annot) <- NULL

    names(credsets) <- annot$studyLocusId
    attr(credsets, "annot") <- annot

    ## variants on alt/patch contigs (eg 17_KI270857v1_alt_811563_C_G) can't match variants in the index credible set,
    ## so their PP counts against colocalisation. Data are left as Open Targets provides them, but flagged
    alt_snps <- lapply(credsets, function(x) grep("_(alt|fix|random)_", x$snp, value = TRUE))
    has_alt <- lengths(alt_snps) > 0
    if (any(has_alt))
        warning(sum(has_alt), " Open Targets credible set(s) contain variants on alt/patch contigs. These cannot match variants on the primary chromosomes. So colocalisation with these credible sets may be underestimated and H3 overstated: assess any lack of colocalisation with them carefully. Affected studyLocusId(s): ",
                paste(utils::head(names(credsets)[has_alt]), collapse = ", "),
                if (sum(has_alt) > 5) ", ...", call. = FALSE)
    credsets
}

################################################################################
## eqtl catalogue queries

## Releases: eQTL Catalogue credible sets split into small files sorted by genomic position
## (data-raw/eqtl_cat_split.R), published on Zenodo with one record per release (data-raw/zenodo_upload.sh).
## To add a release: split, upload, copy its index to inst/extdata, and add a row here.
.ecat_releases <- data.frame(
    version = c("r7", "r8-beta"),
    record = c("23397994", "23398345"),
    source_record = c("23350544", "22846842"), # the original eQTL Catalogue records the files were split from
    metadata = c("dataset_metadata_r7.tsv", "dataset_metadata.tsv"),
    metadata_md5 = c("7020d45bc577bff33f2904671a97eb03", "52b0305f56bd046eb5eeb2b6f2efa874"))
.ecat_default_version <- "r7"

.ecat_release <- function(version) {
    if (length(version) != 1 || !version %in% .ecat_releases$version)
        stop("version must be one of: ", paste(.ecat_releases$version, collapse = ", "), call. = FALSE)
    .ecat_releases[.ecat_releases$version == version, ]
}
.ecat_url <- function(version = .ecat_default_version) sprintf("https://zenodo.org/records/%s", .ecat_release(version)$record)

## Session cache for the file indexes
.ecat_env <- new.env(parent = emptyenv())

## Index of a release's chunk files, one row per chunk and chromosome:
## quant_method, file, chromosome, pos_min, pos_max, n_cs, n_rows, size, md5
.ecat_index <- function(version = .ecat_default_version) {
    .ecat_release(version)
    if (is.null(.ecat_env[[version]])) {
        f <- system.file("extdata", sprintf("eQTL_Catalogue_index_%s.tsv", version), package = "coloc")
        if (!nzchar(f)) stop("eQTL Catalogue file index for ", version, " not found in inst/extdata")
        .ecat_env[[version]] <- utils::read.delim(f, colClasses = c(chromosome = "character", md5 = "character"))
    }
    .ecat_env[[version]]
}

## Check quant_method against what the release contains (eg r7 has no majiq)
.ecat_check_qm <- function(quant_method, version) {
    qms <- unique(.ecat_index(version)$quant_method)
    bad <- setdiff(quant_method, qms)
    if (length(bad))
        stop("quant_method ", paste(bad, collapse = ", "), " not available in eQTL Catalogue ", version, ". Choose from: ", paste(qms, collapse = ", "), call. = FALSE)
    quant_method
}

## Cache directory for one release, so releases never overwrite each other
ecat_dir <- function(version = .ecat_default_version)
    file.path(.cache_dir("eqtl_catalogue"), paste0("eQTL_Catalogue_", .ecat_release(version)$version))

## Make sure the given chunk files (rows of the index) and the release metadata are in the cache, downloading any
## that are missing. Returns local paths named by quant_method, with the metadata file named "metadata"
.ecat_fetch_files <- function(files, version, what) {
    rel <- .ecat_release(version)
    local <- ecat_dir(version)
    base <- paste0(.ecat_url(version), "/files/")
    files <- unique(files[, c("quant_method", "file", "size", "md5")])
    have <- file.exists(file.path(local, files$file)) & file.size(file.path(local, files$file)) == files$size
    .check_download_size(sum(files$size[!have]), what)
    if (any(!have))
        message("Downloading ", sum(!have), " eQTL Catalogue ", version, " file(s) (", round(sum(files$size[!have]) / 1e6), " MB) to ", local)
    getone <- function(f, m, s) .download_once(paste0(base, f, "?download=1"), local, file = f, md5 = m, size = s)
    paths <- stats::setNames(mapply(getone, files$file, files$md5, files$size, USE.NAMES = FALSE), files$quant_method)
    ## c(paths, metadata = .download_once(paste0(base, rel$metadata, "?download=1"), local, file = rel$metadata, md5 = rel$metadata_md5))
   c(paths, metadata = .ecat_metadata_path(version)) 
}

##' Download eQTL Catalogue credible sets
##'
##' Downloads all the eQTL Catalogue credible set files for the chosen
##' quantification methods and release, plus the dataset metadata, to
##' the cache directory. query_ecat() downloads only the files it
##' needs for each query, so you only need this to download
##' everything in advance, for example on an HPC login node before
##' running query_ecat() on compute nodes without internet access,
##' using the same cache directory. Files already present are not
##' downloaded again.
##'
##' @inheritParams query_ecat
##' @return named character vector of local file paths, named by
##'     quant_method, with the dataset metadata file named "metadata"
##' @export
##' @author Chris Wallace
##' @seealso [query_ecat()]
ecat_fetch <- function(quant_method = "ge", version = .ecat_default_version) {
    quant_method <- .ecat_check_qm(quant_method, version)
    idx <- .ecat_index(version)
    .ecat_fetch_files(idx[idx$quant_method %in% quant_method, ], version, "Downloading these eQTL Catalogue files")
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
##' The eQTL Catalogue credible sets are stored, for each release and
##' quantification method, as a series of files each covering a block
##' of the genome. query_ecat() downloads only the file(s)
##' overlapping the query region (typically one, ~25 MB per
##' quantification method) to the cache dir, and reuses them in later
##' queries. To download everything in advance, use
##' [ecat_fetch()]. As for [query_ot()], files are kept in the
##' directory given by `options(coloc.cache_dir = ...)`; see
##' [coloc-options]. These files are the eQTL Catalogue's own
##' credible sets, unchanged but split by genomic position: release r7
##' at <`r .ecat_url("r7")`> and r8-beta at
##' <`r .ecat_url("r8-beta")`>.
##'
##' @inheritParams query_ot
##' @param quant_method the eQTL Catalogue contains results for up to
##'     six different ways of quantifying gene expression: exon, ge,
##'     leafcutter, majiq, tx, txrev. See
##'     <https://www.ebi.ac.uk/eqtl/Methods/>.
##'     [majiq](https://www.biorxiv.org/content/10.1101/2024.07.02.601792v1)
##'     is a recent splicing quantification that was introduced to
##'     eQTL Catalogue in r8, so is not available in r7. The default
##'     is ge, but you can choose more than one if you wish.
##' @param version eQTL Catalogue release: "r7" (default,
##'     <https://zenodo.org/records/23397994>) or "r8-beta"
##'     (<https://zenodo.org/records/23398345>). Each release is kept
##'     in its own subdirectory of the cache directory.
##' @return a list of credible sets with an attribute "annot"
##'     containing a data.frame where one row per element in the list
##'     describes the annotation of the list (most importantly,
##'     dataset_id).
##' @references If you use eQTL Catalogue results, please cite: Kerimov N, Hayhurst JD, et al. (2021) A
##'     compendium of uniformly processed human gene expression and splicing quantitative trait loci. Nature
##'     Genetics 53:1290-1299. \doi{10.1038/s41588-021-00924-w}
##' @export
##' @author Chris Wallace
##' @seealso [coloc_query()], [ecat_fetch()]
##' @examples
##' \dontrun{
##' q <- data.frame(chr = 6,
##'    pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
##'    snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T",
##'            "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
##'    PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))
##'
##' ecat <- query_ecat(q)
##' out2 <- coloc_query(q, ecat)
##' out2[out2$PP.H4 > 0, ]
##'
##' ecat8 <- query_ecat(q, version = "r8-beta") # the newer, beta release
##' }
query_ecat <- function(cs, quant_method = "ge", window = 1e+5, version = .ecat_default_version) {
    chr <- .check_q(cs)
    .check_window(window)
    lower_pos <- as.integer(max(0, min(cs$pos, na.rm = TRUE) - window)) # NA pos: eg variants on alt contigs
    upper_pos <- as.integer(max(cs$pos, na.rm = TRUE) + window)
    quant_method <- .ecat_check_qm(quant_method, version)

    ## only the chunk files overlapping the query region
    idx <- .ecat_index(version)
    files <- idx[idx$quant_method %in% quant_method & idx$chromosome == as.character(chr) & idx$pos_max >= lower_pos & idx$pos_min <= upper_pos, ]
    if (!nrow(files)) {
        warning("No eQTL Catalogue ", version, " credible sets for ", paste(quant_method, collapse = ", "), " near chromosome ", chr, ":", lower_pos, "-", upper_pos)
        return(structure(list(), annot = data.frame()))
    }
    paths <- .ecat_fetch_files(files, version, "This eQTL Catalogue query")
    path_manifest <- paths[["metadata"]]
    ds <- arrow::open_dataset(unname(paths[names(paths) != "metadata"]), format = "parquet")

    pos_ref <- arrow::Expression$field_ref("position")
    filt <- arrow::Expression$field_ref("chromosome") == as.character(chr) & pos_ref >= lower_pos & pos_ref <= upper_pos
    ## releases differ: r8-beta has molecular_trait_object_id, r7 has gene_id instead; take whichever are present
    cols <- intersect(c("dataset_id", "molecular_trait_id", "molecular_trait_object_id", "gene_id", "cs_id", "region", "variant", "pip","beta"), names(ds))

    tbl <- arrow::Scanner$create(ds, projection = cols, filter = filt)$ToTable()  |> as.data.frame()
    names(tbl)[names(tbl) == "variant"] <- "snp"
    names(tbl)[names(tbl) == "pip"] <- "PP"
    tbl$snp <- sub("^chr", "", tbl$snp) # match snp format eg 6_90279406_G_C
    tbl <- unique(tbl)
    credsets <- split(tbl[, c("snp", "PP", "beta")], tbl$cs_id)
    ok <- sapply(credsets, function(x) is_prob(sum(x$PP)))
    if(any(!ok)) {
        warning(sum(!ok), " / ", length(ok), " credsets have sum(PP) not in [0,1]. Removing: ", paste(head(names(credsets)[!ok]), collapse = ", "))
        credsets <- credsets[ok]
    }
    annot <- unique(tbl[, setdiff(names(tbl), c("snp", "PP", "beta"))]) # one row per credible set: variant-level columns stay in credsets
    ## manifest <- data.table::fread(path_manifest, data.table = FALSE)
    manifest <- ecat_metadata(version) 
    annot <- merge(annot, manifest, by = "dataset_id", all.x = TRUE)
    annot <- annot[match(names(credsets), annot$cs_id), ]
    rownames(annot) <- NULL
    attr(credsets, "annot") <- annot
    credsets
}

################################################################################
## colocalise
.check_snp_ids <- function(snps) {
    bad <- !grepl("^([1-9]|1[0-9]|2[0-2]|X|Y|MT)(_[A-Za-z0-9.]+_(alt|fix|random))?_[0-9]+_[ACGT]+_[ACGT]+$", snps) # optional alt/patch contig, as in Open Targets
    if (any(bad)) stop(sum(bad), " snp identifiers are not in the format chrom_position_ref_alt (eg \"6_90279406_G_C\"). First few: ", paste(head(snps[bad], 3), collapse = ", "))
}

##' Colocalisation of an index credible set and a list of credible
##' sets generated by query_ot() or query_ecat()
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
##' @param querylist output of query_ot() or query_ecat(). If
##'     querylist is of length 0 (no credible sets), this function
##'     warns and returns NULL.
##' @return a data.table containing coloc results, one row per
##'     credible set in querylist, with the annotation from querylist
##' @export
##' @author Chris Wallace
##' @seealso [query_ot()], [query_ecat()], [coloc_pip()]
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
    as.data.table(cbind(attr(querylist, "annot"), coloc_res))
}

################################################################################
## many credible sets at once

## Run f(cs, name) for each credible set in cs_list, for the *_many() functions.
## - a credible set whose call fails is skipped (its result is NULL) and its error recorded
## - max_fail failures in a row stop the run, as do failures of every credible set, as that suggests a systematic
##   problem (eg no network) rather than bad input
## - warnings are collected and reported once each at the end (or on stopping), with the credible sets they came from
## Returns a list of results named as cs_list, with attribute "failed": the error message for each skipped credible set
.run_many <- function(cs_list, f, what, max_fail) {
    if (is.data.frame(cs_list)) stop("cs_list should be a list of credible sets. For a single credible set, use the function without _many", call. = FALSE)
    stopifnot("cs_list must be a non-empty list" = is.list(cs_list) && length(cs_list) > 0,
              "max_fail must be a positive number" = is.numeric(max_fail) && length(max_fail) == 1 && max_fail >= 1)
    if (is.null(names(cs_list))) names(cs_list) <- seq_along(cs_list)
    if (anyDuplicated(names(cs_list))) stop("names of cs_list must be unique", call. = FALSE)
    warns <- character(); warn_from <- character(); errs <- character(); in_a_row <- 0
    on.exit({ # report, whether we finish or stop
        for (w in unique(warns))
            warning(w, " [", what, ", credible set(s): ", paste(unique(warn_from[warns == w]), collapse = ", "), "]", call. = FALSE)
        if (length(errs))
            warning(length(errs), " of ", length(cs_list), " credible set(s) failed in ", what, " and were skipped: ",
                    paste0(names(errs), ": ", errs, collapse = "; "), call. = FALSE)
    })
    res <- stats::setNames(vector("list", length(cs_list)), names(cs_list))
    for (nm in names(cs_list)) {
        out <- withCallingHandlers(tryCatch(f(cs_list[[nm]], nm), error = function(e) e),
                                   warning = function(w) {
                                       warns <<- c(warns, conditionMessage(w)); warn_from <<- c(warn_from, nm)
                                       invokeRestart("muffleWarning")
                                   })
        if (inherits(out, "error")) {
            errs[nm] <- conditionMessage(out); in_a_row <- in_a_row + 1
            if (in_a_row >= max_fail)
                stop(in_a_row, " credible sets in a row failed in ", what, ", so stopping, as this suggests a systematic problem (eg no network). Last error, for credible set ", nm, ": ", errs[nm], call. = FALSE)
        } else {
            in_a_row <- 0
            res[nm] <- list(out)
        }
    }
    if (length(errs) == length(cs_list))
        stop("all ", length(cs_list), " credible sets failed in ", what, ". Last error: ", errs[length(errs)], call. = FALSE)
    attr(res, "failed") <- errs
    res
}

##' Query Open Targets or eQTL Catalogue for many credible sets
##'
##' Runs [query_ot()] or [query_ecat()] for each credible set in a
##' list. The results can then be colocalised with the same credible
##' sets using [coloc_query_many()]. Keeping the two steps separate
##' means the queries, which need network access, can be run once
##' (eg on an HPC login node) and saved with saveRDS(), and the
##' colocalisation run later or repeated without downloading again.
##'
##' A credible set whose query fails is skipped, with its error
##' recorded, and the others continue. If `max_fail` credible sets in
##' a row fail, or all of them do, the run stops, as that suggests a
##' systematic problem such as no network access. Warnings are
##' collected and reported once each at the end, with the names of the
##' credible sets they came from.
##'
##' @param cs_list a list of index credible sets, each as described
##'     for [query_ot()]. Names are used to label the results; if
##'     unnamed, the credible sets are numbered.
##' @param ... further arguments passed to [query_ot()] or
##'     [query_ecat()], eg type, quant_method, version, window
##' @param max_fail stop after this many credible sets in a row fail
##' @return a list, named as cs_list, of the query results for each
##'     credible set (NULL for those that failed). Attribute "failed"
##'     gives the error message for each failed credible set.
##' @export
##' @author Chris Wallace
##' @seealso [coloc_query_many()]
##' @examples
##' \dontrun{
##' cs_list <- list(BACH2 = q, other = q2) # index credible sets, as for query_ot()
##' ot <- query_ot_many(cs_list, type = "qtl")
##' ecat <- query_ecat_many(cs_list, quant_method = c("ge", "leafcutter"))
##' res_ot <- coloc_query_many(cs_list, ot)
##' res_ecat <- coloc_query_many(cs_list, ecat)
##' }
query_ot_many <- function(cs_list, ..., max_fail = 5)
    .run_many(cs_list, function(cs, nm) query_ot(cs, ...), "query_ot", max_fail)

##' @rdname query_ot_many
##' @export
query_ecat_many <- function(cs_list, ..., max_fail = 5)
    .run_many(cs_list, function(cs, nm) query_ecat(cs, ...), "query_ecat", max_fail)

##' Colocalise many credible sets with their query results
##'
##' Runs [coloc_query()] for each credible set in cs_list against its
##' own query results from [query_ot_many()] or
##' [query_ecat_many()], matched by name, and combines the results.
##' Credible sets whose query failed, or found nothing, are skipped.
##' Failures and warnings are handled as for [query_ot_many()].
##'
##' @inheritParams query_ot_many
##' @param query_list output of [query_ot_many()] or
##'     [query_ecat_many()] for the same cs_list
##' @return a data.table of coloc results, as from [coloc_query()],
##'     with an extra first column query_cs giving the name of the
##'     index credible set each row refers to. Attribute "failed" gives
##'     the error message for each credible set that failed, in the
##'     query or in the colocalisation.
##' @export
##' @author Chris Wallace
##' @seealso [query_ot_many()], [coloc_query()]
coloc_query_many <- function(cs_list, query_list, max_fail = 5) {
    if (is.null(names(cs_list))) names(cs_list) <- seq_along(cs_list) # as the query functions name them
    unknown <- setdiff(names(query_list), names(cs_list))
    if (is.null(names(query_list)) || length(unknown))
        stop("query_list should be the output of query_ot_many() or query_ecat_many() for the same cs_list: ",
             "names not found in cs_list: ", paste(utils::head(unknown, 5), collapse = ", "), call. = FALSE)
    keep <- names(query_list)[!vapply(query_list, is.null, logical(1))] # skip failed queries
    if (!length(keep)) stop("query_list contains no query results", call. = FALSE)
    res <- .run_many(cs_list[keep], function(cs, nm) coloc_query(cs, query_list[[nm]]), "coloc_query", max_fail)
    out <- data.table::rbindlist(res, idcol = "query_cs", fill = TRUE)
    data.table::setattr(out, "failed", c(attr(query_list, "failed"), attr(res, "failed"))) # failed queries, then failed colocs
    out
}

################################################################################

## ensembl to symbol
.hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
.gene_env <- new.env(parent = emptyenv())

## Ensembl id -> symbol table from HGNC, downloaded once and read once per session
.hgnc_genes <- function() {
    if (is.null(.gene_env$hgnc)) {
        path <- .download_once(.hgnc_url, .cache_dir("eqtl_catalogue")) # shared by all releases
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
##' @return character vector of gene symbols, the same length as ids, with NA where no symbol is found
##' @export
ensembl_to_symbol <- function(ids) {
    g <- .hgnc_genes()
    ids <- sub("\\..*","",ids)
    g$symbol[match(sub("\\.[0-9]+$", "", ids), g$ensembl_gene_id)]
}

## Local path to a release's dataset metadata, downloaded to the cache dir if not already there
.ecat_metadata_path <- function(version = .ecat_default_version) {
    rel <- .ecat_release(version)
    url <- paste0(.ecat_url(version), "/files/", rel$metadata, "?download=1")
    .download_once(url, ecat_dir(version), file = rel$metadata, md5 = rel$metadata_md5)
}

##' Load the eQTL Catalogue dataset metadata
##'
##' Returns the eQTL Catalogue dataset metadata for a release: one row
##' per dataset, giving the study, tissue or cell type (including
##' its ontology term, tissue_id), condition, sample size and
##' quantification method. This is the table that query_ecat()
##' merges into the "annot" attribute of its results, by dataset_id.
##'
##' The file is read from the cache directory if it is already there
##' (eg from an earlier query_ecat() or ecat_fetch() call),
##' and downloaded otherwise. It is read only once per R session. See
##' [coloc-options] for setting a persistent cache directory.
##'
##' @inheritParams query_ecat
##' @return a data.frame, one row per dataset
##' @export
##' @author Chris Wallace
##' @seealso [query_ecat()], [ecat_fetch()]
##' @examples
##' \dontrun{
##' meta <- ecat_metadata()
##' table(meta$tissue_label)
##' meta7 <- ecat_metadata("r7")
##' }
ecat_metadata <- function(version = .ecat_default_version) {
    key <- paste0("metadata_", version) # .ecat_env[[version]] already holds the file index
    if (is.null(.ecat_env[[key]]))
        .ecat_env[[key]] <- data.table::fread(.ecat_metadata_path(version), data.table = FALSE)
    .ecat_env[[key]]
}

################################################################################
## group labels: roll biosamples (tissues, cell types) up to coarse groups

## Open Targets biosample index: every UBERON/CL/EFO term with its ancestors (via is_a and part_of)
.ot_biosample_url <- sprintf("%s/platform/%s/output/biosample/part-00000-b13417e4-8fed-4fbd-b720-70130065798e-c000.snappy.parquet", .ot_s3_url, .ot_release)
.ot_biosample_md5 <- "9ccfc396daf483776ef0994ea1de7115"
.ot_biosample_size <- 8428035

## Terms with no useful ontology ancestors, assigned by hand
.group_overrides <- c(EFO_0005292 = "hematopoietic cell", # lymphoblastoid cell line
                      BTO_0000930 = "neural cell")
## When a term rolls up to several groups, the first of these that matches wins
.group_priority <- c("hematopoietic cell", "epithelial cell", "neural cell", "connective tissue cell", "adipose tissue")

## Biosample lookup, built once per session: ancestors (including the term itself) and names, keyed by id, plus target terms
.biosample_lookup <- function() {
    if (is.null(.ot_env$biosample)) {
        path <- .download_once(.ot_biosample_url, .ot_dir(), file = paste0("biosample_", .ot_release, ".parquet"), md5 = .ot_biosample_md5, size = .ot_biosample_size)
        bs <- as.data.frame(arrow::read_parquet(path, col_select = c("biosampleId", "biosampleName", "ancestors")))
        f <- system.file("extdata", "cellxgene_group_terms.tsv", package = "coloc")
        if (!nzchar(f)) stop("cellxgene_group_terms.tsv not found in inst/extdata")
        .ot_env$biosample <- list(anc = stats::setNames(Map(c, bs$biosampleId, lapply(bs$ancestors, as.character)), bs$biosampleId),
                                  name = stats::setNames(bs$biosampleName, bs$biosampleId),
                                  targets = utils::read.delim(f)$id)
    }
    .ot_env$biosample
}

## Roll one biosample id up to its most specific target term(s)
.rollup_one <- function(id, bs) {
    if (id %in% names(.group_overrides)) return(.group_overrides[[id]])
    hit <- intersect(bs$anc[[id]], bs$targets) # NULL for ids not in the index
    if (!length(hit)) return(NA_character_)
    hit <- hit[!vapply(hit, function(h) any(h %in% unlist(bs$anc[setdiff(hit, h)])), logical(1))] # drop hits that are ancestors of other hits
    nm <- unname(bs$name[hit])
    for (kp in .group_priority) if (kp %in% nm) return(kp)
    paste(sort(nm), collapse = "|") # unresolved tie: add to .group_overrides or .group_priority
}

## Group label for each of a vector of biosample ids (UBERON_0001157 or UBERON:0001157 form)
.group_label <- function(ids) {
    ids <- sub(":", "_", as.character(ids), fixed = TRUE)
    u <- unique(ids[!is.na(ids)])
    if (!length(u)) return(rep(NA_character_, length(ids)))
    bs <- .biosample_lookup()
    unname(stats::setNames(vapply(u, .rollup_one, character(1), bs = bs), u)[ids])
}

## Open Targets disease index: EFO/MONDO trait terms with their ancestors and therapeutic areas, used to group GWAS
.ot_disease_url <- sprintf("%s/platform/%s/output/disease/disease.parquet", .ot_s3_url, .ot_release)
.ot_disease_md5 <- "91474ea14a5792bc2d43d3fd3da88595"
.ot_disease_size <- 7150787

## GWAS are grouped by the Open Targets therapeutic area of their trait. Where a trait has several, the first of these
## wins: causes (cancer, infection) and specific categories (psychiatric, immune) before organ systems
.gwas_area_priority <- c("cancer or benign tumor", "infectious disease", "psychiatric disorder", "immune system disorder",
                         "cardiovascular disorder", "nervous system disorder", "endocrine system disorder", "nutritional or metabolic disease",
                         "respiratory or thoracic disease", "gastrointestinal disease", "hepatobiliary disorder", "urinary system disorder",
                         "reproductive system or breast disease", "pregnancy or perinatal disease", "musculoskeletal or connective tissue disease",
                         "integumentary system disorder", "hematologic disorder", "disorder of visual system", "disorder of ear",
                         "injury, poisoning or other complication")
## Traits in the "measurement" area (eg protein or metabolite levels) are grouped by the first of these they descend from,
## or "other measurement"
.gwas_measurement_priority <- c("protein measurement", "lipid or lipoprotein measurement", "metabolite measurement", "neuroimaging measurement",
                                "brain attribute", "microbiome measurement", "seropositivity measurement", "cardiovascular measurement",
                                "hematological measurement", "anthropometric measurement", "hormone measurement", "diet measurement",
                                "facial morphology trait")
## Areas used only when a trait has no disease area and is not a measurement
.gwas_other_areas <- c("phenotype", "biological_process", "genetic, familial or congenital disease", "medical procedure")

## Disease lookup, built once per session: names of each term's ancestors (including itself) and of its therapeutic areas
.disease_lookup <- function() {
    if (is.null(.ot_env$disease)) {
        path <- .download_once(.ot_disease_url, .ot_dir(), file = paste0("disease_", .ot_release, ".parquet"), md5 = .ot_disease_md5, size = .ot_disease_size)
        ds <- as.data.frame(arrow::read_parquet(path, col_select = c("id", "name", "ancestors", "therapeuticAreas")))
        nm <- stats::setNames(ds$name, ds$id)
        .ot_env$disease <- list(anc = stats::setNames(lapply(seq_along(ds$id), function(i) unname(nm[c(ds$id[i], as.character(ds$ancestors[[i]]))])), ds$id),
                                area = stats::setNames(lapply(ds$therapeuticAreas, function(a) unname(nm[as.character(a)])), ds$id))
    }
    .ot_env$disease
}

## Group label for one GWAS study, from its trait ids (diseaseIds)
.gwas_rollup_one <- function(ids, dl) {
    ids <- as.character(ids)
    areas <- unique(unlist(dl$area[ids]))
    for (a in .gwas_area_priority) if (a %in% areas) return(a)
    if ("measurement" %in% areas) {
        anc <- unique(unlist(dl$anc[ids]))
        for (m in .gwas_measurement_priority) if (m %in% anc) return(m)
        return("other measurement")
    }
    for (a in .gwas_other_areas) if (a %in% areas) return(a)
    NA_character_
}

## Group label for each of a vector of GWAS studyIds; study is a table of studyId and diseaseIds (list column)
.gwas_group_label <- function(study_ids, study) {
    u <- unique(study_ids)
    dl <- .disease_lookup()
    lab <- vapply(study$diseaseIds[match(u, study$studyId)], .gwas_rollup_one, character(1), dl = dl)
    unname(stats::setNames(lab, u)[study_ids])
}

#' Add column symbol for gene symbol
#'
#' for the gene id in the first of cols present in res, or the original id where there is none
##'  
##'
##' @param res coloc_query result
##' @param cols name of column holding ensembl_gene_id
##' @return res with additional column "symbol". if no symbol is found, the ensembl_gene_id is used in its place
##' @export 
##' @author Chris Wallace
add_symbol <- function(res, cols=c("ensembl_gene_id", "molecular_trait_object_id", "gene_id","geneId")) {
    gcol <- intersect(cols, names(res))
    if (!length(gcol)) {
        warning("no gene id column found, so symbol not added", call. = FALSE)
        return(res)
    }
    ids <- as.character(res[[gcol[1]]])
    res$symbol <- ensembl_to_symbol(ids)
    res$symbol[is.na(res$symbol)] <- ids[is.na(res$symbol)] # no symbol: keep the original id
    res
}

##' Add tissue / cell type group labels to colocalisation results
##'
##' Adds a column `group_label` that rolls each result's biosample (tissue or cell type) up to a coarse group, such as
##' "colon", "brain" or "hematopoietic cell", for summarising and plotting results (eg with [plot_query()]).
##'
##' Biosamples are UBERON, Cell Ontology or EFO terms. Each is mapped to its most specific ancestor (following is_a
##' and part_of relations, as given by the Open Targets biosample index) among the CZ CELLxGENE curated high-level
##' terms (tissue_general and cell_class lists from <https://github.com/chanzuckerberg/cellxgene-ontology-guide>).
##' A few terms with no useful ancestors, such as lymphoblastoid cell lines, are assigned by hand, and where a term
##' matches several groups, cell classes such as "hematopoietic cell" take precedence. Any remaining ties are
##' returned joined by "|". The Open Targets biosample index (~8 MB) is downloaded to the cache directory on first use.
##'
##' For eQTL Catalogue results, the biosample is the dataset's tissue_id. For Open Targets QTL studies it is the
##' study's biosampleId.
##'
##' Open Targets GWAS studies have no biosample, so annotate_ot() groups them by trait instead: by the Open Targets
##' therapeutic area of the study's trait (diseaseIds), such as "immune system disorder". Where a trait has several
##' areas, causes (cancer, infection) and specific categories (psychiatric, immune) take precedence over organ systems.
##' Traits in the "measurement" area, the large majority of Open Targets GWAS, are split further by type, such as
##' "protein measurement" or "metabolite measurement". The Open Targets disease index (~7 MB) is downloaded to the
##' cache directory on first use.
##'
##' Both functions also add a column `symbol`, the HGNC gene symbol for the Ensembl gene id (see
##' [ensembl_to_symbol()]), or the original id where no symbol is found. For Open Targets this is the study's geneId,
##' which is NA for GWAS studies. Column `label` is the point label used by [plot_query()]: the gene symbol,
##' or for Open Targets GWAS, the trait (traitFromSource, without any "(UKB data field ...)" suffix). annotate_ot() also adds `study_class`, "gwas" or "qtl", for
##' faceting plots.
##'
##' @param res output of [coloc_query()] (or [coloc_query_many()]) for credible sets from [query_ecat()] or
##'     [query_ot()] respectively
##' @inheritParams query_ecat
##' @return res, with columns group_label, symbol and label added, and for annotate_ot(), study_class (plus tissue_id,
##'     or the Open Targets study columns studyType, biosampleId, geneId and traitFromSource, if not already present)
##' @export
##' @author Chris Wallace
##' @seealso [plot_query()]
##' @examples
##' \dontrun{
##' res <- annotate_ecat(coloc_query(q, query_ecat(q)))
##' table(res$group_label)
##' res_ot <- annotate_ot(coloc_query(q, query_ot(q, type = "qtl")))
##' }
annotate_ecat <- function(res, version = .ecat_default_version) {
    stopifnot("res must be a data.frame from coloc_query() with column dataset_id" = is.data.frame(res) && "dataset_id" %in% names(res))
    if (!"tissue_id" %in% names(res)) {
        meta <- ecat_metadata(version)
        res$tissue_id <- meta$tissue_id[match(res$dataset_id, meta$dataset_id)]
    }
    res$group_label <- .group_label(res$tissue_id)
    res <- add_symbol(res) # r8-beta: molecular_trait_object_id, r7: gene_id
    res$label <- res$symbol
    res
}

##' @rdname annotate_ecat
##' @export
annotate_ot <- function(res) {
    stopifnot("res must be a data.frame from coloc_query() with column studyId" = is.data.frame(res) && "studyId" %in% names(res))
    miss <- setdiff(c("studyType", "biosampleId", "geneId", "traitFromSource"), names(res))
    study <- .ot_study(unique(res$studyId), cols = c(miss, "diseaseIds")) # only the study columns not already there, plus traits
    for (m in miss) res[[m]] <- study[[m]][match(res$studyId, study$studyId)]
    gw <- !is.na(res$studyType) & res$studyType == "gwas"
    res$study_class <- ifelse(gw, "gwas", "qtl")
    res$group_label <- NA_character_
    res$group_label[!gw] <- .group_label(res$biosampleId[!gw]) # QTL: by tissue / cell type
    if (any(gw)) res$group_label[gw] <- .gwas_group_label(res$studyId[gw], study) # GWAS: by trait
    res <- add_symbol(res, "geneId")
    res$label <- ifelse(gw, sub("\\s*\\(UKB data field [0-9]+\\)", "", res$traitFromSource), res$symbol)
    res
}
