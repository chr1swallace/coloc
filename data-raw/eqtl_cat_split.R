## Split eQTL Catalogue credible set files into smaller files sorted by genomic position, for upload to Zenodo.
##
## The eQTL Catalogue publishes one large parquet file per quantification method, each covering the whole
## genome, so every query_eqtl_cat() call needs the whole file. This script rewrites each file as a series of
## smaller files ("chunks"), so a query only needs the chunk(s) overlapping its region, as for Open Targets.
##
## - credible sets are sorted by chromosome (1-22, X, ...) then by their first variant's position, and rows
##   within a credible set by position
## - a credible set is never split across chunks
## - chunks may run across a chromosome boundary, to keep the number of files per release under Zenodo's limit
##   (100 files per record by default)
## - an index records, for each chunk and each chromosome it contains, the range of variant positions, so
##   query_eqtl_cat() can pick the chunks overlapping a query window
##
## Output, for each release, in <work>/<release>/:
##   eQTL_Catalogue_<release>_<quant_method>_<NNN>.parquet  the chunks
##   eQTL_Catalogue_index_<release>.tsv                     the index (also to be copied to inst/extdata)
##   <metadata file>                                  the dataset metadata, unchanged
## Upload the contents of each <work>/<release>/ directory to its own Zenodo record.
##
## Needs R with arrow and data.table, network access to zenodo.org, ~5 GB of disk for the source files plus
## about the same again for the output, and enough memory to hold one chromosome of the largest file at once.
## Run from anywhere: Rscript eqtl_cat_split.R [work_dir]

library(arrow)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[1] else "eqtl_cat_split"

releases <- list(
  r7 = list(record = "23350544",
            files = c(exon = "eQTL_Catalogue_r7_cs_exon_071026.parquet",
                      ge = "eQTL_Catalogue_r7_cs_ge_071026.parquet",
                      leafcutter = "eQTL_Catalogue_r7_cs_leafcutter_071026.parquet",
                      tx = "eQTL_Catalogue_r7_cs_tx_071026.parquet",
                      txrev = "eQTL_Catalogue_r7_cs_txrev_071026.parquet"),
            md5 = c(exon = "e1487d912bcd5acf2339f99cb802d920", ge = "9552ccd8ac3cbb0cd39c5714e3ff46c8",
                    leafcutter = "40e4b76d1df75ce7a1dd603916671fb8", tx = "e037630600b8bcf0ca17d3e3aac0fce3",
                    txrev = "0c9086bd8580774f60c1a03e73e740e2"),
            metadata = c(file = "dataset_metadata_r7.tsv", md5 = "7020d45bc577bff33f2904671a97eb03")),
  `r8-beta` = list(record = "22846842",
            files = c(exon = "eQTL_Catalogue_r8-beta_cs_exon_190926.parquet",
                      ge = "eQTL_Catalogue_r8-beta_cs_ge_190926.parquet",
                      leafcutter = "eQTL_Catalogue_r8-beta_cs_leafcutter_190926.parquet",
                      majiq = "eQTL_Catalogue_r8-beta_cs_majiq_190926.parquet",
                      tx = "eQTL_Catalogue_r8-beta_cs_tx_190926.parquet",
                      txrev = "eQTL_Catalogue_r8-beta_cs_txrev_190926.parquet"),
            md5 = c(exon = "5085e0b6c61b57ff8acd3aad9f7fc297", ge = "1c54f8dfab4ced90fe0da3154b9b03f7",
                    leafcutter = "e83fbe0b283c63c28f843db4e155060b", majiq = "3963dbb8ffc347b46d25f13983ef9d0a",
                    tx = "c11a60ba590e6cf24417de8935cfb699", txrev = "277bf7c2f592a2e48dbfab59618a4e27"),
            metadata = c(file = "dataset_metadata.tsv", md5 = "52b0305f56bd046eb5eeb2b6f2efa874"))
)

max_files <- 95        # per release, leaving room under Zenodo's 100-file limit for the index and metadata
min_chunk_bytes <- 15e6 # don't make chunks smaller than this, even if the file limit would allow it
row_group_rows <- 5e4  # row groups within a chunk, so local scans can also skip data by position

## Download url to path unless already there with the right md5
fetch <- function(url, path, md5) {
  if (file.exists(path) && unname(tools::md5sum(path)) == md5) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old))
  message("Downloading ", basename(path))
  utils::download.file(url, tmp <- paste0(path, ".part"), mode = "wb", quiet = TRUE)
  if (unname(tools::md5sum(tmp)) != md5) stop("Checksum mismatch for ", basename(path))
  file.rename(tmp, path)
  path
}

## Chromosome order used throughout: 1-22, X, Y, MT, then anything else alphabetically
chr_order <- function(chrs) {
  known <- c(1:22, "X", "Y", "MT")
  c(intersect(known, chrs), sort(setdiff(chrs, known)))
}

## Sort one chromosome's rows: credible sets by their first position, then rows within a set by position
sort_chromosome <- function(d) {
  d[, .cs_start := min(position), by = cs_id]
  setorderv(d, c(".cs_start", "cs_id", "position"))
  d[, .cs_start := NULL]
  d
}

## Row index at which to cut buffer d (sorted, credible sets contiguous) to take about n rows without splitting
## a credible set: the end of the first credible set ending at or after row n. Returns nrow(d) if none.
cut_point <- function(d, n) {
  ends <- which(c(d$cs_id[-1] != d$cs_id[-nrow(d)], TRUE)) # last row of each credible set
  cut <- ends[ends >= n][1]
  if (is.na(cut)) nrow(d) else cut
}

## Write one chunk, returning its index rows (one per chromosome in the chunk)
write_chunk <- function(d, path) {
  arrow::write_parquet(d, path, chunk_size = row_group_rows)
  idx <- d[, .(pos_min = min(position), pos_max = max(position), n_cs = uniqueN(cs_id), n_rows = .N), by = chromosome]
  idx[, `:=`(file = basename(path), size = file.size(path), md5 = unname(tools::md5sum(path)))]
  idx[]
}

## Split one source file into chunks of about chunk_rows rows; returns index rows
split_file <- function(src, out_dir, prefix, chunk_rows) {
  ds <- arrow::open_dataset(src, format = "parquet")
  chrs <- chr_order(unique(as.character(arrow::Scanner$create(ds, projection = "chromosome")$ToTable()$chromosome$as_vector())))
  buffer <- NULL; k <- 0; index <- list()
  emit <- function(d) {
    k <<- k + 1
    index[[k]] <<- write_chunk(d, file.path(out_dir, sprintf("%s_%03d.parquet", prefix, k)))
  }
  for (chr in chrs) {
    filt <- arrow::Expression$field_ref("chromosome") == chr
    d <- as.data.table(arrow::Scanner$create(ds, filter = filt)$ToTable())
    d[, chromosome := as.character(chromosome)]
    buffer <- rbind(buffer, sort_chromosome(d))
    while (!is.null(buffer) && nrow(buffer) >= chunk_rows) {
      cut <- cut_point(buffer, chunk_rows)
      emit(buffer[seq_len(cut)])
      buffer <- if (cut < nrow(buffer)) buffer[-seq_len(cut)] else NULL
    }
  }
  if (!is.null(buffer) && nrow(buffer)) emit(buffer)
  rbindlist(index)
}

## Check a release's output: every credible set in exactly one chunk, all rows present, chunks within the file limit
check_release <- function(out_dir, index, srcs) {
  n_src <- sum(vapply(srcs, function(s) nrow(arrow::open_dataset(s)), numeric(1)))
  stopifnot("rows lost or duplicated" = sum(index$n_rows) == n_src)
  files <- unique(index$file)
  cs <- rbindlist(lapply(files, function(f) data.table(file = f, cs_id = unique(arrow::read_parquet(file.path(out_dir, f), col_select = "cs_id")$cs_id))))
  stopifnot("a credible set is split across chunks" = !anyDuplicated(cs$cs_id))
  message(length(files), " chunks, ", round(sum(unique(index[, .(file, size)])$size) / 1e9, 2), " GB")
  if (length(files) + 2 > 100) warning("More than 100 files: Zenodo's default per-record limit. Increase min_chunk_bytes or request a higher limit.")
}

for (rel in names(releases)) {
  r <- releases[[rel]]
  src_dir <- file.path(work, "source", rel); out_dir <- file.path(work, rel)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- sprintf("https://zenodo.org/records/%s/files/", r$record)

  srcs <- mapply(function(f, m) fetch(paste0(base, f, "?download=1"), file.path(src_dir, f), m), r$files, r$md5)
  meta <- fetch(paste0(base, r$metadata[["file"]], "?download=1"), file.path(out_dir, r$metadata[["file"]]), r$metadata[["md5"]])

  ## one chunk size for the whole release, chosen to keep the number of files under max_files
  chunk_bytes <- max(min_chunk_bytes, sum(file.size(srcs)) / (max_files - length(srcs)))
  index <- rbindlist(lapply(names(srcs), function(qm) {
    message(rel, " ", qm)
    bytes_per_row <- file.size(srcs[[qm]]) / nrow(arrow::open_dataset(srcs[[qm]]))
    idx <- split_file(srcs[[qm]], out_dir, sprintf("eQTL_Catalogue_%s_%s", rel, qm), chunk_rows = chunk_bytes / bytes_per_row)
    idx[, quant_method := qm][]
  }))
  setcolorder(index, c("quant_method", "file", "chromosome", "pos_min", "pos_max", "n_cs", "n_rows", "size", "md5"))
  fwrite(index, file.path(out_dir, sprintf("eQTL_Catalogue_index_%s.tsv", rel)), sep = "\t")
  check_release(out_dir, index, srcs)
}
