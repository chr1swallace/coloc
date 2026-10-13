## Regenerate inst/extdata/ot_credset_index_<release>.tsv for a new Open Targets release.
## For each credible set parquet file, records its size, md5 and the first and last (chromosome, position),
## which query_ot() uses to download only the files overlapping a query and to check the downloads.
## Run from the package root. Needs network access to the Open Targets S3 bucket.
## After running, update .ot_release (and .ot_manifest_url) in R/coloc-query.R to match.

release <- "26.06"
region <- "eu-west-1"
bucket <- "open-targets-public-data-releases"
prefix <- sprintf("platform/%s/output/credible_set/", release)

## file names, sizes and md5s from the bucket listing (plain https); a single page holds up to 1000 files
xml <- paste(readLines(sprintf("https://%s.s3.%s.amazonaws.com/?list-type=2&prefix=%s", bucket, region, prefix), warn = FALSE), collapse = "")
stopifnot(!grepl("<IsTruncated>true</IsTruncated>", xml, fixed = TRUE))
contents <- regmatches(xml, gregexpr("<Contents>.*?</Contents>", xml, perl = TRUE))[[1]]
tag <- function(x, name) sub(sprintf(".*<%s>(.*?)</%s>.*", name, name), "\\1", x, perl = TRUE)
files <- data.frame(key = tag(contents, "Key"), size = as.numeric(tag(contents, "Size")), md5 = gsub("&quot;|\"", "", tag(contents, "ETag")))
files <- files[grepl("\\.parquet$", files$key), ]
files$md5[grepl("-", files$md5)] <- NA # ETags of multipart uploads are not md5s; those files are checked by size only

## first and last (chromosome, position) of each file
fs <- arrow::S3FileSystem$create(anonymous = TRUE, region = region)
ends <- do.call(rbind, lapply(files$key, function(k) {
  t <- arrow::read_parquet(fs$path(paste0(bucket, "/", k)), col_select = c("chromosome", "position"), as_data_frame = FALSE)
  chr <- t$chromosome$as_vector(); pos <- t$position$as_vector(); n <- length(chr)
  data.frame(chr_first = chr[1], pos_first = pos[1], chr_last = chr[n], pos_last = pos[n])
}))
index <- cbind(data.frame(file = basename(files$key), size = files$size), ends, md5 = files$md5)

## query_ot() relies on files being sorted by (chromosome as text, position) and in order across files
chrs <- unique(c(rbind(index$chr_first, index$chr_last)))
k1 <- match(index$chr_first, chrs) * 1e10 + index$pos_first
k2 <- match(index$chr_last, chrs) * 1e10 + index$pos_last
stopifnot(all(k1 <= k2), all(head(k2, -1) <= tail(k1, -1)))

write.table(index, sprintf("inst/extdata/ot_credset_index_%s.tsv", release), sep = "\t", row.names = FALSE, quote = FALSE)
