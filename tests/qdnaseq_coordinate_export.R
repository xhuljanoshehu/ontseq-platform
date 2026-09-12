# Synthetic boundary regression for the actual R export functions; no analysis is run.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
functions <- new.env(parent=globalenv())
functions$fail <- function(message, code=2L) stop(message)
for (expression in parse(file=args[[1L]])) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      as.character(expression[[2L]]) %in% c("qdnaseq_export_intervals", "collapse_segments")) {
    eval(expression, envir=functions)
  }
}

native <- data.frame(
  chromosome=c("7", "7", "7", "8", "8", "8"),
  start=c(1L, 501L, 1000L, 1L, 501L, 1000L),
  end=c(500L, 999L, 1000L, 500L, 999L, 1000L)
)
called <- data.frame(
  bin=seq_len(nrow(native)), chr=native$chromosome,
  segments=c(1, 3, 1, 3, 3, 3), calls=c(-1, 1, -1, 1, 1, 1)
)
exported <- functions$qdnaseq_export_intervals(native)
segments <- functions$collapse_segments(called, native)
stopifnot(identical(native$start, c(1L, 501L, 1000L, 1L, 501L, 1000L)))
stopifnot(identical(exported$start, c(0L, 500L, 999L, 0L, 500L, 999L)))
stopifnot(identical(exported$end, native$end))
stopifnot(all(exported$end - exported$start == native$end - native$start + 1L))
stopifnot(all(exported$coordinate_system == "zero_based_half_open"))
stopifnot(identical(as.integer(segments$start), c(0L, 500L, 999L, 0L)))
stopifnot(identical(as.integer(segments$end), c(500L, 999L, 1000L, 1000L)))
stopifnot(identical(as.integer(segments$bin_count), c(1L, 1L, 1L, 3L)))
stopifnot(all(segments$coordinate_system == "zero_based_half_open"))
for (bad_start in c(0, -1, 1.5, NA_real_, Inf)) {
  bad <- native
  bad$start[[1L]] <- bad_start
  stopifnot(inherits(try(functions$qdnaseq_export_intervals(bad), silent=TRUE), "try-error"))
}
bad <- native
bad$end[[1L]] <- 0L
stopifnot(inherits(try(functions$qdnaseq_export_intervals(bad), silent=TRUE), "try-error"))
stopifnot(inherits(try(functions$qdnaseq_export_intervals(exported), silent=TRUE), "try-error"))

# If installed, verify both pinned annotation lanes use the expected native convention.
if (requireNamespace("QDNAseq", quietly=TRUE)) {
  for (genome in c("hg19", "hg38")) {
    if (!requireNamespace(paste0("QDNAseq.", genome), quietly=TRUE)) next
    bins <- QDNAseq::getBinAnnotations(binSize=500, genome=genome)
    fd <- Biobase::pData(bins)
    corrected <- functions$qdnaseq_export_intervals(fd)
    chr7 <- which(as.character(fd$chromosome) == "7")
    stopifnot(fd$start[chr7[[1L]]] == 1L, corrected$start[chr7[[1L]]] == 0L)
    stopifnot(all(corrected$end == fd$end))
    stopifnot(all(corrected$end - corrected$start == fd$end - fd$start + 1L))
  }
}

dir.create(args[[2L]], recursive=TRUE, showWarnings=FALSE)
write.table(exported, file.path(args[[2L]], "bins.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(segments, file.path(args[[2L]], "segments.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")
cat("QDNAseq coordinate export regression PASS\n")
