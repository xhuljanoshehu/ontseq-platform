#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(QDNAseq)
  library(ACE)
  library(Biobase)
  library(ggplot2)
  library(jsonlite)
})

# QDNAseq::segmentBins() delegates to DNAcopy::segment(), which uses RNG and explicitly
# documents that identical inputs can otherwise yield slightly different segmentations.
# Keep this technical seed fixed and record it in the emitted summary. The R script SHA is
# part of the ONTSeq CNV stage signature, so changing this constant invalidates resume.
SEGMENTATION_SEED <- 104729L

fail <- function(message, code = 2L) {
  cat(sprintf("ERROR: %s\n", message), file = stderr())
  quit(save = "no", status = code)
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      fail(sprintf("invalid argument sequence near '%s'", key))
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

required <- function(opts, name) {
  value <- opts[[name]]
  if (is.null(value) || !nzchar(value)) fail(sprintf("missing --%s", name))
  value
}

as_int_vector <- function(value, name) {
  parsed <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
  if (!length(parsed) || any(is.na(parsed)) || any(parsed <= 0)) {
    fail(sprintf("--%s must be a comma-separated list of positive integers", name))
  }
  unique(parsed)
}

safe_pkg_version <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

collapse_segments <- function(called_template, feature_data) {
  df <- called_template
  if (!nrow(df)) return(data.frame())
  idx <- as.integer(df$bin)
  if (any(is.na(idx)) || any(idx < 1L) || any(idx > nrow(feature_data))) {
    fail("ACE returned bin indices outside the QDNAseq feature table")
  }
  df$start <- as.integer(feature_data$start[idx])
  df$end <- as.integer(feature_data$end[idx])
  df$chr <- gsub("^chr", "", as.character(df$chr), ignore.case = TRUE)
  keep <- !is.na(df$segments) & !is.na(df$start) & !is.na(df$end) & df$chr %in% c(as.character(1:22))
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) return(data.frame())

  # ACEcall names its purity/ploidy-adjusted absolute segment estimate `segments`.
  # Group consecutive bins with the same chromosome, call and adjusted estimate.
  seg_key <- round(as.numeric(df$segments), 6)
  call_key <- ifelse(is.na(df$calls), 0, as.numeric(df$calls))
  boundary <- c(
    TRUE,
    df$chr[-1L] != df$chr[-nrow(df)] |
      seg_key[-1L] != seg_key[-nrow(df)] |
      call_key[-1L] != call_key[-nrow(df)]
  )
  group <- cumsum(boundary)
  parts <- split(seq_len(nrow(df)), group)
  rows <- lapply(parts, function(ix) {
    qvals <- if ("qnorm_log10" %in% names(df)) as.numeric(df$qnorm_log10[ix]) else rep(NA_real_, length(ix))
    data.frame(
      chromosome = paste0("chr", df$chr[ix[[1L]]]),
      start = min(df$start[ix], na.rm = TRUE),
      end = max(df$end[ix], na.rm = TRUE),
      bin_count = length(ix),
      absolute_copy_number = stats::median(as.numeric(df$segments[ix]), na.rm = TRUE),
      call = stats::median(call_key[ix], na.rm = TRUE),
      qnorm_log10 = if (all(is.na(qvals))) NA_real_ else stats::median(qvals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

weighted_chromosome_summary <- function(segments) {
  if (!nrow(segments)) return(data.frame())
  chromosomes <- unique(segments$chromosome)
  rows <- lapply(chromosomes, function(chr) {
    part <- segments[segments$chromosome == chr, , drop = FALSE]
    widths <- pmax(1, part$end - part$start)
    expanded <- rep(part$absolute_copy_number, pmax(1L, round(widths / max(1, min(widths)))))
    data.frame(
      chromosome = chr,
      copy_number = stats::median(expanded, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
bam <- normalizePath(required(opts, "bam"), mustWork = TRUE)
out_dir <- required(opts, "output-dir")
sample_id <- required(opts, "sample-id")
genome_build <- required(opts, "genome-build")
bin_sizes <- as_int_vector(required(opts, "bin-sizes-kbp"), "bin-sizes-kbp")
primary_bin <- suppressWarnings(as.integer(required(opts, "primary-bin-kbp")))
penalty <- suppressWarnings(as.numeric(required(opts, "ace-penalty")))
ploidy_min <- suppressWarnings(as.numeric(required(opts, "ploidy-min")))
ploidy_max <- suppressWarnings(as.numeric(required(opts, "ploidy-max")))
ploidy_step <- suppressWarnings(as.numeric(required(opts, "ploidy-step")))
threads <- suppressWarnings(as.integer(required(opts, "threads")))

if (is.na(primary_bin) || !(primary_bin %in% bin_sizes)) fail("primary bin must be one of the configured bin sizes")
if (is.na(penalty) || penalty < 0 || penalty > 1) fail("ACE penalty must be in [0,1]")
if (any(is.na(c(ploidy_min, ploidy_max, ploidy_step))) || ploidy_min <= 0 || ploidy_max <= ploidy_min || ploidy_step <= 0) {
  fail("invalid ploidy search range")
}
if (is.na(threads) || threads < 1L) fail("threads must be >= 1")

genome <- switch(
  genome_build,
  "GRCh37" = "hg19",
  "GRCh38" = "hg38",
  fail(sprintf("unsupported genome build '%s'", genome_build))
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

if (requireNamespace("future", quietly = TRUE)) {
  workers <- max(1L, threads)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
}

runs <- list()
chromosome_by_bin <- list()

for (bin_size in bin_sizes) {
  message(sprintf("QDNAseq/ACE: %s, bin=%s kbp", sample_id, bin_size))
  bins <- tryCatch(
    QDNAseq::getBinAnnotations(binSize = bin_size, genome = genome),
    error = function(e) fail(sprintf("cannot load QDNAseq %s/%s-kbp annotations: %s", genome, bin_size, conditionMessage(e)))
  )
  # QDNAseq defines chunkSize as an optional integer number of nucleotides. NULL is
  # the package default and avoids accidentally coercing a logical TRUE to a 1-nt chunk.
  read_counts <- QDNAseq::binReadCounts(bins, bamfiles = bam, chunkSize = NULL)
  filtered <- QDNAseq::applyFilters(read_counts, residual = TRUE, blacklist = TRUE)
  filtered <- QDNAseq::estimateCorrection(filtered)
  copy_numbers <- QDNAseq::correctBins(filtered)
  copy_numbers <- QDNAseq::normalizeBins(copy_numbers)
  copy_numbers <- QDNAseq::smoothOutlierBins(copy_numbers)

  # Reset before every resolution so a result depends on its own data and this explicit
  # technical seed, not on bin-size iteration order or RNG consumed by an earlier run.
  set.seed(SEGMENTATION_SEED)
  segmented <- QDNAseq::segmentBins(copy_numbers, transformFun = "sqrt")
  segmented <- QDNAseq::normalizeSegmentedBins(segmented)

  rds_path <- file.path(out_dir, sprintf("%s.%skbp.segmented.rds", sample_id, bin_size))
  saveRDS(segmented, rds_path)

  template <- ACE::objectsampletotemplate(segmented, index = 1)
  autosomal_segments <- as.numeric(template$segments[!template$chr %in% c("X", "Y")])
  autosomal_segments <- autosomal_segments[is.finite(autosomal_segments)]
  if (!length(autosomal_segments)) fail(sprintf("QDNAseq produced no finite autosomal segments for %s kbp", bin_size))
  normalized_segment_median <- stats::median(autosomal_segments)
  message(sprintf("QDNAseq normalized segment median: %.6f", normalized_segment_median))

  prows <- max(1L, as.integer(round((ploidy_max - ploidy_min) / ploidy_step)))
  model <- ACE::squaremodel(
    template,
    ptop = ploidy_max,
    pbottom = ploidy_min,
    prows = prows,
    method = "RMSE",
    exclude = c("X", "Y"),
    penalty = penalty,
    penploidy = 0,
    cellularities = seq(5, 100),
    highlightminima = TRUE,
    standard = 1
  )
  candidates <- model$minimadf
  candidates <- candidates[is.finite(candidates$error), , drop = FALSE]
  if (!nrow(candidates)) fail(sprintf("ACE produced no finite model candidate for %s kbp", bin_size))
  candidates$distance_to_diploid <- abs(candidates$ploidy - 2)
  candidates <- candidates[order(candidates$error, -candidates$cellularity, candidates$distance_to_diploid), , drop = FALSE]
  chosen <- candidates[1L, , drop = FALSE]

  # squaremodel() accepts cellularities as percentages (5..100) but stores them in
  # errordf/minimadf as fractions (0.05..1.00). Do not divide the returned value again.
  cellularity_fraction <- as.numeric(chosen$cellularity[[1L]])
  chosen_ploidy <- as.numeric(chosen$ploidy[[1L]])
  if (!is.finite(cellularity_fraction) || cellularity_fraction <= 0 || cellularity_fraction > 1) {
    fail(sprintf("ACE returned invalid cellularity %.8g for %s kbp", cellularity_fraction, bin_size))
  }
  if (!is.finite(chosen_ploidy) || chosen_ploidy <= 0) {
    fail(sprintf("ACE returned invalid ploidy %.8g for %s kbp", chosen_ploidy, bin_size))
  }
  message(
    sprintf(
      "ACE selected fit: cellularity=%.3f, ploidy=%.3f, error=%.6g",
      cellularity_fraction,
      chosen_ploidy,
      as.numeric(chosen$error[[1L]])
    )
  )

  called <- ACE::ACEcall(
    segmented,
    QDNAseqobjectsample = 1,
    cellularity = cellularity_fraction,
    ploidy = chosen_ploidy,
    standard = 1,
    plot = TRUE,
    onlyautosomes = TRUE
  )
  adjusted_segments <- as.numeric(called$calledtemplate$segments)
  adjusted_segments <- adjusted_segments[is.finite(adjusted_segments)]
  if (!length(adjusted_segments)) fail(sprintf("ACE returned no finite adjusted segments for %s kbp", bin_size))
  message(
    sprintf(
      "ACE adjusted segment range: %.3f..%.3f",
      min(adjusted_segments),
      max(adjusted_segments)
    )
  )

  fd <- Biobase::fData(segmented)
  bins_path <- file.path(out_dir, sprintf("%s.%skbp.bins.tsv", sample_id, bin_size))
  utils::write.table(fd, bins_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  segments <- collapse_segments(called$calledtemplate, fd)
  segment_path <- file.path(out_dir, sprintf("%s.%skbp.segments.tsv", sample_id, bin_size))
  utils::write.table(segments, segment_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  chromosome_summary <- weighted_chromosome_summary(segments)
  chromosome_path <- file.path(out_dir, sprintf("%s.%skbp.chromosomes.tsv", sample_id, bin_size))
  utils::write.table(chromosome_summary, chromosome_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  chromosome_by_bin[[as.character(bin_size)]] <- chromosome_summary

  fit_plot <- file.path(out_dir, sprintf("%s.%skbp.ace-fit.png", sample_id, bin_size))
  cn_plot <- file.path(out_dir, sprintf("%s.%skbp.copy-number.png", sample_id, bin_size))
  ggplot2::ggsave(fit_plot, model$matrixplot, width = 12, height = 7, dpi = 140)
  ggplot2::ggsave(cn_plot, called$calledplot, width = 14, height = 7, dpi = 140)

  candidate_rows <- head(candidates, 12L)
  model_path <- file.path(out_dir, sprintf("%s.%skbp.ace-models.tsv", sample_id, bin_size))
  utils::write.table(
    candidates,
    model_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  runs[[length(runs) + 1L]] <- list(
    bin_size_kbp = bin_size,
    cellularity = cellularity_fraction,
    ploidy = chosen_ploidy,
    fit_error = as.numeric(chosen$error[[1L]]),
    candidate_count = nrow(candidates),
    alternatives = lapply(seq_len(nrow(candidate_rows)), function(i) list(
      cellularity = as.numeric(candidate_rows$cellularity[[i]]),
      ploidy = as.numeric(candidate_rows$ploidy[[i]]),
      fit_error = as.numeric(candidate_rows$error[[i]])
    )),
    segment_file = basename(segment_path),
    chromosome_file = basename(chromosome_path),
    bins_file = basename(bins_path),
    model_file = basename(model_path),
    fit_plot = basename(fit_plot),
    copy_number_plot = basename(cn_plot),
    rds_file = basename(rds_path),
    segment_count = nrow(segments)
  )
}

all_chr <- sort(unique(unlist(lapply(chromosome_by_bin, function(x) x$chromosome))))
consensus_rows <- lapply(all_chr, function(chr) {
  values <- vapply(chromosome_by_bin, function(tbl) {
    hit <- tbl[tbl$chromosome == chr, "copy_number", drop = TRUE]
    if (!length(hit)) NA_real_ else as.numeric(hit[[1L]])
  }, numeric(1))
  values <- values[is.finite(values)]
  if (!length(values)) return(NULL)
  data.frame(
    chromosome = chr,
    median_copy_number = stats::median(values),
    rounded_copy_number = round(stats::median(values)),
    agreeing_bins = sum(round(values) == round(stats::median(values))),
    contributing_bins = length(values),
    min_copy_number = min(values),
    max_copy_number = max(values),
    stringsAsFactors = FALSE
  )
})
consensus_rows <- Filter(Negate(is.null), consensus_rows)
consensus <- if (length(consensus_rows)) do.call(rbind, consensus_rows) else data.frame()
consensus_path <- file.path(out_dir, sprintf("%s.consensus.chromosomes.tsv", sample_id))
utils::write.table(consensus, consensus_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

primary <- runs[[which(vapply(runs, function(x) x$bin_size_kbp == primary_bin, logical(1)))[1L]]]
summary <- list(
  schema_version = "0.1.0",
  sample_id = sample_id,
  genome_build = genome_build,
  genome_annotation = genome,
  primary_bin_size_kbp = primary_bin,
  segmentation_seed = SEGMENTATION_SEED,
  ace_penalty = penalty,
  ploidy_search = list(min = ploidy_min, max = ploidy_max, step = ploidy_step),
  consensus_strategy = "median_rounded_across_bins",
  package_versions = list(
    R = paste(R.version$major, R.version$minor, sep = "."),
    QDNAseq = safe_pkg_version("QDNAseq"),
    ACE = safe_pkg_version("ACE"),
    DNAcopy = safe_pkg_version("DNAcopy"),
    QDNAseq_hg19 = safe_pkg_version("QDNAseq.hg19"),
    QDNAseq_hg38 = safe_pkg_version("QDNAseq.hg38")
  ),
  runs = runs,
  primary = primary,
  consensus_file = basename(consensus_path)
)
summary_path <- file.path(out_dir, sprintf("%s.qdnaseq-ace.summary.json", sample_id))
jsonlite::write_json(summary, summary_path, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 8)
message(sprintf("Wrote %s", summary_path))
