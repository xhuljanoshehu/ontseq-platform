#!/usr/bin/env Rscript

# ONTSeq compatibility wrapper for the ichorCNA 0.5.1 R package.
#
# Scientific implementation: ichorCNA::run_ichorCNA
# This file only exposes a stable argument-vector boundary for ONTSeq.
# Research Use Only; analytical validation remains separate.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--id", type = "character", dest = "id"),
  make_option("--WIG", type = "character", dest = "WIG"),
  make_option("--ploidy", type = "character", dest = "ploidy"),
  make_option("--normal", type = "character", dest = "normal"),
  make_option("--maxCN", type = "integer", dest = "maxCN"),
  make_option("--gcWig", type = "character", dest = "gcWig"),
  make_option("--mapWig", type = "character", dest = "mapWig"),
  make_option("--centromere", type = "character", dest = "centromere"),
  make_option("--minMapScore", type = "double", dest = "minMapScore"),
  make_option("--includeHOMD", type = "logical", dest = "includeHOMD"),
  make_option("--chrs", type = "character", dest = "chrs"),
  make_option("--chrTrain", type = "character", dest = "chrTrain"),
  make_option("--chrNormalize", type = "character", dest = "chrNormalize"),
  make_option("--genomeBuild", type = "character", dest = "genomeBuild"),
  make_option("--genomeStyle", type = "character", dest = "genomeStyle"),
  make_option("--estimateNormal", type = "logical", dest = "estimateNormal"),
  make_option("--estimatePloidy", type = "logical", dest = "estimatePloidy"),
  make_option(
    "--estimateScPrevalence",
    type = "logical",
    dest = "estimateScPrevalence"
  ),
  make_option("--txnE", type = "double", dest = "txnE"),
  make_option("--txnStrength", type = "double", dest = "txnStrength"),
  make_option("--outDir", type = "character", dest = "outDir"),
  make_option(
    "--normalPanel",
    type = "character",
    dest = "normalPanel",
    default = NULL
  )
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

required <- c(
  "id",
  "WIG",
  "ploidy",
  "normal",
  "maxCN",
  "gcWig",
  "mapWig",
  "centromere",
  "minMapScore",
  "includeHOMD",
  "chrs",
  "chrTrain",
  "chrNormalize",
  "genomeBuild",
  "genomeStyle",
  "estimateNormal",
  "estimatePloidy",
  "estimateScPrevalence",
  "txnE",
  "txnStrength",
  "outDir"
)

missing <- required[vapply(required, function(name) is.null(opt[[name]]), logical(1))]
if (length(missing) > 0) {
  stop(
    paste0(
      "Missing required ichorCNA wrapper option(s): ",
      paste(missing, collapse = ", ")
    )
  )
}

ichorCNA::run_ichorCNA(
  tumor_wig = opt$WIG,
  gcWig = opt$gcWig,
  mapWig = opt$mapWig,
  normal_panel = opt$normalPanel,
  id = opt$id,
  centromere = opt$centromere,
  minMapScore = opt$minMapScore,
  normal = opt$normal,
  ploidy = opt$ploidy,
  maxCN = opt$maxCN,
  estimateNormal = opt$estimateNormal,
  estimatePloidy = opt$estimatePloidy,
  estimateScPrevalence = opt$estimateScPrevalence,
  chrNormalize = opt$chrNormalize,
  chrTrain = opt$chrTrain,
  chrs = opt$chrs,
  genomeBuild = opt$genomeBuild,
  genomeStyle = opt$genomeStyle,
  includeHOMD = opt$includeHOMD,
  txnE = opt$txnE,
  txnStrength = opt$txnStrength,
  outDir = opt$outDir
)
