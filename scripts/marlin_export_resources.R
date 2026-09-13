args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: marlin_export_resources.R <marlin_v1.features.RData> <output.txt>")
}
load(args[1])
if (!exists("betas_sub_names")) {
  stop("marlin_v1.features.RData does not define betas_sub_names")
}
if (length(betas_sub_names) != 357340) {
  stop("MARLIN v1 feature artifact must contain exactly 357340 identifiers")
}
if (length(unique(betas_sub_names)) != 357340) {
  stop("MARLIN v1 feature identifiers must be unique")
}
if (anyNA(betas_sub_names) || any(trimws(as.character(betas_sub_names)) == "")) {
  stop("MARLIN v1 feature identifiers must be non-empty")
}
writeLines(as.character(betas_sub_names), con = args[2], sep = "\n", useBytes = TRUE)
