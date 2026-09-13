args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: marlin_infer_locked.R <feature-vector.txt> <model.hdf5> <scores.tsv>")
}

suppressPackageStartupMessages(library(keras))

values <- scan(args[1], what = double(), quiet = TRUE)
if (length(values) != 357340) {
  stop("MARLIN v1 inference requires exactly 357340 features")
}
if (any(!values %in% c(-1, 0, 1))) {
  stop("MARLIN feature vector contains values outside {-1,0,1}")
}

model <- load_model_hdf5(args[2])
pred <- as.numeric(predict(model, matrix(values, nrow = 1), verbose = 0))
if (length(pred) != 42) {
  stop("MARLIN v1 model must emit exactly 42 scores")
}
if (any(!is.finite(pred)) || any(pred < 0) || any(pred > 1)) {
  stop("MARLIN model emitted invalid probability scores")
}

write.table(
  data.frame(model_id = seq_along(pred), score = pred),
  file = args[3],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
