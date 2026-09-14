args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: marlin_runtime_probe.R <runtime-probe.tsv>")
}

suppressPackageStartupMessages({
  library(keras)
  library(tensorflow)
  library(reticulate)
})

keras_backend <- keras::k_backend()
if (!identical(keras_backend, "tensorflow")) {
  stop(sprintf("expected keras tensorflow backend, observed %s", keras_backend))
}

# Force the Python TensorFlow binding to initialize rather than only loading the R package.
probe_tensor <- tensorflow::tf$constant(1L)
invisible(probe_tensor$numpy())

r_version <- paste(R.version$major, R.version$minor, sep = ".")
keras_version <- as.character(utils::packageVersion("keras"))
tensorflow_version <- as.character(tensorflow::tf$`__version__`)
python_version <- as.character(
  reticulate::py_eval("'.'.join(map(str, __import__('sys').version_info[:3]))")
)
gpus <- tensorflow::tf$config$list_physical_devices("GPU")
execution_backend <- if (length(gpus) > 0) "gpu" else "cpu"

probe <- data.frame(
  key = c(
    "R_version",
    "keras_version",
    "tensorflow_version",
    "python_version",
    "execution_backend"
  ),
  value = c(
    r_version,
    keras_version,
    tensorflow_version,
    python_version,
    execution_backend
  ),
  stringsAsFactors = FALSE
)

write.table(
  probe,
  file = args[[1]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  eol = "\n"
)
