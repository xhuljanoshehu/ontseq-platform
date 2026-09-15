args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: marlin_runtime_probe.R <runtime-probe.tsv>")
}

runtime_prefix <- normalizePath(file.path(R.home(), "..", ".."), mustWork = TRUE)
runtime_python <- file.path(runtime_prefix, "bin", "python")
if (!file.exists(runtime_python)) {
  stop(sprintf("MARLIN runtime Python is missing: %s", runtime_python))
}
Sys.setenv(RETICULATE_PYTHON = runtime_python)

suppressPackageStartupMessages({
  library(reticulate)
})
reticulate::use_python(runtime_python, required = TRUE)
python_config <- reticulate::py_config()
if (normalizePath(python_config$python, mustWork = TRUE) != normalizePath(runtime_python, mustWork = TRUE)) {
  stop("reticulate did not bind to the MARLIN runtime Python")
}

suppressPackageStartupMessages({
  library(keras)
  library(tensorflow)
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
