options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

now_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE)
sha256_value <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

pkg_version_or_na <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
}

write_session_info <- function(path) {
  capture.output(print(utils::sessionInfo()), file = path)
}

phase_status <- function(phase, status, started_at, input_hash, parameter_hash,
                         n_expected, n_completed, producer_functions,
                         biological_verdict, blocker, output_root) {
  data.frame(
    phase = phase,
    status = status,
    started_at = started_at,
    completed_at = now_utc(),
    input_hash = input_hash,
    parameter_hash = parameter_hash,
    n_samples_expected = n_expected,
    n_samples_completed = n_completed,
    inference_unit = "physical_tissue/biological_sample",
    spot_level_inference = FALSE,
    producer_functions = producer_functions,
    biological_verdict = biological_verdict,
    blocker = blocker,
    output_root = output_root,
    stringsAsFactors = FALSE
  )
}

write_phase_status <- function(out, ...) {
  write_tsv(phase_status(...), file.path(out, "STATUS.tsv"))
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, colnames(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

empty_df <- function(columns) {
  out <- setNames(replicate(length(columns), character(), simplify = FALSE), columns)
  as.data.frame(out, stringsAsFactors = FALSE)
}
