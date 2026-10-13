options(stringsAsFactors = FALSE, width = 140)

phase10_root <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L && nzchar(args[[1L]])) {
    normalizePath(args[[1L]], mustWork = TRUE)
  } else {
    normalizePath("/home/zhenhanli/myocardial/scop/GSE214611_spatial", mustWork = TRUE)
  }
}

phase10_dirs <- function(root) {
  dirs <- list(
    objects = file.path(root, "objects", "phase_10"),
    metadata = file.path(root, "metadata"),
    results = file.path(root, "results", "phase_10"),
    logs = file.path(root, "logs"),
    reports = file.path(root, "reports")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

phase10_log_start <- function(root, stem) {
  dir.create(file.path(root, "logs"), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, "logs", paste0(stem, ".log"))
  con <- file(path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  assign(".phase10_log_connection", con, envir = .GlobalEnv)
  cat(sprintf("[%s] START %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), stem))
  invisible(path)
}

phase10_log_stop <- function(stem, status = "PASS") {
  cat(sprintf("[%s] END %s status=%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), stem, status))
  sink(type = "message")
  sink()
  con <- get0(".phase10_log_connection", envir = .GlobalEnv)
  if (!is.null(con)) close(con)
  invisible(status)
}

stamp <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), ..., "\n", sep = "")
}

require_file <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Required file is absent or empty: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

logcpm_sparse <- function(counts, scale_factor = 1e4) {
  lib <- Matrix::colSums(counts)
  scaling <- scale_factor / pmax(lib, 1)
  out <- counts %*% Matrix::Diagonal(x = scaling)
  out@x <- log1p(out@x)
  out
}

safe_row_z <- function(x) {
  x <- as.matrix(x)
  mu <- rowMeans(x)
  sig <- matrixStats::rowSds(x)
  sig[!is.finite(sig) | sig == 0] <- 1
  sweep(sweep(x, 1L, mu, "-"), 1L, sig, "/")
}

scop_plot_object <- function(x) {
  if (is.list(x) && !is.null(x$plot)) x$plot else x
}
