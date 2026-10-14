#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200, warn = 1)
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: install_spacet.R LIB RUNTIME_ROOT")
lib <- args[[1]]
runtime <- args[[2]]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(runtime, "objects"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

write_status <- function(status, reason = "") {
  x <- data.frame(
    component = "SpaCET", status = status,
    expected_version = "1.4.0",
    expected_commit = "0199d93c23c6fce142beb47b49c76e8b673100ad",
    library = normalizePath(lib, mustWork = FALSE), reason = reason,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write.table(x, file.path(runtime, "STATUS.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
  capture.output(sessionInfo(), file = file.path(runtime, "sessionInfo.txt"))
}

tryCatch({
  if (!requireNamespace("BiocManager", quietly = TRUE)) stop("BiocManager is unavailable")
  if (!requireNamespace("remotes", quietly = TRUE)) stop("remotes is unavailable")
  if (!requireNamespace("BiRewire", quietly=TRUE)) BiocManager::install("BiRewire", lib = lib, version = "3.22", update = FALSE,
                       ask = FALSE, dependencies = NA)
  if (!identical(as.character(packageVersion("BiRewire", lib.loc = lib)), "3.41.0")) {
    stop("BiRewire did not resolve to 3.41.0")
  }
  if (!requireNamespace("MUDAN", quietly=TRUE)) remotes::install_github(
    "JEFworks/MUDAN", ref = "3636ed2f9a404cbdfcfdce5500384b59a72de666",
    lib = lib, dependencies = NA, upgrade = "never", force = TRUE
  )
  mudan_desc <- packageDescription("MUDAN", lib.loc = lib)
  if (!identical(mudan_desc$RemoteSha, "3636ed2f9a404cbdfcfdce5500384b59a72de666")) {
    stop("MUDAN installed commit does not match the fixed commit")
  }
  remotes::install_github(
    "data2intelligence/SpaCET", ref = "0199d93c23c6fce142beb47b49c76e8b673100ad",
    lib = lib, dependencies = NA, upgrade = "never", force = TRUE
  )
  suppressPackageStartupMessages(library(SpaCET, lib.loc = lib))
  if (!identical(as.character(packageVersion("SpaCET")), "1.4.0")) stop("SpaCET version is not 1.4.0")
  desc <- packageDescription("SpaCET")
  if (!identical(desc$RemoteSha, "0199d93c23c6fce142beb47b49c76e8b673100ad")) {
    stop("SpaCET installed commit does not match the fixed commit")
  }
  if (!is.function(SpaCET::create.SpaCET.object)) stop("create.SpaCET.object is not callable")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 0, 3, 1, 2, 1, 0), nrow = 3,
           dimnames = list(c("A1BG", "A2M", "ACE"), c("spot1", "spot2", "spot3"))),
    sparse = TRUE
  )
  coords <- data.frame(
    pixel_row = c(10, 20, 30), pixel_col = c(10, 20, 10),
    coordinate_x_um = c(0, 100, 50), coordinate_y_um = c(0, 0, 86.60254),
    row.names = colnames(counts)
  )
  meta <- data.frame(smoke_group = c("a", "a", "b"), row.names = colnames(counts))
  obj <- SpaCET::create.SpaCET.object(
    counts = counts, spotCoordinates = coords, metaData = meta,
    imagePath = NA, platform = "Visium", organism = "human"
  )
  if (!inherits(obj, "SpaCET")) stop("SpaCET smoke constructor returned the wrong class")
  smoke_path <- file.path(runtime, "objects/spacet_smoke.rds")
  saveRDS(obj, smoke_path)
  if (!inherits(readRDS(smoke_path), "SpaCET")) stop("SpaCET smoke object failed reload")
  formals_text <- paste(deparse(formals(SpaCET::create.SpaCET.object)), collapse = " ")
  write.table(data.frame(
    check = c("library", "version", "commit", "MUDAN_commit", "constructor", "save_reload", "formals"),
    result = c("PASS", "1.4.0", desc$RemoteSha, mudan_desc$RemoteSha, "PASS", "PASS", formals_text)
  ), file.path(runtime, "interface_check.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA")
  write_status("PASS")
  file.create(file.path(runtime, "_SUCCESS"))
  message("SPACET_EVENT status=PASS version=1.4.0 smoke=PASS")
}, error = function(e) {
  write_status("BLOCKED_DEPENDENCY", conditionMessage(e))
  writeLines(conditionMessage(e), file.path(runtime, "error_message.txt"))
  message("SPACET_EVENT status=BLOCKED_DEPENDENCY error=", conditionMessage(e))
  quit(save = "no", status = 20L)
})
