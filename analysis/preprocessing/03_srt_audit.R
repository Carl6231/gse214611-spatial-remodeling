#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
dir.create(file.path(root, "logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(root, "logs", "srt_audit.txt")
sink(log_file, append = FALSE, split = TRUE)
on.exit(sink(), add = TRUE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

srt_path <- file.path(root, "reference", "srt.rds")
requested_metadata <- c("celltype", "sample_id", "timepoint")

log_msg("read-only Seurat audit started")
log_msg("path:", srt_path)
log_msg("symlink:", Sys.readlink(srt_path))
log_msg("exists:", file.exists(srt_path))
if (!file.exists(srt_path)) {
  log_msg("readRDS_status: FAIL(file_missing)")
  quit(save = "no", status = 1)
}

info <- file.info(srt_path)
log_msg("size_bytes:", info$size)
log_msg("mtime:", format(info$mtime, "%Y-%m-%d %H:%M:%S %z"))
log_msg("R:", R.version.string)
log_msg("Seurat_version:", as.character(packageVersion("Seurat")))
log_msg("SeuratObject_version:", as.character(packageVersion("SeuratObject")))

srt <- tryCatch(
  readRDS(srt_path),
  error = function(e) e
)
if (inherits(srt, "error")) {
  log_msg("readRDS_status: FAIL")
  log_msg("readRDS_error:", conditionMessage(srt))
  quit(save = "no", status = 1)
}

log_msg("readRDS_status: PASS")
log_msg("class:", paste(class(srt), collapse = ","))
log_msg("is_Seurat:", inherits(srt, "Seurat"))
if (!inherits(srt, "Seurat")) {
  log_msg("object is not a Seurat object; stopping after class audit")
  quit(save = "no", status = 1)
}

assays <- SeuratObject::Assays(srt)
log_msg("Assays:", paste(assays, collapse = ","))
for (assay_name in assays) {
  layers <- tryCatch(SeuratObject::Layers(srt[[assay_name]]), error = function(e) paste0("ERROR:", conditionMessage(e)))
  log_msg(paste0("Layers[", assay_name, "]:"), paste(layers, collapse = ","))
}

metadata <- srt[[]]
log_msg("metadata_columns_n:", ncol(metadata))
log_msg("metadata_columns:", paste(colnames(metadata), collapse = ","))
for (column_name in requested_metadata) {
  log_msg(paste0("metadata_has_", column_name, ":"), column_name %in% colnames(metadata))
}

if ("celltype" %in% colnames(metadata)) {
  cat("\n===== table(celltype, useNA = 'ifany') =====\n")
  print(table(metadata$celltype, useNA = "ifany"))
} else {
  cat("\n===== table(celltype) unavailable: column missing =====\n")
}

if (all(c("sample_id", "timepoint") %in% colnames(metadata))) {
  cat("\n===== table(sample_id, timepoint, useNA = 'ifany') =====\n")
  print(table(metadata$sample_id, metadata$timepoint, useNA = "ifany"))
} else {
  cat("\n===== table(sample_id, timepoint) unavailable: one or both columns missing =====\n")
}

log_msg("cells:", ncol(srt))
log_msg("features_default_assay:", nrow(srt))
log_msg("read-only Seurat audit finished")
