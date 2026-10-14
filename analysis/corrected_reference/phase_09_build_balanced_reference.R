args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath("/home/zhenhanli/myocardial/scop/GSE214611_spatial", mustWork = TRUE)
.libPaths(c(file.path(root, "Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

input_path <- file.path(root, "objects", "phase_09", "reference_expanded_immune.rds")
output_path <- file.path(root, "objects", "phase_09", "reference_expanded_immune_balanced.rds")
out_dir <- file.path(root, "results", "phase_09", "reference_build")
log_file <- file.path(root, "logs", "phase_09_build_balanced_reference.log")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

stamp <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "|", ..., "\n")
}
file_sha256 <- function(path) {
  sub(" .*", "", system2("sha256sum", shQuote(path), stdout = TRUE))
}

if (!file.exists(input_path)) stop("Missing expanded reference: ", input_path)
stamp("Loading full expanded reference")
full <- readRDS(input_path)
required_meta <- c(
  "sample_id", "timepoint", "P09_audited_state", "P09_expanded_celltype"
)
if (!inherits(full, "Seurat") ||
    !"RNA" %in% Assays(full) ||
    !"counts" %in% Layers(full[["RNA"]]) ||
    !all(required_meta %in% colnames(full[[]]))) {
  stop("Full expanded reference is incompatible")
}

set.seed(214611)
cap_per_large_class <- 3000L
meta <- full[[]]
celltypes <- sort(unique(as.character(meta$P09_expanded_celltype)))
selected_by_class <- lapply(celltypes, function(celltype) {
  cells <- rownames(meta)[as.character(meta$P09_expanded_celltype) == celltype]
  if (length(cells) <= cap_per_large_class) return(cells)
  audited <- as.character(meta[cells, "P09_audited_state"])
  audited[is.na(audited) | !nzchar(audited)] <- "not_applicable"
  strata <- paste(
    as.character(meta[cells, "sample_id"]),
    as.character(meta[cells, "timepoint"]),
    audited,
    sep = "|"
  )
  by_stratum <- split(cells, strata)
  anchors <- vapply(by_stratum, function(x) sample(x, 1L), character(1L))
  remaining <- setdiff(cells, anchors)
  c(anchors, sample(remaining, cap_per_large_class - length(anchors)))
})
names(selected_by_class) <- celltypes
selected_cells <- unlist(selected_by_class, use.names = FALSE)
balanced <- subset(full, cells = selected_cells)

full_counts <- LayerData(full, assay = "RNA", layer = "counts")
balanced_counts <- LayerData(balanced, assay = "RNA", layer = "counts")
profile_rows <- list()
for (celltype in celltypes) {
  full_cells <- rownames(meta)[as.character(meta$P09_expanded_celltype) == celltype]
  balanced_cells <- colnames(balanced)[
    as.character(balanced$P09_expanded_celltype) == celltype
  ]
  full_sum <- Matrix::rowSums(full_counts[, full_cells, drop = FALSE])
  balanced_sum <- Matrix::rowSums(balanced_counts[, balanced_cells, drop = FALSE])
  full_cpm <- 1e6 * full_sum / sum(full_sum)
  balanced_cpm <- 1e6 * balanced_sum / sum(balanced_sum)
  profile_rows[[celltype]] <- data.table(
    celltype = celltype,
    full_cells = length(full_cells),
    balanced_cells = length(balanced_cells),
    pearson_log1p_cpm = cor(log1p(full_cpm), log1p(balanced_cpm), method = "pearson"),
    spearman_cpm = cor(full_cpm, balanced_cpm, method = "spearman"),
    median_abs_log1p_cpm_difference = median(abs(log1p(full_cpm) - log1p(balanced_cpm))),
    q95_abs_log1p_cpm_difference = unname(quantile(
      abs(log1p(full_cpm) - log1p(balanced_cpm)), 0.95
    ))
  )
}
profile_audit <- rbindlist(profile_rows)

balanced@misc$phase_09_balanced_reference <- list(
  source_path = input_path,
  source_sha256 = file_sha256(input_path),
  random_seed = 214611L,
  cap_per_large_class = cap_per_large_class,
  selected_cells_by_class = lapply(selected_by_class, length),
  profile_audit = profile_audit
)
saveRDS(balanced, output_path, compress = FALSE)
readback <- readRDS(output_path)

class_counts <- as.data.table(table(
  celltype = as.character(readback$P09_expanded_celltype)
))
setnames(class_counts, "N", "balanced_cells")
full_class_counts <- as.data.table(table(
  celltype = as.character(full$P09_expanded_celltype)
))
setnames(full_class_counts, "N", "full_cells")
class_counts <- merge(full_class_counts, class_counts, by = "celltype", all = TRUE)
fwrite(
  class_counts,
  file.path(root, "metadata", "phase_09_balanced_reference_class_counts.tsv"),
  sep = "\t"
)
fwrite(
  profile_audit,
  file.path(root, "metadata", "phase_09_balanced_reference_profile_audit.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    reference_variant = c("full_audit", "balanced_primary"),
    path = c(input_path, output_path),
    sha256 = c(file_sha256(input_path), file_sha256(output_path)),
    cells = c(ncol(full), ncol(readback)),
    features = c(nrow(full), nrow(readback)),
    intended_samples = c(
      "Sham_Day3_rep1;Sham_Day3_rep2",
      "MI_1HR_rep1;MI_4HR_rep1;MI_Day3_rep1;MI_Day3_rep2;MI_Day3_rep3;MI_Day7_rep1;MI_Day7_rep2;MI_Day7_rep3;NeedlePass_Day3_rep3;NeedlePass_Day3_rep4"
    )
  ),
  file.path(root, "metadata", "phase_09_expanded_reference_variants.tsv"),
  sep = "\t"
)

validation <- data.table(
  check = c(
    "readback_seurat", "features_3000", "all_11_classes",
    "rare_classes_retained_complete", "large_classes_capped_3000",
    "minimum_profile_pearson_gt_0_98", "minimum_profile_spearman_gt_0_98",
    "source_reference_unchanged"
  ),
  status = c(
    ifelse(inherits(readback, "Seurat"), "PASS", "FAIL"),
    ifelse(nrow(readback) == 3000L, "PASS", "FAIL"),
    ifelse(length(unique(as.character(readback$P09_expanded_celltype))) == 11L, "PASS", "FAIL"),
    ifelse(all(class_counts[full_cells <= cap_per_large_class, full_cells == balanced_cells]), "PASS", "FAIL"),
    ifelse(all(class_counts[full_cells > cap_per_large_class, balanced_cells == cap_per_large_class]), "PASS", "FAIL"),
    ifelse(min(profile_audit$pearson_log1p_cpm) > 0.98, "PASS", "FAIL"),
    ifelse(min(profile_audit$spearman_cpm) > 0.98, "PASS", "FAIL"),
    ifelse(file_sha256(input_path) == balanced@misc$phase_09_balanced_reference$source_sha256, "PASS", "FAIL")
  ),
  value = c(
    inherits(readback, "Seurat"), nrow(readback),
    length(unique(as.character(readback$P09_expanded_celltype))),
    sum(class_counts$full_cells <= cap_per_large_class),
    sum(class_counts$full_cells > cap_per_large_class),
    min(profile_audit$pearson_log1p_cpm),
    min(profile_audit$spearman_cpm),
    file_sha256(input_path)
  )
)
fwrite(
  validation,
  file.path(root, "metadata", "phase_09_balanced_reference_validation.tsv"),
  sep = "\t"
)
stamp(
  "Balanced reference cells=", ncol(readback),
  " minimum Pearson=", round(min(profile_audit$pearson_log1p_cpm), 4),
  " minimum Spearman=", round(min(profile_audit$spearman_cpm), 4)
)
stamp("Validation:", sum(validation$status == "PASS"), "/", nrow(validation), "PASS")
if (any(validation$status == "FAIL")) quit(save = "no", status = 2L)
