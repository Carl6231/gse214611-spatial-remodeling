#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else {
  "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
}
.libPaths(c(file.path(root, "Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
})

dir.create(file.path(root, "metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results", "phase_09", "expanded_deconvolution", "summary"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "logs"), recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(root, "logs", "phase_09_aggregate_deconvolution.log")
log_connection <- file(log_path, open = "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  while (sink.number(type = "message") > 2L) sink(type = "message")
  while (sink.number() > 0L) sink()
  close(log_connection)
}, add = TRUE)

stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
log_msg <- function(...) cat(stamp(), "|", ..., "\n")
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
exact_permutation_test <- function(x, y) {
  observed <- mean(x) - mean(y)
  combined <- c(x, y)
  n_x <- length(x)
  combinations <- combn(seq_along(combined), n_x)
  null <- apply(combinations, 2L, function(index) {
    mean(combined[index]) - mean(combined[-index])
  })
  list(
    effect = observed,
    p_two_sided = mean(abs(null) >= abs(observed) - 1e-12),
    null_min = min(null),
    null_max = max(null),
    permutations = length(null)
  )
}

samples <- c(
  "Sham_Day3_rep1", "Sham_Day3_rep2", "MI_1HR_rep1", "MI_4HR_rep1",
  "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
  "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3",
  "NeedlePass_Day3_rep3", "NeedlePass_Day3_rep4"
)
celltypes <- c(
  "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage",
  "Monocyte", "Neutrophil", "Dendritic", "B_cell", "T_NK",
  "Pericyte", "Schwann_cell"
)

status_files <- file.path(
  root, "metadata", "phase_09_expanded_sample_status",
  paste0(samples, ".tsv")
)
if (!all(file.exists(status_files))) {
  stop("Expanded sample status files are incomplete")
}
status <- do.call(rbind, lapply(status_files, read.delim, check.names = FALSE))
if (nrow(status) != 12L || any(status$status != "PASS")) {
  stop("Expanded deconvolution is not 12/12 PASS")
}
write_tsv(status, file.path(root, "metadata", "phase_09_expanded_sample_status.tsv"))

spot_files <- file.path(
  root, "results", "phase_09", "expanded_deconvolution", samples,
  "spot_proportions.tsv"
)
if (!all(file.exists(spot_files))) stop("Spot proportion tables are incomplete")
spots <- do.call(rbind, lapply(spot_files, function(path) {
  read.delim(path, check.names = FALSE)
}))
if (anyDuplicated(paste(spots$sample_key, spots$barcode, sep = "::"))) {
  stop("Duplicated sample-barcode rows")
}
write_tsv_gz(
  spots,
  file.path(root, "metadata", "phase_09_expanded_spot_proportions.tsv.gz")
)

row_sum_audit <- do.call(rbind, lapply(samples, function(sample) {
  current <- spots[spots$sample_key == sample, , drop = FALSE]
  rctd <- as.matrix(current[, paste0("RCTD_", celltypes), drop = FALSE])
  dwls <- as.matrix(
    current[, paste0("SpatialDWLS_", celltypes), drop = FALSE]
  )
  data.frame(
    sample_key = sample,
    spots = nrow(current),
    rctd_min_row_sum = min(rowSums(rctd)),
    rctd_max_row_sum = max(rowSums(rctd)),
    dwls_min_row_sum = min(rowSums(dwls)),
    dwls_max_row_sum = max(rowSums(dwls)),
    rctd_all_finite_nonnegative = all(is.finite(rctd) & rctd >= 0),
    dwls_all_finite_nonnegative = all(is.finite(dwls) & dwls >= 0),
    stringsAsFactors = FALSE
  )
}))
write_tsv(
  row_sum_audit,
  file.path(root, "metadata", "phase_09_expanded_row_sum_audit.tsv")
)

summary_rows <- list()
for (sample in samples) {
  current <- spots[spots$sample_key == sample, , drop = FALSE]
  for (method in c("RCTD", "SpatialDWLS")) {
    for (celltype in celltypes) {
      column <- paste0(method, "_", celltype)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        sample_key = sample,
        GSM = current$GSM[[1L]],
        sample_id = current$sample_id[[1L]],
        condition = current$condition[[1L]],
        timepoint = current$timepoint[[1L]],
        replicate = current$replicate[[1L]],
        analysis_role = current$analysis_role[[1L]],
        include_primary_MI_analysis =
          current$include_primary_MI_analysis[[1L]],
        method = method,
        component = celltype,
        mean_proportion = mean(current[[column]]),
        median_proportion = median(current[[column]]),
        q25 = unname(quantile(current[[column]], 0.25)),
        q75 = unname(quantile(current[[column]], 0.75)),
        stringsAsFactors = FALSE
      )
    }
  }
}
sample_summary <- do.call(rbind, summary_rows)
write_tsv(
  sample_summary,
  file.path(root, "metadata", "phase_09_expanded_sample_summary.tsv")
)

reallocation <- do.call(rbind, lapply(split(spots, spots$sample_key), function(x) {
  data.frame(
    sample_key = x$sample_key[[1L]],
    condition = x$condition[[1L]],
    timepoint = x$timepoint[[1L]],
    replicate = x$replicate[[1L]],
    analysis_role = x$analysis_role[[1L]],
    include_primary_MI_analysis = x$include_primary_MI_analysis[[1L]],
    spots = nrow(x),
    old_P4_macrophage_mean = mean(x$old_P4_macrophage),
    corrected_macrophage_mean = mean(x$RCTD_Macrophage),
    corrected_monocyte_mean = mean(x$RCTD_Monocyte),
    corrected_neutrophil_mean = mean(x$RCTD_Neutrophil),
    corrected_dendritic_mean = mean(x$RCTD_Dendritic),
    corrected_B_cell_mean = mean(x$RCTD_B_cell),
    corrected_T_NK_mean = mean(x$RCTD_T_NK),
    corrected_total_myeloid_mean = mean(
      x$RCTD_Macrophage + x$RCTD_Monocyte +
        x$RCTD_Neutrophil + x$RCTD_Dendritic
    ),
    macrophage_delta_new_minus_old = mean(
      x$RCTD_Macrophage - x$old_P4_macrophage
    ),
    old_new_macrophage_spearman = suppressWarnings(cor(
      x$old_P4_macrophage, x$RCTD_Macrophage,
      method = "spearman", use = "pairwise.complete.obs"
    )),
    stringsAsFactors = FALSE
  )
}))
rownames(reallocation) <- NULL
write_tsv(
  reallocation,
  file.path(root, "metadata", "phase_09_macrophage_reallocation_by_slice.tsv")
)

contrasts <- list(
  Day3_MI_vs_Sham = list(
    numerator = c("MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3"),
    denominator = c("Sham_Day3_rep1", "Sham_Day3_rep2")
  ),
  Day7_MI_vs_Sham = list(
    numerator = c("MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3"),
    denominator = c("Sham_Day3_rep1", "Sham_Day3_rep2")
  ),
  Day7_MI_vs_Day3_MI = list(
    numerator = c("MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3"),
    denominator = c("MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3")
  )
)
comparison_rows <- list()
for (contrast_name in names(contrasts)) {
  contrast <- contrasts[[contrast_name]]
  for (method in c("RCTD", "SpatialDWLS")) {
    for (component in celltypes) {
      values <- sample_summary[
        sample_summary$method == method &
          sample_summary$component == component, ,
        drop = FALSE
      ]
      numerator <- values$mean_proportion[
        match(contrast$numerator, values$sample_key)
      ]
      denominator <- values$mean_proportion[
        match(contrast$denominator, values$sample_key)
      ]
      test <- exact_permutation_test(numerator, denominator)
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        contrast = contrast_name,
        method = method,
        component = component,
        numerator_n = length(numerator),
        denominator_n = length(denominator),
        numerator_mean = mean(numerator),
        denominator_mean = mean(denominator),
        mean_difference = test$effect,
        log2_fold_change = log2(
          (mean(numerator) + 1e-6) / (mean(denominator) + 1e-6)
        ),
        exact_p_two_sided = test$p_two_sided,
        permutations = test$permutations,
        stringsAsFactors = FALSE
      )
    }
  }
}
comparisons <- do.call(rbind, comparison_rows)
write_tsv(
  comparisons,
  file.path(root, "metadata", "phase_09_expanded_slice_comparisons.tsv")
)

direction_table <- reshape(
  comparisons[, c(
    "contrast", "method", "component", "mean_difference",
    "log2_fold_change", "exact_p_two_sided"
  )],
  idvar = c("contrast", "component"),
  timevar = "method",
  direction = "wide"
)
direction_table$direction_concordant <- sign(
  direction_table$mean_difference.RCTD
) == sign(direction_table$mean_difference.SpatialDWLS)
write_tsv(
  direction_table,
  file.path(root, "metadata", "phase_09_expanded_method_sensitivity.tsv")
)

plot_data <- reallocation[
  reallocation$sample_key %in% c(
    "Sham_Day3_rep1", "Sham_Day3_rep2",
    "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
    "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3"
  ), ,
]
plot_data$group <- ifelse(
  grepl("^Sham", plot_data$sample_key), "Sham Day3",
  ifelse(grepl("Day3", plot_data$sample_key), "MI Day3", "MI Day7")
)
plot_long <- reshape(
  plot_data[, c(
    "sample_key", "group", "old_P4_macrophage_mean",
    "corrected_macrophage_mean", "corrected_monocyte_mean",
    "corrected_neutrophil_mean", "corrected_dendritic_mean"
  )],
  varying = c(
    "old_P4_macrophage_mean", "corrected_macrophage_mean",
    "corrected_monocyte_mean", "corrected_neutrophil_mean",
    "corrected_dendritic_mean"
  ),
  v.names = "mean_proportion",
  timevar = "component",
  times = c(
    "Old six-class Macrophage", "Corrected Macrophage", "Monocyte",
    "Neutrophil", "Dendritic"
  ),
  direction = "long"
)
plot_long$group <- factor(
  plot_long$group, levels = c("Sham Day3", "MI Day3", "MI Day7")
)
plot_long$component <- factor(
  plot_long$component,
  levels = c(
    "Old six-class Macrophage", "Corrected Macrophage",
    "Monocyte", "Neutrophil", "Dendritic"
  )
)
p <- ggplot(
  plot_long,
  aes(group, mean_proportion, color = component, group = component)
) +
  geom_line(aes(group = interaction(component, sample_key)),
            linewidth = 0.35, alpha = 0.35) +
  geom_point(position = position_jitter(width = 0.07), size = 2.3) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 4) +
  scale_color_manual(values = c(
    "Old six-class Macrophage" = "#8C8C8C",
    "Corrected Macrophage" = "#CC79A7",
    "Monocyte" = "#C81D25",
    "Neutrophil" = "#E15759",
    "Dendritic" = "#7B2CBF"
  )) +
  labs(
    title = "Expanded-reference immune reallocation",
    subtitle = "Dots are physical slices; no spot-level inference",
    x = NULL, y = "Mean slice-level RCTD proportion", color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(
  file.path(
    root, "results", "phase_09", "expanded_deconvolution", "summary",
    "phase_09_reference_reallocation.png"
  ),
  p, width = 9, height = 6, dpi = 320, bg = "white"
)
ggsave(
  file.path(
    root, "results", "phase_09", "expanded_deconvolution", "summary",
    "phase_09_reference_reallocation.pdf"
  ),
  p, width = 9, height = 6, device = cairo_pdf
)

immune_maps <- file.path(
  root, "results", "phase_09", "expanded_deconvolution", samples,
  "expanded_immune_rctd.png"
)
dominant_maps <- file.path(
  root, "results", "phase_09", "expanded_deconvolution", samples,
  "expanded_consensus_dominant.png"
)
if (!all(file.exists(c(immune_maps, dominant_maps)))) {
  stop("Per-slice native scop plots are incomplete")
}
montage_status <- c(
  system2("montage", c(
    shQuote(immune_maps), "-tile", "3x4", "-geometry", "900x620+8+8",
    shQuote(file.path(
      root, "results", "phase_09", "expanded_deconvolution", "summary",
      "phase_09_expanded_immune_contact_sheet.jpg"
    ))
  )),
  system2("montage", c(
    shQuote(dominant_maps), "-tile", "3x4", "-geometry", "650x560+8+8",
    shQuote(file.path(
      root, "results", "phase_09", "expanded_deconvolution", "summary",
      "phase_09_expanded_dominant_contact_sheet.jpg"
    ))
  ))
)
if (any(montage_status != 0L)) stop("Contact sheet creation failed")

validation <- data.frame(
  check = c(
    "sample_status_12_pass",
    "spot_rows_25548",
    "unique_sample_barcodes",
    "rctd_row_sums_one",
    "dwls_row_sums_one",
    "finite_nonnegative_proportions",
    "phase4_macrophage_present",
    "slice_level_comparisons_complete",
    "method_sensitivity_complete",
    "native_scop_immune_maps_12",
    "native_scop_dominant_maps_12",
    "contact_sheets_exported"
  ),
  value = c(
    nrow(status) == 12L && all(status$status == "PASS"),
    nrow(spots) == 25548L,
    !anyDuplicated(paste(spots$sample_key, spots$barcode)),
    all(abs(row_sum_audit$rctd_min_row_sum - 1) < 1e-8) &&
      all(abs(row_sum_audit$rctd_max_row_sum - 1) < 1e-8),
    all(abs(row_sum_audit$dwls_min_row_sum - 1) < 1e-8) &&
      all(abs(row_sum_audit$dwls_max_row_sum - 1) < 1e-8),
    all(row_sum_audit$rctd_all_finite_nonnegative) &&
      all(row_sum_audit$dwls_all_finite_nonnegative),
    all(is.finite(spots$old_P4_macrophage)),
    nrow(comparisons) == 3L * 2L * length(celltypes),
    nrow(direction_table) == 3L * length(celltypes),
    sum(file.exists(immune_maps)) == 12L,
    sum(file.exists(dominant_maps)) == 12L,
    all(file.exists(file.path(
      root, "results", "phase_09", "expanded_deconvolution", "summary",
      c(
        "phase_09_expanded_immune_contact_sheet.jpg",
        "phase_09_expanded_dominant_contact_sheet.jpg"
      )
    )))
  ),
  status = "",
  stringsAsFactors = FALSE
)
validation$status <- ifelse(as.logical(validation$value), "PASS", "FAIL")
write_tsv(
  validation,
  file.path(root, "metadata", "phase_09_expanded_aggregate_validation.tsv")
)
log_msg(
  "Aggregate complete: spots=", nrow(spots),
  " validation=", paste(table(validation$status), collapse = ",")
)
if (any(validation$status != "PASS")) stop("Expanded aggregate validation failed")
