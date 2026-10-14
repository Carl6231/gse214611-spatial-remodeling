#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else {
  "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
}
.libPaths(c(file.path(root, "Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(scop)
})
if (!"SpatialSpotPlot" %in% getNamespaceExports("scop")) {
  stop("scop::SpatialSpotPlot is unavailable")
}

for (path in c(
  file.path(root, "results", "phase_09", "macrophage_mapping", "summary"),
  file.path(root, "results", "phase_09", "macrophage_mapping", "visual_qc"),
  file.path(root, "metadata"),
  file.path(root, "logs")
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
log_path <- file.path(root, "logs", "phase_09_aggregate_macrophage.log")
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

samples <- c(
  "Sham_Day3_rep1", "Sham_Day3_rep2", "MI_1HR_rep1", "MI_4HR_rep1",
  "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
  "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3",
  "NeedlePass_Day3_rep3", "NeedlePass_Day3_rep4"
)
states <- c(
  "Resident_Folr2_macrophage",
  "Trem2_Spp1_macrophage",
  "MHCII_macrophage",
  "Cycling_macrophage",
  "Stress_response_macrophage"
)
state_palette <- c(
  Resident_Folr2_macrophage = "#2F6FB0",
  Trem2_Spp1_macrophage = "#D9822B",
  MHCII_macrophage = "#7B2CBF",
  Cycling_macrophage = "#2A8C82",
  Stress_response_macrophage = "#C81D25",
  Uncertain = "#8C8C8C",
  Not_macrophage_candidate = "#D0D0D0"
)

status_files <- file.path(
  root, "metadata", "phase_09_macrophage_sample_status",
  paste0(samples, ".tsv")
)
if (!all(file.exists(status_files))) stop("Macrophage status files incomplete")
status <- do.call(rbind, lapply(status_files, read.delim, check.names = FALSE))
if (nrow(status) != 12L || any(status$status != "PASS")) {
  stop("Refined macrophage mapping is not 12/12 PASS")
}
write_tsv(
  status,
  file.path(root, "metadata", "phase_09_macrophage_sample_status.tsv")
)

mapping_files <- file.path(
  root, "results", "phase_09", "macrophage_mapping", samples,
  "spot_mapping.tsv"
)
if (!all(file.exists(mapping_files))) stop("Macrophage spot tables incomplete")
mapping <- do.call(rbind, lapply(mapping_files, function(path) {
  read.delim(path, check.names = FALSE)
}))
if (anyDuplicated(paste(mapping$sample_key, mapping$barcode, sep = "::"))) {
  stop("Duplicated sample-barcode mapping rows")
}
write_tsv_gz(
  mapping,
  file.path(root, "metadata", "phase_09_macrophage_spot_mapping.tsv.gz")
)

row_sum_rows <- list()
sample_summary_rows <- list()
for (sample in samples) {
  current <- mapping[mapping$sample_key == sample, , drop = FALSE]
  rctd <- as.matrix(current[, paste0("RCTD_", states), drop = FALSE])
  dwls <- as.matrix(current[, paste0("DWLS_", states), drop = FALSE])
  colnames(rctd) <- states
  colnames(dwls) <- states
  row_sum_rows[[sample]] <- data.frame(
    sample_key = sample,
    mapped_spots = nrow(current),
    rctd_min_row_sum = min(rowSums(rctd)),
    rctd_max_row_sum = max(rowSums(rctd)),
    dwls_min_row_sum = min(rowSums(dwls)),
    dwls_max_row_sum = max(rowSums(dwls)),
    rctd_nonnegative_finite = all(is.finite(rctd) & rctd >= 0),
    dwls_nonnegative_finite = all(is.finite(dwls) & dwls >= 0),
    stringsAsFactors = FALSE
  )
  for (method in c("RCTD", "DWLS")) {
    matrix <- if (method == "RCTD") rctd else dwls
    for (state in states) {
      sample_summary_rows[[length(sample_summary_rows) + 1L]] <- data.frame(
        sample_key = sample,
        method = method,
        state = state,
        mapped_spots = nrow(current),
        mean_state_proportion = mean(matrix[, state]),
        median_state_proportion = median(matrix[, state]),
        stringsAsFactors = FALSE
      )
    }
  }
}
row_sum_audit <- do.call(rbind, row_sum_rows)
sample_summary <- do.call(rbind, sample_summary_rows)
write_tsv(
  row_sum_audit,
  file.path(root, "metadata", "phase_09_macrophage_row_sum_audit.tsv")
)
write_tsv(
  sample_summary,
  file.path(root, "metadata", "phase_09_macrophage_state_summary.tsv")
)

dominant_summary <- as.data.frame(table(
  sample_key = mapping$sample_key,
  dominant_state = mapping$P09M_dominant_state
), stringsAsFactors = FALSE)
dominant_summary <- dominant_summary[dominant_summary$Freq > 0L, ]
total_by_sample <- aggregate(Freq ~ sample_key, dominant_summary, sum)
names(total_by_sample)[2L] <- "mapped_spots"
dominant_summary <- merge(
  dominant_summary, total_by_sample, by = "sample_key", all.x = TRUE
)
dominant_summary$fraction <- dominant_summary$Freq / dominant_summary$mapped_spots
write_tsv(
  dominant_summary,
  file.path(root, "metadata", "phase_09_macrophage_dominant_summary.tsv")
)

reference_composition <- read.delim(
  file.path(root, "metadata", "phase_09_refined_macrophage_composition.tsv"),
  check.names = FALSE
)
state_support <- do.call(rbind, lapply(split(
  reference_composition, reference_composition$state
), function(x) {
  data.frame(
    state = x$state[[1L]],
    cells = unique(x$state_cells)[[1L]],
    samples_detected = length(unique(x$sample_id)),
    timepoints_detected = length(unique(x$timepoint)),
    max_sample_fraction = max(x$fraction_within_state),
    stringsAsFactors = FALSE
  )
}))
rownames(state_support) <- NULL
write_tsv(
  state_support,
  file.path(root, "metadata", "phase_09_macrophage_state_support.tsv")
)

state_maps <- character(length(samples))
confidence_maps <- character(length(samples))
object_audit_rows <- list()
for (i in seq_along(samples)) {
  sample <- samples[[i]]
  object_path <- file.path(
    root, "objects", "phase_09", "macrophage_mapping",
    paste0(sample, "_phase09_macrophage.rds")
  )
  object <- readRDS(object_path)
  if (!inherits(object, "Seurat") || length(Images(object)) != 1L) {
    stop("Invalid mapped object: ", sample)
  }
  image_name <- Images(object)[[1L]]
  state_values <- as.character(object$P09M_dominant_state)
  state_values[!nzchar(state_values) | is.na(state_values)] <-
    "Not_macrophage_candidate"
  state_maps[[i]] <- file.path(
    root, "results", "phase_09", "macrophage_mapping",
    sample, "refined_state_map.png"
  )
  confidence_maps[[i]] <- file.path(
    root, "results", "phase_09", "macrophage_mapping",
    sample, "refined_confidence_map.png"
  )
  p_state <- SpatialSpotPlot(
    object,
    values = setNames(state_values, colnames(object)),
    plot_type = "point",
    image = image_name,
    crop = FALSE,
    coord.cols = c("x", "y"),
    overlay_image = TRUE,
    image.alpha = 0.85,
    pt.alpha = 0.9,
    palcolor = state_palette[unique(state_values)],
    legend.title = "Refined state"
  ) + ggtitle(
    paste0(sample, ": corrected macrophage states"),
    subtitle = "Only corrected Macrophage >=0.10 spots were mapped"
  )
  ggsave(
    state_maps[[i]], p_state,
    width = 9, height = 8, dpi = 260, bg = "white"
  )

  confidence <- as.numeric(object$P09M_method_cosine)
  confidence[!as.logical(object$P09M_candidate)] <- NA_real_
  p_confidence <- SpatialSpotPlot(
    object,
    values = setNames(confidence, colnames(object)),
    plot_type = "point",
    image = image_name,
    crop = FALSE,
    coord.cols = c("x", "y"),
    overlay_image = TRUE,
    image.alpha = 0.85,
    pt.alpha = 0.9,
    palcolor = c("#440154", "#21908C", "#FDE725"),
    legend.title = "RCTD-DWLS cosine"
  ) + ggtitle(paste0(sample, ": mapping agreement"))
  ggsave(
    confidence_maps[[i]], p_confidence,
    width = 9, height = 8, dpi = 260, bg = "white"
  )
  object_audit_rows[[sample]] <- data.frame(
    sample_key = sample,
    spots = ncol(object),
    candidate_spots = sum(object$P09M_candidate),
    images = length(Images(object)),
    exact_mapping_rows =
      sum(object$P09M_candidate) ==
        sum(mapping$sample_key == sample),
    state_map_exported = file.exists(state_maps[[i]]),
    confidence_map_exported = file.exists(confidence_maps[[i]]),
    stringsAsFactors = FALSE
  )
  rm(object)
  gc()
}
object_audit <- do.call(rbind, object_audit_rows)
write_tsv(
  object_audit,
  file.path(root, "metadata", "phase_09_macrophage_object_audit.tsv")
)

montage_status <- c(
  system2("montage", c(
    shQuote(state_maps), "-tile", "3x4", "-geometry", "650x560+8+8",
    shQuote(file.path(
      root, "results", "phase_09", "macrophage_mapping", "visual_qc",
      "phase_09_refined_state_contact_sheet.jpg"
    ))
  )),
  system2("montage", c(
    shQuote(confidence_maps), "-tile", "3x4", "-geometry", "650x560+8+8",
    shQuote(file.path(
      root, "results", "phase_09", "macrophage_mapping", "visual_qc",
      "phase_09_refined_confidence_contact_sheet.jpg"
    ))
  ))
)
if (any(montage_status != 0L)) stop("Macrophage contact sheet failed")

validation <- data.frame(
  check = c(
    "sample_status_12_pass",
    "mapping_rows_match_status",
    "unique_sample_barcodes",
    "rctd_rows_sum_one",
    "dwls_rows_sum_one",
    "finite_nonnegative_proportions",
    "low_confidence_is_uncertain",
    "states_equal_5",
    "state_min_reference_cells_25",
    "state_multisample_support",
    "object_readback_12",
    "spatial_maps_12",
    "confidence_maps_12",
    "contact_sheets_2"
  ),
  value = c(
    nrow(status) == 12L && all(status$status == "PASS"),
    nrow(mapping) == sum(status$mapped_spots),
    !anyDuplicated(paste(mapping$sample_key, mapping$barcode)),
    all(abs(row_sum_audit$rctd_min_row_sum - 1) < 1e-8) &&
      all(abs(row_sum_audit$rctd_max_row_sum - 1) < 1e-8),
    all(abs(row_sum_audit$dwls_min_row_sum - 1) < 1e-8) &&
      all(abs(row_sum_audit$dwls_max_row_sum - 1) < 1e-8),
    all(row_sum_audit$rctd_nonnegative_finite) &&
      all(row_sum_audit$dwls_nonnegative_finite),
    all(
      mapping$P09M_confidence != "Low" |
        mapping$P09M_dominant_state == "Uncertain"
    ),
    setequal(
      sub("^RCTD_", "", grep("^RCTD_", colnames(mapping), value = TRUE)),
      states
    ),
    min(state_support$cells) >= 25L,
    min(state_support$samples_detected) >= 2L,
    nrow(object_audit) == 12L &&
      all(object_audit$exact_mapping_rows) &&
      all(object_audit$images == 1L),
    sum(file.exists(state_maps)) == 12L,
    sum(file.exists(confidence_maps)) == 12L,
    all(file.exists(file.path(
      root, "results", "phase_09", "macrophage_mapping", "visual_qc",
      c(
        "phase_09_refined_state_contact_sheet.jpg",
        "phase_09_refined_confidence_contact_sheet.jpg"
      )
    )))
  ),
  status = "",
  stringsAsFactors = FALSE
)
validation$status <- ifelse(as.logical(validation$value), "PASS", "FAIL")
write_tsv(
  validation,
  file.path(root, "metadata", "phase_09_macrophage_aggregate_validation.tsv")
)
log_msg(
  "Macrophage aggregate complete: mapped_spots=", nrow(mapping),
  " validation=", paste(table(validation$status), collapse = ",")
)
if (any(validation$status != "PASS")) {
  stop("Macrophage aggregate validation failed")
}
