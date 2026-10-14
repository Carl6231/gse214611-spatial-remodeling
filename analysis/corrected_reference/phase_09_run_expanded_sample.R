#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else {
  "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
}
if (length(args) < 2L || !nzchar(args[[2L]])) {
  stop("Usage: phase_09_run_expanded_sample.R <project_root> <sample_key> [primary|balanced_shadow]")
}
sample_key <- args[[2L]]
run_mode <- if (length(args) >= 3L) args[[3L]] else "primary"
if (!run_mode %in% c("primary", "balanced_shadow")) {
  stop("Unsupported run_mode: ", run_mode)
}
.libPaths(c(file.path(root, "Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scop)
})

required_scop <- c(
  "RunRCTD", "RunSpatialDWLS", "SpatialDeconvolutionPlot",
  "SpatialSpotPlot", "GetSpatialResult"
)
if (!all(required_scop %in% getNamespaceExports("scop"))) {
  stop("Required scop exports are unavailable")
}

reference_variant <- if (run_mode == "balanced_shadow") {
  "balanced_shadow"
} else {
  "balanced_primary"
}
reference_path <- file.path(
  root, "objects", "phase_09",
  if (reference_variant == "full_audit") {
    "reference_expanded_immune.rds"
  } else {
    "reference_expanded_immune_balanced.rds"
  }
)
feature_path <- file.path(
  root, "metadata", "phase_09_expanded_feature_selection.tsv"
)
registry_path <- file.path(root, "metadata", "spatial_qc_summary.tsv")
input_path <- file.path(
  root, "objects", "deconvolution", "samples",
  paste0(sample_key, "_phase04.rds")
)
result_dir <- if (run_mode == "balanced_shadow") {
  file.path(
    root, "results", "phase_09", "expanded_deconvolution",
    "balanced_shadow", sample_key
  )
} else {
  file.path(root, "results", "phase_09", "expanded_deconvolution", sample_key)
}
object_dir <- if (run_mode == "balanced_shadow") {
  file.path(root, "objects", "phase_09", "deconvolution_balanced_shadow")
} else {
  file.path(root, "objects", "phase_09", "deconvolution")
}
status_dir <- if (run_mode == "balanced_shadow") {
  file.path(root, "metadata", "phase_09_balanced_shadow_status")
} else {
  file.path(root, "metadata", "phase_09_expanded_sample_status")
}
log_dir <- if (run_mode == "balanced_shadow") {
  file.path(root, "logs", "phase_09_balanced_shadow_samples")
} else {
  file.path(root, "logs", "phase_09_expanded_samples")
}
output_object <- file.path(
  object_dir,
  if (run_mode == "balanced_shadow") {
    paste0(sample_key, "_phase09_expanded_balanced_shadow.rds")
  } else {
    paste0(sample_key, "_phase09_expanded.rds")
  }
)
sample_status_path <- file.path(status_dir, paste0(sample_key, ".tsv"))
log_path <- file.path(log_dir, paste0(sample_key, ".log"))
for (path in c(result_dir, object_dir, status_dir, log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

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
sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  sub(" .*", "", output)
}
normalize_proportions <- function(
    x, spots, celltypes, label, fallback = NULL,
    force_fallback_spots = character()) {
  x <- as.matrix(x)
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(label, " proportions lack dimnames")
  }
  if (all(spots %in% colnames(x)) && !all(spots %in% rownames(x))) {
    x <- t(x)
  }
  if (!all(spots %in% rownames(x))) {
    stop(label, " proportions do not cover all spots")
  }
  missing_types <- setdiff(celltypes, colnames(x))
  if (length(missing_types)) {
    stop(label, " proportions miss: ", paste(missing_types, collapse = ", "))
  }
  x <- x[spots, celltypes, drop = FALSE]
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  fallback_rows <- union(
    rownames(x)[!is.finite(totals) | totals <= 0],
    intersect(force_fallback_spots, rownames(x))
  )
  if (length(fallback_rows)) {
    if (is.null(fallback) ||
        !all(fallback_rows %in% rownames(fallback)) ||
        !all(celltypes %in% colnames(fallback))) {
      stop(label, " fallback does not cover zero-information spots")
    }
    replacement <- as.matrix(
      fallback[fallback_rows, celltypes, drop = FALSE]
    )
    replacement[!is.finite(replacement) | replacement < 0] <- 0
    replacement_totals <- rowSums(replacement)
    if (any(replacement_totals <= 0)) {
      stop(label, " fallback contains zero-sum rows")
    }
    x[fallback_rows, ] <- replacement / replacement_totals
    totals <- rowSums(x)
  }
  if (any(!is.finite(totals) | totals <= 0)) {
    stop(label, " has unresolved zero-sum rows")
  }
  x <- x / totals
  attr(x, "fallback_rows") <- fallback_rows
  x
}
build_phase4_fallback <- function(object, method, spots, celltypes) {
  broad_types <- c(
    "Cardiomyocyte", "Endothelial", "Fibroblast",
    "Macrophage", "Pericyte", "Schwann_cell"
  )
  fields <- paste0("P4_", method, "_prop_", broad_types)
  metadata <- object[[]]
  if (!all(spots %in% rownames(metadata)) ||
      !all(fields %in% colnames(metadata))) {
    stop("Phase 4 ", method, " fallback fields are incomplete")
  }
  fallback <- matrix(
    0, nrow = length(spots), ncol = length(celltypes),
    dimnames = list(spots, celltypes)
  )
  fallback[, broad_types] <- as.matrix(
    metadata[spots, fields, drop = FALSE]
  )
  fallback[!is.finite(fallback) | fallback < 0] <- 0
  totals <- rowSums(fallback)
  if (any(totals <= 0)) {
    stop("Phase 4 ", method, " fallback contains zero-sum rows")
  }
  fallback / totals
}
dominant_type <- function(x) colnames(x)[max.col(x, ties.method = "first")]
top_margin <- function(x) apply(x, 1L, function(z) {
  ordered <- sort(z, decreasing = TRUE)
  ordered[[1L]] - ordered[[2L]]
})
cosine_rows <- function(x, y) {
  denominator <- sqrt(rowSums(x^2) * rowSums(y^2))
  value <- rowSums(x * y) / denominator
  value[!is.finite(value)] <- 0
  pmax(0, pmin(1, value))
}

status_lines <- readLines(file.path(root, "docs", "STATUS.md"), warn = FALSE)
phase9_start <- match("## Phase 9", status_lines)
if (is.na(phase9_start) ||
    !any(status_lines[phase9_start:length(status_lines)] == "- status: RUNNING")) {
  stop("Phase 9 is not RUNNING")
}
for (path in c(reference_path, feature_path, registry_path, input_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

registry <- read.delim(registry_path, check.names = FALSE)
registry_row <- registry[registry$sample_key == sample_key, , drop = FALSE]
if (nrow(registry_row) != 1L) stop("Unknown sample_key: ", sample_key)

celltypes <- c(
  "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage",
  "Monocyte", "Neutrophil", "Dendritic", "B_cell", "T_NK",
  "Pericyte", "Schwann_cell"
)
immune_types <- c(
  "Macrophage", "Monocyte", "Neutrophil", "Dendritic", "B_cell", "T_NK"
)
marker_panel <- list(
  Cardiomyocyte = c("Tnnt2", "Myh6", "Myh7", "Actc1"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr", "Emcn"),
  Fibroblast = c("Col1a1", "Col1a2", "Dcn", "Pdgfra"),
  Macrophage = c("Adgre1", "C1qa", "C1qb", "Csf1r"),
  Monocyte = c("Ccr2", "Vcan", "Ly6c2", "Chil3"),
  Neutrophil = c("S100a8", "S100a9", "Csf3r", "Retnlg"),
  Dendritic = c("Flt3", "Clec9a", "Itgax", "Zbtb46"),
  B_cell = c("Cd79a", "Cd19", "Ms4a1", "Cd22"),
  T_NK = c("Cd3d", "Cd3e", "Trbc1", "Nkg7"),
  Pericyte = c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8"),
  Schwann_cell = c("S100b", "Plp1", "Sox10", "Mpz")
)
palette <- c(
  Cardiomyocyte = "#D55E00",
  Endothelial = "#0072B2",
  Fibroblast = "#009E73",
  Macrophage = "#CC79A7",
  Monocyte = "#C81D25",
  Neutrophil = "#E15759",
  Dendritic = "#7B2CBF",
  B_cell = "#4E79A7",
  T_NK = "#2F4B7C",
  Pericyte = "#E69F00",
  Schwann_cell = "#56B4E9",
  Uncertain = "#8C8C8C"
)

started <- stamp()
running_status <- data.frame(
  sample_key = sample_key,
  status = "RUNNING",
  start_time = started,
  finish_time = "",
  input_path = input_path,
  input_sha256 = "",
  output_path = output_object,
  output_sha256 = "",
  spots = registry_row$spots_after[[1L]],
  rctd_runtime_seconds = NA_real_,
  dwls_runtime_seconds = NA_real_,
  dominant_agreement_rate = NA_real_,
  median_cosine_similarity = NA_real_,
  median_celltype_spearman = NA_real_,
  uncertain_fraction = NA_real_,
  error = "",
  stringsAsFactors = FALSE
)
write_tsv(running_status, sample_status_path)

run_sample <- function() {
  input_hash_before <- sha256(input_path)
  features_table <- read.delim(feature_path, check.names = FALSE)
  features <- features_table$feature[as.logical(features_table$selected)]
  if (length(features) != 3000L) stop("Expected 3000 selected features")

  log_msg("Loading expanded reference")
  reference <- readRDS(reference_path)
  if (!inherits(reference, "Seurat") ||
      (reference_variant == "full_audit" && ncol(reference) != 176526L) ||
      (reference_variant != "full_audit" && ncol(reference) >= 50000L) ||
      !setequal(
        unique(as.character(reference$P09_expanded_celltype)), celltypes
      )) {
    stop("Expanded reference validation failed")
  }
  log_msg("Loading Phase 4 spatial object")
  object <- readRDS(input_path)
  if (!inherits(object, "Seurat") ||
      ncol(object) != registry_row$spots_after[[1L]] ||
      length(Images(object)) != 1L) {
    stop("Phase 4 object validation failed")
  }
  image_name <- Images(object)[[1L]]
  if (!all(features %in% rownames(object)) ||
      !all(features %in% rownames(reference))) {
    stop("Selected features are not shared")
  }
  spatial_counts <- LayerData(
    object[["Spatial"]], layer = "counts", fast = FALSE
  )
  selected_umi <- Matrix::colSums(
    spatial_counts[features, , drop = FALSE]
  )
  zero_information_spots <- names(selected_umi)[
    !is.finite(selected_umi) | selected_umi <= 0
  ]
  phase4_rctd_fallback <- build_phase4_fallback(
    object, "RCTD", colnames(object), celltypes
  )
  phase4_dwls_fallback <- build_phase4_fallback(
    object, "DWLS", colnames(object), celltypes
  )
  if (length(zero_information_spots)) {
    log_msg(
      "Expanded reference has zero selected-feature UMI for ",
      length(zero_information_spots),
      " spot(s); using explicitly flagged Phase 4 broad fallback"
    )
  }

  log_msg(
    "Starting expanded RCTD; spots=", ncol(object),
    " reference_cells=", ncol(reference),
    " reference_variant=", reference_variant,
    " classes=", length(celltypes)
  )
  write_tsv(
    data.frame(
      sample_key = sample_key,
      reference_variant = reference_variant,
      reference_path = reference_path,
      reference_sha256 = sha256(reference_path),
      reference_cells = ncol(reference),
      reference_features = nrow(reference),
      stringsAsFactors = FALSE
    ),
    file.path(result_dir, "reference_usage.tsv")
  )
  rctd_start <- proc.time()[["elapsed"]]
  object <- RunRCTD(
    srt = object,
    reference = reference,
    reference_label = "P09_expanded_celltype",
    assay = "Spatial",
    reference_assay = "RNA",
    layer = "counts",
    reference_layer = "counts",
    features = features,
    image = image_name,
    coord.cols = c("x", "y"),
    rctd_mode = "full",
    max_cores = 8,
    min_cells = 25,
    prefix = "P09_RCTD",
    store_results = TRUE,
    round_counts = TRUE,
    create_rctd_params = list(
      pixel_count_min = 1,
      UMI_min = 1,
      ref_UMI_min = 1,
      ref_n_cells_min = 25,
      ref_n_cells_max = 100000
    ),
    verbose = TRUE,
    coordinate_space = "raw",
    tool_name = "P09_RCTD"
  )
  rctd_seconds <- proc.time()[["elapsed"]] - rctd_start
  if (is.null(object@tools$P09_RCTD$weights)) {
    stop("Expanded RCTD weights were not stored")
  }
  rctd <- normalize_proportions(
    object@tools$P09_RCTD$proportions,
    colnames(object), celltypes, "P09_RCTD",
    fallback = phase4_rctd_fallback,
    force_fallback_spots = zero_information_spots
  )
  rctd_fallback_spots <- attr(rctd, "fallback_rows")
  object@tools$P09_RCTD$object <- NULL
  object@tools$P09_RCTD$proportions <- rctd
  build_result <- get("spatial_result_build", envir = asNamespace("scop"))
  object@tools$P09_RCTD <- build_result(
    bundle = object@tools$P09_RCTD,
    method = "RCTD",
    result_type = "deconvolution",
    provenance = list(
      producer = "RunRCTD",
      backend_id = "spacexr",
      phase = "Phase 9 expanded immune reference"
    )
  )
  log_msg("Expanded RCTD completed in ", round(rctd_seconds, 1), " seconds")
  invisible(gc())

  log_msg("Starting expanded SpatialDWLS")
  dwls_start <- proc.time()[["elapsed"]]
  object <- RunSpatialDWLS(
    srt = object,
    reference = reference,
    reference_label = "P09_expanded_celltype",
    assay = "Spatial",
    reference_assay = "RNA",
    layer = "counts",
    reference_layer = "counts",
    features = features,
    image = image_name,
    coord.cols = c("x", "y"),
    min_cells = 25,
    prefix = "P09_DWLS",
    tool_name = "P09_DWLS",
    normalize = TRUE,
    store_results = TRUE,
    verbose = TRUE,
    coordinate_space = "raw"
  )
  dwls_seconds <- proc.time()[["elapsed"]] - dwls_start
  dwls_result <- GetSpatialResult(object, tool_name = "P09_DWLS")
  dwls <- normalize_proportions(
    if (!is.null(dwls_result$proportions)) {
      dwls_result$proportions
    } else {
      dwls_result$weights
    },
    colnames(object), celltypes, "P09_DWLS",
    fallback = phase4_dwls_fallback,
    force_fallback_spots = zero_information_spots
  )
  dwls_fallback_spots <- attr(dwls, "fallback_rows")
  object@tools$P09_DWLS$proportions <- dwls
  log_msg("Expanded SpatialDWLS completed in ", round(dwls_seconds, 1), " seconds")

  fallback_spots <- union(rctd_fallback_spots, dwls_fallback_spots)
  expanded_inference_available <- !(colnames(object) %in% fallback_spots)
  names(expanded_inference_available) <- colnames(object)
  valid_index <- expanded_inference_available
  valid_spots <- names(expanded_inference_available)[
    expanded_inference_available
  ]
  if (length(valid_spots) < 2L) {
    stop("Fewer than two spots have valid expanded-reference inference")
  }
  rctd_dominant <- dominant_type(rctd)
  dwls_dominant <- dominant_type(dwls)
  methods_agree <- rctd_dominant == dwls_dominant
  methods_agree[!expanded_inference_available] <- FALSE
  rctd_top <- apply(rctd, 1L, max)
  rctd_margin <- top_margin(rctd)
  cosine <- cosine_rows(rctd, dwls)
  confidence <- ifelse(
    methods_agree & rctd_top >= 0.50 & rctd_margin >= 0.20 & cosine >= 0.85,
    "High",
    ifelse(
      methods_agree & rctd_top >= 0.35 & rctd_margin >= 0.10 &
        cosine >= 0.75,
      "Moderate", "Low"
    )
  )
  confidence[!expanded_inference_available] <- "Low"
  consensus_dominant <- ifelse(confidence == "Low", "Uncertain", rctd_dominant)
  confidence_score <- (
    0.35 * rctd_top +
      0.25 * pmin(1, rctd_margin / 0.5) +
      0.30 * cosine +
      0.10 * as.numeric(methods_agree)
  )
  confidence_score[!expanded_inference_available] <- 0
  inference_source <- ifelse(
    expanded_inference_available,
    "expanded_reference_RCTD_and_SpatialDWLS",
    "Phase4_broad_fallback_zero_selected_UMI"
  )

  meta_add <- data.frame(
    P09E_RCTD_dominant = rctd_dominant,
    P09E_DWLS_dominant = dwls_dominant,
    P09E_methods_agree = methods_agree,
    P09E_RCTD_max_prop = rctd_top,
    P09E_RCTD_margin = rctd_margin,
    P09E_method_cosine = cosine,
    P09E_confidence = factor(confidence, levels = c("Low", "Moderate", "High")),
    P09E_confidence_score = confidence_score,
    P09E_expanded_inference_available = expanded_inference_available,
    P09E_inference_source = inference_source,
    P09E_dominant_type = factor(
      consensus_dominant, levels = c(celltypes, "Uncertain")
    ),
    P09E_total_myeloid = rowSums(rctd[, c(
      "Macrophage", "Monocyte", "Neutrophil", "Dendritic"
    ), drop = FALSE]),
    P09E_total_lymphoid = rowSums(rctd[, c("B_cell", "T_NK"), drop = FALSE]),
    row.names = colnames(object),
    stringsAsFactors = FALSE
  )
  for (celltype in celltypes) {
    meta_add[[paste0("P09E_RCTD_prop_", celltype)]] <- rctd[, celltype]
    meta_add[[paste0("P09E_DWLS_prop_", celltype)]] <- dwls[, celltype]
  }

  total_umi <- Matrix::colSums(spatial_counts)
  marker_rows <- list()
  for (celltype in celltypes) {
    genes <- intersect(marker_panel[[celltype]], rownames(spatial_counts))
    marker_sum <- Matrix::colSums(spatial_counts[genes, , drop = FALSE])
    marker_score <- log1p(1e4 * marker_sum / pmax(total_umi, 1))
    meta_add[[paste0("P09E_marker_score_", celltype)]] <-
      marker_score[rownames(meta_add)]
    marker_rows[[celltype]] <- data.frame(
      sample_key = sample_key,
      celltype = celltype,
      markers_available = paste(genes, collapse = ";"),
      marker_positive_fraction = mean(marker_sum > 0),
      rctd_marker_spearman = suppressWarnings(cor(
        rctd[valid_spots, celltype], marker_score[valid_spots],
        method = "spearman", use = "pairwise.complete.obs"
      )),
      dwls_marker_spearman = suppressWarnings(cor(
        dwls[valid_spots, celltype], marker_score[valid_spots],
        method = "spearman", use = "pairwise.complete.obs"
      )),
      stringsAsFactors = FALSE
    )
  }
  object <- AddMetaData(object, meta_add)
  fallback_audit <- data.frame(
    sample_key = rep(sample_key, length(fallback_spots)),
    barcode = fallback_spots,
    reason = rep(
      "zero selected-feature UMI; no expanded-reference inference",
      length(fallback_spots)
    ),
    selected_feature_UMI = as.numeric(selected_umi[fallback_spots]),
    total_UMI = as.numeric(total_umi[fallback_spots]),
    rctd_fallback_used = fallback_spots %in% rctd_fallback_spots,
    dwls_fallback_used = fallback_spots %in% dwls_fallback_spots,
    phase4_rctd_dominant = as.character(
      object$P4_RCTD_dominant[match(fallback_spots, colnames(object))]
    ),
    phase4_dwls_dominant = as.character(
      object$P4_DWLS_dominant[match(fallback_spots, colnames(object))]
    ),
    phase4_consensus = as.character(
      object$P4_dominant_type[match(fallback_spots, colnames(object))]
    ),
    phase4_confidence = as.character(
      object$P4_confidence[match(fallback_spots, colnames(object))]
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    fallback_audit,
    file.path(result_dir, "expanded_inference_fallback.tsv")
  )

  method_cor <- vapply(celltypes, function(celltype) {
    suppressWarnings(cor(
      rctd[valid_spots, celltype], dwls[valid_spots, celltype],
      method = "spearman", use = "pairwise.complete.obs"
    ))
  }, numeric(1L))
  write_tsv(
    data.frame(
      sample_key = sample_key,
      celltype = celltypes,
      spearman_rctd_vs_dwls = method_cor,
      mean_RCTD = colMeans(rctd[valid_spots, , drop = FALSE]),
      mean_SpatialDWLS = colMeans(dwls[valid_spots, , drop = FALSE]),
      stringsAsFactors = FALSE
    ),
    file.path(result_dir, "method_agreement_by_celltype.tsv")
  )
  write_tsv(
    do.call(rbind, marker_rows),
    file.path(result_dir, "marker_concordance.tsv")
  )

  spot_table <- data.frame(
    sample_key = sample_key,
    GSM = registry_row$GSM[[1L]],
    sample_id = registry_row$sample_id[[1L]],
    condition = registry_row$condition[[1L]],
    timepoint = registry_row$timepoint[[1L]],
    replicate = registry_row$replicate[[1L]],
    analysis_role = registry_row$analysis_role[[1L]],
    include_primary_MI_analysis =
      registry_row$include_primary_MI_analysis[[1L]],
    barcode = colnames(object),
    rctd_dominant = rctd_dominant,
    dwls_dominant = dwls_dominant,
    methods_agree = methods_agree,
    rctd_max_prop = rctd_top,
    rctd_margin = rctd_margin,
    method_cosine = cosine,
    confidence = confidence,
    dominant_type = consensus_dominant,
    total_myeloid = meta_add$P09E_total_myeloid,
    total_lymphoid = meta_add$P09E_total_lymphoid,
    old_P4_macrophage = as.numeric(object$P4_RCTD_prop_Macrophage),
    stringsAsFactors = FALSE
  )
  for (celltype in celltypes) {
    spot_table[[paste0("RCTD_", celltype)]] <- rctd[, celltype]
    spot_table[[paste0("SpatialDWLS_", celltype)]] <- dwls[, celltype]
  }
  write_tsv(spot_table, file.path(result_dir, "spot_proportions.tsv"))

  p_immune <- SpatialDeconvolutionPlot(
    object,
    tool_name = "P09_RCTD",
    cell_types = immune_types,
    plot_type = "point",
    ncol = 3,
    image = image_name,
    crop = FALSE,
    coord.cols = c("x", "y"),
    overlay_image = TRUE,
    image.alpha = 0.8,
    pt.alpha = 0.9,
    legend.position = "right"
  ) + plot_annotation(
    title = paste0(sample_key, ": expanded immune RCTD")
  )
  ggsave(
    file.path(result_dir, "expanded_immune_rctd.png"),
    p_immune, width = 13, height = 8, dpi = 240, bg = "white"
  )

  p_dominant <- SpatialSpotPlot(
    object,
    values = setNames(consensus_dominant, colnames(object)),
    plot_type = "point",
    image = image_name,
    crop = FALSE,
    coord.cols = c("x", "y"),
    overlay_image = TRUE,
    image.alpha = 0.85,
    pt.alpha = 0.85,
    palcolor = palette,
    legend.title = "Expanded component"
  ) + ggtitle(
    paste0(sample_key, ": expanded-reference consensus"),
    subtitle = "Low-confidence spots remain Uncertain"
  )
  ggsave(
    file.path(result_dir, "expanded_consensus_dominant.png"),
    p_dominant, width = 9, height = 8, dpi = 240, bg = "white"
  )

  object@misc$phase_09_expanded_deconvolution <- list(
    main_method = "RCTD full mode",
    sensitivity_method = "SpatialDWLS",
    reference_variant = reference_variant,
    reference_path = reference_path,
    reference_cells = ncol(reference),
    reference_features = length(features),
    celltypes = celltypes,
    expanded_inference_fallback_count = length(fallback_spots),
    expanded_inference_fallback_barcodes = fallback_spots,
    expanded_inference_fallback_reason =
      "zero selected-feature UMI; Phase 4 broad fallback; consensus Uncertain",
    spots_are_mixtures = TRUE,
    prior_phase4_preserved = TRUE
  )
  saveRDS(object, output_object, compress = FALSE)
  readback <- readRDS(output_object)
  if (!inherits(readback, "Seurat") ||
      ncol(readback) != ncol(object) ||
      !all(c(
        "P09E_dominant_type", "P09E_confidence", "P09E_total_myeloid",
        "P09E_expanded_inference_available", "P09E_inference_source"
      ) %in% colnames(readback[[]]))) {
    stop("Expanded sample readback failed")
  }
  input_hash_after <- sha256(input_path)
  if (!identical(input_hash_before, input_hash_after)) {
    stop("Phase 4 input object changed")
  }
  output_hash <- sha256(output_object)
  final_status <- running_status
  final_status$status <- "PASS"
  final_status$finish_time <- stamp()
  final_status$input_sha256 <- input_hash_before
  final_status$output_sha256 <- output_hash
  final_status$rctd_runtime_seconds <- rctd_seconds
  final_status$dwls_runtime_seconds <- dwls_seconds
  final_status$dominant_agreement_rate <- mean(methods_agree[valid_index])
  final_status$median_cosine_similarity <- median(cosine[valid_index])
  final_status$median_celltype_spearman <- median(method_cor, na.rm = TRUE)
  final_status$uncertain_fraction <- mean(consensus_dominant == "Uncertain")
  write_tsv(final_status, sample_status_path)
  log_msg(
    "PASS valid_spots=", length(valid_spots),
    " fallback_spots=", length(fallback_spots),
    " agreement=", round(mean(methods_agree[valid_index]), 3),
    " cosine=", round(median(cosine[valid_index]), 3),
    " median_spearman=", round(median(method_cor, na.rm = TRUE), 3)
  )
}

tryCatch(
  run_sample(),
  error = function(error) {
    failed <- running_status
    failed$status <- "FAIL"
    failed$finish_time <- stamp()
    failed$error <- conditionMessage(error)
    write_tsv(failed, sample_status_path)
    log_msg("FAIL: ", conditionMessage(error))
    stop(error)
  }
)
