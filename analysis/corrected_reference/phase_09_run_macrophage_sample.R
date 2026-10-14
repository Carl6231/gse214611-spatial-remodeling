#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: phase_09_run_macrophage_sample.R <root> <sample_key>")
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
sample_key <- args[[2L]]
.libPaths(c(file.path(root, "Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(scop)
})
set.seed(20260731L)

for (path in c(
  file.path(root, "logs", "phase_09_macrophage_samples"),
  file.path(root, "metadata", "phase_09_macrophage_sample_status"),
  file.path(root, "objects", "phase_09", "macrophage_mapping"),
  file.path(root, "results", "phase_09", "macrophage_mapping", sample_key)
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
log_path <- file.path(
  root, "logs", "phase_09_macrophage_samples", paste0(sample_key, ".log")
)
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
normalize_props <- function(x, spots, types, label) {
  x <- as.matrix(x)
  if (all(spots %in% colnames(x)) && !all(spots %in% rownames(x))) x <- t(x)
  if (!all(spots %in% rownames(x)) || !all(types %in% colnames(x))) {
    stop(label, " proportions do not cover spots/types")
  }
  x <- x[spots, types, drop = FALSE]
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  if (any(totals <= 0)) stop(label, " contains zero-sum rows")
  x / totals
}
cosine_rows <- function(x, y) {
  value <- rowSums(x * y) / sqrt(rowSums(x^2) * rowSums(y^2))
  value[!is.finite(value)] <- 0
  pmax(0, pmin(1, value))
}
top_margin <- function(x) apply(x, 1L, function(z) {
  ordered <- sort(z, decreasing = TRUE)
  ordered[[1L]] - ordered[[2L]]
})

status_path <- file.path(
  root, "metadata", "phase_09_macrophage_sample_status",
  paste0(sample_key, ".tsv")
)
input_path <- file.path(
  root, "objects", "phase_09", "deconvolution",
  paste0(sample_key, "_phase09_expanded.rds")
)
output_path <- file.path(
  root, "objects", "phase_09", "macrophage_mapping",
  paste0(sample_key, "_phase09_macrophage.rds")
)
reference_path <- file.path(
  root, "objects", "phase_09", "reference_macrophage_refined.rds"
)
for (path in c(input_path, reference_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

reference <- readRDS(reference_path)
if (!inherits(reference, "Seurat")) stop("Refined macrophage reference unreadable")
features <- rownames(reference)
types <- levels(reference$P09_macrophage_state)
if (length(types) != 5L || any(table(reference$P09_macrophage_state) < 25L)) {
  stop("Refined macrophage state reference failed validation")
}

started <- stamp()
sample_status <- data.frame(
  sample_key = sample_key,
  status = "RUNNING",
  start_time = started,
  finish_time = "",
  input_path = input_path,
  input_sha256 = "",
  output_path = output_path,
  output_sha256 = "",
  total_spots = NA_integer_,
  candidate_spots = NA_integer_,
  mapped_spots = NA_integer_,
  median_cosine = NA_real_,
  dominant_agreement = NA_real_,
  uncertain_fraction = NA_real_,
  error = "",
  stringsAsFactors = FALSE
)
write_tsv(sample_status, status_path)

run_sample <- function() {
  input_hash_before <- sha256(input_path)
  object <- readRDS(input_path)
  if (!inherits(object, "Seurat") ||
      !"P09E_RCTD_prop_Macrophage" %in% colnames(object[[]])) {
    stop("Expanded deconvolution input is incompatible")
  }
  candidate_spots <- colnames(object)[
    as.numeric(object$P09E_RCTD_prop_Macrophage) >= 0.10
  ]
  if (length(candidate_spots) < 20L) {
    stop("Fewer than 20 corrected-macrophage candidate spots")
  }
  candidate <- subset(object, cells = candidate_spots)
  image_name <- Images(candidate)[[1L]]
  if (!all(features %in% rownames(candidate))) stop("Mapping features are missing")

  log_msg(
    "Starting refined macrophage RCTD; candidate_spots=",
    length(candidate_spots), " states=", length(types)
  )
  candidate <- RunRCTD(
    srt = candidate,
    reference = reference,
    reference_label = "P09_macrophage_state",
    assay = "Spatial",
    reference_assay = "RNA",
    layer = "counts",
    reference_layer = "counts",
    features = features,
    image = image_name,
    coord.cols = c("x", "y"),
    rctd_mode = "full",
    max_cores = 8L,
    min_cells = 25L,
    prefix = "P09M_RCTD",
    tool_name = "P09M_RCTD",
    create_rctd_params = list(
      pixel_count_min = 1,
      UMI_min = 1,
      ref_UMI_min = 1,
      UMI_min_sigma = 1,
      ref_n_cells_min = 25,
      ref_n_cells_max = 100000
    ),
    verbose = TRUE,
    coordinate_space = "raw"
  )
  if (is.null(candidate@tools$P09M_RCTD$weights)) {
    stop("Refined macrophage RCTD weights were not stored")
  }
  rctd <- normalize_props(
    candidate@tools$P09M_RCTD$proportions,
    candidate_spots, types, "P09M_RCTD"
  )
  candidate@tools$P09M_RCTD$object <- NULL
  candidate@tools$P09M_RCTD$proportions <- rctd

  log_msg("Starting refined macrophage SpatialDWLS")
  candidate <- RunSpatialDWLS(
    srt = candidate,
    reference = reference,
    reference_label = "P09_macrophage_state",
    assay = "Spatial",
    reference_assay = "RNA",
    layer = "counts",
    reference_layer = "counts",
    features = features,
    image = image_name,
    coord.cols = c("x", "y"),
    min_cells = 25L,
    prefix = "P09M_DWLS",
    tool_name = "P09M_DWLS",
    normalize = TRUE,
    store_results = TRUE,
    verbose = TRUE,
    coordinate_space = "raw"
  )
  dwls_result <- GetSpatialResult(candidate, tool_name = "P09M_DWLS")
  dwls <- normalize_props(
    if (!is.null(dwls_result$proportions)) {
      dwls_result$proportions
    } else {
      dwls_result$weights
    },
    candidate_spots, types, "P09M_DWLS"
  )
  candidate@tools$P09M_DWLS$proportions <- dwls

  rctd_dominant <- types[max.col(rctd, ties.method = "first")]
  dwls_dominant <- types[max.col(dwls, ties.method = "first")]
  methods_agree <- rctd_dominant == dwls_dominant
  top <- apply(rctd, 1L, max)
  margin <- top_margin(rctd)
  cosine <- cosine_rows(rctd, dwls)
  confidence <- ifelse(
    methods_agree & top >= 0.50 & margin >= 0.20 & cosine >= 0.85,
    "High",
    ifelse(
      methods_agree & top >= 0.35 & margin >= 0.10 & cosine >= 0.75,
      "Moderate", "Low"
    )
  )
  final_state <- ifelse(confidence == "Low", "Uncertain", rctd_dominant)

  metadata <- data.frame(
    P09M_candidate = colnames(object) %in% candidate_spots,
    P09M_RCTD_dominant = NA_character_,
    P09M_DWLS_dominant = NA_character_,
    P09M_methods_agree = NA,
    P09M_RCTD_max_prop = NA_real_,
    P09M_RCTD_margin = NA_real_,
    P09M_method_cosine = NA_real_,
    P09M_confidence = "Not_mapped",
    P09M_dominant_state = "Not_macrophage_candidate",
    row.names = colnames(object),
    stringsAsFactors = FALSE
  )
  metadata[candidate_spots, "P09M_RCTD_dominant"] <- rctd_dominant
  metadata[candidate_spots, "P09M_DWLS_dominant"] <- dwls_dominant
  metadata[candidate_spots, "P09M_methods_agree"] <- methods_agree
  metadata[candidate_spots, "P09M_RCTD_max_prop"] <- top
  metadata[candidate_spots, "P09M_RCTD_margin"] <- margin
  metadata[candidate_spots, "P09M_method_cosine"] <- cosine
  metadata[candidate_spots, "P09M_confidence"] <- confidence
  metadata[candidate_spots, "P09M_dominant_state"] <- final_state
  for (type in types) {
    metadata[[paste0("P09M_RCTD_prop_", type)]] <- NA_real_
    metadata[[paste0("P09M_DWLS_prop_", type)]] <- NA_real_
    metadata[candidate_spots, paste0("P09M_RCTD_prop_", type)] <- rctd[, type]
    metadata[candidate_spots, paste0("P09M_DWLS_prop_", type)] <- dwls[, type]
  }
  object <- AddMetaData(object, metadata)
  object@tools$P09M_RCTD <- candidate@tools$P09M_RCTD
  object@tools$P09M_DWLS <- candidate@tools$P09M_DWLS
  object@misc$phase_09_macrophage_mapping <- list(
    candidate_threshold = "P09E_RCTD_prop_Macrophage >= 0.10",
    candidate_spots = length(candidate_spots),
    reference_path = reference_path,
    states = types,
    excluded_phase6_clusters = c("0", "4", "5", "6", "7", "9", "10"),
    spots_are_mixtures = TRUE
  )
  saveRDS(object, output_path, compress = FALSE)
  readback <- readRDS(output_path)
  if (!inherits(readback, "Seurat") ||
      ncol(readback) != ncol(object) ||
      sum(readback$P09M_candidate) != length(candidate_spots)) {
    stop("Refined macrophage output readback failed")
  }

  spot_table <- cbind(
    data.frame(
      sample_key = sample_key,
      barcode = candidate_spots,
      corrected_macrophage_prop =
        as.numeric(object$P09E_RCTD_prop_Macrophage[candidate_spots]),
      old_P4_macrophage_prop =
        as.numeric(object$P4_RCTD_prop_Macrophage[candidate_spots]),
      P09M_RCTD_dominant = rctd_dominant,
      P09M_DWLS_dominant = dwls_dominant,
      P09M_methods_agree = methods_agree,
      P09M_RCTD_max_prop = top,
      P09M_RCTD_margin = margin,
      P09M_method_cosine = cosine,
      P09M_confidence = confidence,
      P09M_dominant_state = final_state,
      stringsAsFactors = FALSE
    ),
    setNames(as.data.frame(rctd), paste0("RCTD_", types)),
    setNames(as.data.frame(dwls), paste0("DWLS_", types))
  )
  write_tsv(
    spot_table,
    file.path(
      root, "results", "phase_09", "macrophage_mapping",
      sample_key, "spot_mapping.tsv"
    )
  )

  input_hash_after <- sha256(input_path)
  if (!identical(input_hash_before, input_hash_after)) {
    stop("Expanded deconvolution input changed")
  }
  sample_status$status <- "PASS"
  sample_status$finish_time <- stamp()
  sample_status$input_sha256 <- input_hash_before
  sample_status$output_sha256 <- sha256(output_path)
  sample_status$total_spots <- ncol(object)
  sample_status$candidate_spots <- length(candidate_spots)
  sample_status$mapped_spots <- nrow(rctd)
  sample_status$median_cosine <- median(cosine)
  sample_status$dominant_agreement <- mean(methods_agree)
  sample_status$uncertain_fraction <- mean(final_state == "Uncertain")
  write_tsv(sample_status, status_path)
  log_msg(
    "PASS candidate_spots=", length(candidate_spots),
    " agreement=", round(mean(methods_agree), 3),
    " cosine=", round(median(cosine), 3)
  )
}

tryCatch(
  run_sample(),
  error = function(error) {
    sample_status$status <- "FAIL"
    sample_status$finish_time <- stamp()
    sample_status$error <- conditionMessage(error)
    write_tsv(sample_status, status_path)
    log_msg("FAIL: ", conditionMessage(error))
    stop(error)
  }
)
