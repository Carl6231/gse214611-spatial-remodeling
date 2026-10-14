#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else {
  "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
}
.libPaths(c(file.path(root, "Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(scop)
  library(cluster)
})
if (!"SpatialSpotPlot" %in% getNamespaceExports("scop")) {
  stop("scop::SpatialSpotPlot is unavailable")
}
set.seed(20260731L)

for (path in c(
  file.path(root, "objects", "phase_09", "niches"),
  file.path(root, "metadata"),
  file.path(root, "results", "phase_09", "niche_robustness"),
  file.path(root, "results", "phase_09", "niche_robustness", "visual_qc"),
  file.path(root, "logs")
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
log_path <- file.path(root, "logs", "phase_09_niche_robustness.log")
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
robust_z <- function(x) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) spread <- sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) return(rep(0, length(x)))
  value <- (x - center) / spread
  value[!is.finite(value)] <- 0
  value
}
row_mean <- function(...) {
  matrix <- cbind(...)
  rowMeans(matrix, na.rm = TRUE)
}
adjusted_rand <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- as.character(x[keep])
  y <- as.character(y[keep])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  index <- sum(choose2(tab))
  sum_a <- sum(choose2(rowSums(tab)))
  sum_b <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- sum_a * sum_b / total
  maximum <- (sum_a + sum_b) / 2
  denominator <- maximum - expected
  if (!is.finite(denominator) || abs(denominator) < 1e-12) {
    return(as.numeric(identical(x, y)))
  }
  (index - expected) / denominator
}
jaccard <- function(x, y, label) {
  x_set <- which(as.character(x) == label)
  y_set <- which(as.character(y) == label)
  union <- length(union(x_set, y_set))
  if (union == 0L) return(NA_real_)
  length(intersect(x_set, y_set)) / union
}
exact_effect <- function(x, y) mean(x) - mean(y)

samples <- c(
  "Sham_Day3_rep1", "Sham_Day3_rep2", "MI_1HR_rep1", "MI_4HR_rep1",
  "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
  "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3",
  "NeedlePass_Day3_rep3", "NeedlePass_Day3_rep4"
)
primary_samples <- c(
  "Sham_Day3_rep1", "Sham_Day3_rep2",
  "MI_1HR_rep1", "MI_4HR_rep1",
  "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
  "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3"
)
broad_types <- c(
  "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage",
  "Monocyte", "Neutrophil", "Dendritic", "B_cell", "T_NK",
  "Pericyte", "Schwann_cell"
)
states <- c(
  "Resident_Folr2_macrophage",
  "Trem2_Spp1_macrophage",
  "MHCII_macrophage",
  "Cycling_macrophage",
  "Stress_response_macrophage"
)
niche_keys <- c(
  "inflammatory_infarct",
  "macrophage_fibroblast_remodeling",
  "vascular_repair",
  "remote_resident"
)
niche_labels <- c(
  inflammatory_infarct = "inflammatory-infarct niche",
  macrophage_fibroblast_remodeling =
    "macrophage-fibroblast remodeling niche",
  vascular_repair = "vascular-repair niche",
  remote_resident = "remote resident niche"
)
niche_palette <- c(
  "inflammatory-infarct niche" = "#C81D25",
  "macrophage-fibroblast remodeling niche" = "#D9822B",
  "vascular-repair niche" = "#2A8C82",
  "remote resident niche" = "#2F6FB0",
  "Uncertain" = "#8C8C8C",
  "Not_macrophage_center" = "#D0D0D0"
)
module_sets <- list(
  inflammation = c(
    "Il1b", "Tnf", "Nfkbia", "Ccl2", "Ccl3", "Ccl4", "Cxcl2",
    "S100a8", "S100a9", "Lcn2", "Ptgs2", "Nlrp3"
  ),
  ECM = c(
    "Col1a1", "Col1a2", "Col3a1", "Fn1", "Postn", "Tnc",
    "Sparc", "Thbs1", "Lox", "Ctgf"
  ),
  vascular = c(
    "Pecam1", "Cdh5", "Kdr", "Emcn", "Klf2", "Klf4",
    "Egfl7", "Eng", "Ramp2", "Rbp7", "Vwf", "Esam"
  ),
  cardiomyocyte_stress = c(
    "Nppa", "Nppb", "Myh7", "Ankrd1", "Xirp2", "Acta1",
    "Atf3", "Fos", "Jun", "Hspa1a", "Hspa1b", "Ddit3"
  )
)
thresholds <- c(0.15, 0.20, 0.25, 0.30)
radius_multipliers <- c(1.0, 1.2, 1.5)
program_gates <- c(0.05, 0.10, 0.15)
margin_gates <- c(0.10, 0.15, 0.20)
grid <- expand.grid(
  center_threshold = thresholds,
  radius_multiplier = radius_multipliers,
  program_gate = program_gates,
  margin_gate = margin_gates,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$setting_id <- sprintf("S%03d", seq_len(nrow(grid)))
baseline_id <- grid$setting_id[
  grid$center_threshold == 0.20 &
    grid$radius_multiplier == 1.2 &
    grid$program_gate == 0.10 &
    grid$margin_gate == 0.15
]
if (nrow(grid) != 108L || length(baseline_id) != 1L) {
  stop("Robustness grid definition failed")
}
write_tsv(grid, file.path(root, "metadata", "phase_09_niche_parameter_grid.tsv"))

registry <- read.delim(
  file.path(root, "metadata", "phase_03_sample_registry.tsv"),
  check.names = FALSE
)
feature_rows <- list()
edge_rows <- list()
module_audit_rows <- list()

for (sample in samples) {
  log_msg("Building corrected neighborhood features for ", sample)
  object_path <- file.path(
    root, "objects", "phase_09", "macrophage_mapping",
    paste0(sample, "_phase09_macrophage.rds")
  )
  object <- readRDS(object_path)
  if (!inherits(object, "Seurat") || length(Images(object)) != 1L) {
    stop("Invalid Phase 9 macrophage object: ", sample)
  }
  spots <- colnames(object)
  image_name <- Images(object)[[1L]]
  coordinates <- GetTissueCoordinates(object[[image_name]])
  if ("cell" %in% colnames(coordinates)) {
    rownames(coordinates) <- as.character(coordinates$cell)
  }
  if (!all(spots %in% rownames(coordinates))) {
    stop("Coordinates do not cover all spots: ", sample)
  }
  coordinates <- as.matrix(
    coordinates[spots, c("x", "y"), drop = FALSE]
  )
  storage.mode(coordinates) <- "double"
  distance <- as.matrix(dist(coordinates))
  diag(distance) <- Inf
  lattice_spacing <- median(apply(distance, 1L, min))
  diag(distance) <- 0

  meta <- object[[]]
  broad_columns <- paste0("P09E_RCTD_prop_", broad_types)
  state_columns <- paste0("P09M_RCTD_prop_", states)
  required_columns <- c(
    broad_columns, state_columns, "P09E_RCTD_prop_Macrophage",
    "P09M_candidate"
  )
  if (!all(required_columns %in% colnames(meta))) {
    stop("Required corrected metadata missing: ", sample)
  }
  total_myeloid <- rowSums(meta[, paste0(
    "P09E_RCTD_prop_",
    c("Macrophage", "Monocyte", "Neutrophil", "Dendritic")
  ), drop = FALSE])
  candidate_centers <- spots[
    total_myeloid >= min(thresholds) &
      as.logical(meta$P09M_candidate)
  ]
  if (length(candidate_centers) < 20L) {
    stop("Fewer than 20 candidate centers: ", sample)
  }
  center_indices <- match(candidate_centers, spots)

  counts <- LayerData(object[["Spatial"]], layer = "counts", fast = FALSE)
  module_genes <- unique(unlist(module_sets, use.names = FALSE))
  available_module_genes <- intersect(module_genes, rownames(counts))
  total_umi <- Matrix::colSums(counts)
  log_expression <- log1p(
    1e4 * t(t(counts[available_module_genes, , drop = FALSE]) /
                pmax(total_umi, 1))
  )
  gene_z <- t(scale(t(as.matrix(log_expression))))
  gene_z[!is.finite(gene_z)] <- 0

  for (radius_multiplier in radius_multipliers) {
    radius <- lattice_spacing * radius_multiplier
    neighbor_indices <- lapply(center_indices, function(index) {
      candidates <- which(distance[index, ] <= radius)
      candidates[order(distance[index, candidates], method = "radix")]
    })
    edge_table <- do.call(rbind, lapply(seq_along(candidate_centers), function(i) {
      index <- neighbor_indices[[i]]
      data.frame(
        sample_key = sample,
        neighbor_sample_key = sample,
        radius_multiplier = radius_multiplier,
        center_barcode = candidate_centers[[i]],
        neighbor_barcode = spots[index],
        distance_pixels = distance[center_indices[[i]], index],
        is_center = spots[index] == candidate_centers[[i]],
        stringsAsFactors = FALSE
      )
    }))
    edge_rows[[paste(sample, radius_multiplier)]] <- edge_table

    feature <- data.frame(
      sample_key = sample,
      GSM = registry$GSM[match(sample, registry$sample_key)],
      condition = registry$condition[match(sample, registry$sample_key)],
      timepoint = registry$timepoint[match(sample, registry$sample_key)],
      replicate = registry$replicate[match(sample, registry$sample_key)],
      analysis_role =
        registry$analysis_role[match(sample, registry$sample_key)],
      include_primary_MI_analysis =
        registry$include_primary_MI_analysis[
          match(sample, registry$sample_key)
        ],
      center_barcode = candidate_centers,
      center_macrophage_prop =
        as.numeric(meta[candidate_centers, "P09E_RCTD_prop_Macrophage"]),
      center_total_myeloid_prop =
        as.numeric(total_myeloid[candidate_centers]),
      radius_multiplier = radius_multiplier,
      lattice_spacing_pixels = lattice_spacing,
      radius_pixels = radius,
      neighborhood_size = lengths(neighbor_indices),
      stringsAsFactors = FALSE
    )

    for (type in broad_types) {
      values <- as.numeric(meta[[paste0("P09E_RCTD_prop_", type)]])
      feature[[paste0("local_cell_", type)]] <- vapply(
        neighbor_indices, function(index) mean(values[index]), numeric(1L)
      )
    }
    macrophage_weight <- as.numeric(meta$P09E_RCTD_prop_Macrophage)
    mapped <- as.logical(meta$P09M_candidate)
    for (state in states) {
      values <- as.numeric(meta[[paste0("P09M_RCTD_prop_", state)]])
      feature[[paste0("local_state_", state)]] <- vapply(
        neighbor_indices,
        function(index) {
          keep <- index[mapped[index] & is.finite(values[index])]
          if (!length(keep) || sum(macrophage_weight[keep]) <= 0) return(0)
          weighted.mean(values[keep], macrophage_weight[keep])
        },
        numeric(1L)
      )
    }
    feature$local_mapped_macrophage_neighbors <- vapply(
      neighbor_indices,
      function(index) sum(mapped[index]),
      integer(1L)
    )

    for (module in names(module_sets)) {
      genes <- intersect(module_sets[[module]], rownames(gene_z))
      if (length(genes) < 4L) {
        stop("Fewer than four genes for module ", module, " in ", sample)
      }
      local_gene_means <- vapply(
        genes,
        function(gene) vapply(
          neighbor_indices,
          function(index) mean(gene_z[gene, index]),
          numeric(1L)
        ),
        numeric(length(neighbor_indices))
      )
      if (is.null(dim(local_gene_means))) {
        local_gene_means <- matrix(
          local_gene_means,
          nrow = length(neighbor_indices),
          dimnames = list(NULL, genes)
        )
      }
      feature[[paste0("module_", module)]] <- rowMeans(local_gene_means)
      feature[[paste0("module_", module, "_positive_genes")]] <-
        rowSums(local_gene_means > 0)
      module_audit_rows[[length(module_audit_rows) + 1L]] <- data.frame(
        sample_key = sample,
        radius_multiplier = radius_multiplier,
        module = module,
        genes_defined = length(module_sets[[module]]),
        genes_detected = length(genes),
        genes = paste(genes, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
    feature_rows[[paste(sample, radius_multiplier)]] <- feature
  }
  rm(object, counts, log_expression, gene_z, distance)
  gc()
}

features <- do.call(rbind, feature_rows)
rownames(features) <- NULL
edges <- do.call(rbind, edge_rows)
rownames(edges) <- NULL
module_audit <- do.call(rbind, module_audit_rows)
write_tsv_gz(
  features,
  file.path(root, "metadata", "phase_09_corrected_neighborhood_features.tsv.gz")
)
write_tsv_gz(
  edges,
  file.path(root, "metadata", "phase_09_corrected_neighborhood_edges.tsv.gz")
)
write_tsv(
  module_audit,
  file.path(root, "metadata", "phase_09_corrected_module_audit.tsv")
)

scale_columns <- c(
  paste0("local_cell_", c(
    "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage",
    "Monocyte", "Neutrophil"
  )),
  paste0("local_state_", states),
  paste0("module_", names(module_sets))
)
for (radius_multiplier in radius_multipliers) {
  index <- features$radius_multiplier == radius_multiplier
  for (column in scale_columns) {
    features[index, paste0("z_", column)] <- robust_z(features[index, column])
  }
}

features$score_inflammatory_infarct <- row_mean(
  features$z_local_cell_Macrophage,
  row_mean(
    features$z_local_cell_Monocyte,
    features$z_local_cell_Neutrophil
  ),
  features$z_module_inflammation,
  features$z_module_cardiomyocyte_stress
)
features$score_macrophage_fibroblast_remodeling <- row_mean(
  features$z_local_cell_Macrophage,
  features$z_local_cell_Fibroblast,
  features$z_local_state_Trem2_Spp1_macrophage,
  features$z_module_ECM
)
features$score_vascular_repair <- row_mean(
  features$z_local_cell_Macrophage,
  features$z_local_cell_Endothelial,
  row_mean(
    features$z_local_state_Resident_Folr2_macrophage,
    features$z_local_state_MHCII_macrophage
  ),
  features$z_module_vascular
)
features$score_remote_resident <- row_mean(
  features$z_local_state_Resident_Folr2_macrophage,
  features$z_local_cell_Cardiomyocyte,
  -features$z_module_inflammation,
  -features$z_module_ECM,
  -features$z_module_cardiomyocyte_stress
)

assignment_rows <- list()
summary_rows <- list()
for (setting_index in seq_len(nrow(grid))) {
  setting <- grid[setting_index, ]
  current <- features[
    features$radius_multiplier == setting$radius_multiplier &
      features$center_total_myeloid_prop >= setting$center_threshold, ,
    drop = FALSE
  ]
  radius_feature <- features[
    features$radius_multiplier == setting$radius_multiplier, ,
    drop = FALSE
  ]
  medians <- vapply(
    c("local_cell_Cardiomyocyte", "local_cell_Endothelial",
      "local_cell_Fibroblast"),
    function(column) median(radius_feature[[column]], na.rm = TRUE),
    numeric(1L)
  )
  gate <- cbind(
    inflammatory_infarct = (
      (current$local_cell_Monocyte + current$local_cell_Neutrophil) >=
        setting$program_gate &
        current$module_inflammation > 0 &
        current$module_inflammation_positive_genes >= 2L &
        current$module_cardiomyocyte_stress_positive_genes >= 2L
    ),
    macrophage_fibroblast_remodeling = (
      current$local_cell_Fibroblast >=
        medians[["local_cell_Fibroblast"]] &
        current$local_state_Trem2_Spp1_macrophage >=
          setting$program_gate &
        current$module_ECM > 0 &
        current$module_ECM_positive_genes >= 2L
    ),
    vascular_repair = (
      current$local_cell_Endothelial >=
        medians[["local_cell_Endothelial"]] &
        (current$local_state_Resident_Folr2_macrophage +
           current$local_state_MHCII_macrophage) >=
          setting$program_gate &
        current$module_vascular > 0 &
        current$module_vascular_positive_genes >= 2L
    ),
    remote_resident = (
      current$local_cell_Cardiomyocyte >=
        medians[["local_cell_Cardiomyocyte"]] &
        current$local_state_Resident_Folr2_macrophage >=
          setting$program_gate &
        current$module_inflammation < 0 &
        current$module_ECM < 0 &
        current$module_cardiomyocyte_stress < 0
    )
  )
  gate[is.na(gate)] <- FALSE
  score_matrix <- cbind(
    inflammatory_infarct = current$score_inflammatory_infarct,
    macrophage_fibroblast_remodeling =
      current$score_macrophage_fibroblast_remodeling,
    vascular_repair = current$score_vascular_repair,
    remote_resident = current$score_remote_resident
  )
  score_matrix[!gate] <- -Inf
  best_index <- max.col(score_matrix, ties.method = "first")
  best_score <- score_matrix[cbind(seq_len(nrow(score_matrix)), best_index)]
  ordered <- t(apply(score_matrix, 1L, sort, decreasing = TRUE))
  margin <- ordered[, 1L] - ordered[, 2L]
  margin[!is.finite(margin)] <- Inf
  assignment <- colnames(score_matrix)[best_index]
  assignment[
    !is.finite(best_score) |
      best_score <= 0 |
      margin < setting$margin_gate
  ] <- "Uncertain"
  assignment <- ifelse(
    assignment == "Uncertain",
    "Uncertain",
    unname(niche_labels[assignment])
  )
  current_assignment <- data.frame(
    setting_id = setting$setting_id,
    center_threshold = setting$center_threshold,
    radius_multiplier = setting$radius_multiplier,
    program_gate = setting$program_gate,
    margin_gate = setting$margin_gate,
    sample_key = current$sample_key,
    condition = current$condition,
    timepoint = current$timepoint,
    analysis_role = current$analysis_role,
    center_barcode = current$center_barcode,
    center_macrophage_prop = current$center_macrophage_prop,
    center_total_myeloid_prop = current$center_total_myeloid_prop,
    niche = assignment,
    niche_score = ifelse(is.finite(best_score), best_score, NA_real_),
    niche_margin = ifelse(is.finite(margin), margin, NA_real_),
    score_inflammatory_infarct = current$score_inflammatory_infarct,
    score_macrophage_fibroblast_remodeling =
      current$score_macrophage_fibroblast_remodeling,
    score_vascular_repair = current$score_vascular_repair,
    score_remote_resident = current$score_remote_resident,
    stringsAsFactors = FALSE
  )
  assignment_rows[[setting$setting_id]] <- current_assignment
  count <- as.data.frame(table(
    sample_key = current_assignment$sample_key,
    niche = current_assignment$niche
  ), stringsAsFactors = FALSE)
  count <- count[count$Freq > 0L, ]
  total <- aggregate(Freq ~ sample_key, count, sum)
  names(total)[2L] <- "centers"
  count <- merge(count, total, by = "sample_key", all.x = TRUE)
  count$fraction <- count$Freq / count$centers
  count$setting_id <- setting$setting_id
  summary_rows[[setting$setting_id]] <- count
}
assignments <- do.call(rbind, assignment_rows)
rownames(assignments) <- NULL
grid_summary <- do.call(rbind, summary_rows)
rownames(grid_summary) <- NULL
write_tsv_gz(
  assignments,
  file.path(root, "metadata", "phase_09_niche_grid_assignments.tsv.gz")
)
write_tsv(
  grid_summary,
  file.path(root, "metadata", "phase_09_niche_grid_slice_summary.tsv")
)

baseline <- assignments[assignments$setting_id == baseline_id, , drop = FALSE]
baseline_key <- paste(baseline$sample_key, baseline$center_barcode, sep = "::")
stability_rows <- list()
for (setting_id in grid$setting_id) {
  current <- assignments[assignments$setting_id == setting_id, , drop = FALSE]
  current_key <- paste(current$sample_key, current$center_barcode, sep = "::")
  shared <- intersect(baseline_key, current_key)
  baseline_label <- baseline$niche[match(shared, baseline_key)]
  current_label <- current$niche[match(shared, current_key)]
  row <- data.frame(
    setting_id = setting_id,
    shared_centers = length(shared),
    fraction_of_baseline_centers_shared =
      length(shared) / length(baseline_key),
    adjusted_rand_vs_baseline = adjusted_rand(baseline_label, current_label),
    assigned_concordance = mean(baseline_label == current_label),
    stringsAsFactors = FALSE
  )
  for (label in unname(niche_labels)) {
    safe <- gsub("[^A-Za-z0-9]+", "_", label)
    row[[paste0("jaccard_", safe)]] <- jaccard(
      baseline_label, current_label, label
    )
  }
  stability_rows[[setting_id]] <- row
}
stability <- merge(
  grid,
  do.call(rbind, stability_rows),
  by = "setting_id", all.x = TRUE
)
write_tsv(
  stability,
  file.path(root, "metadata", "phase_09_niche_grid_stability.tsv")
)

parameter_summary <- do.call(rbind, lapply(grid$setting_id, function(setting_id) {
  current <- assignments[assignments$setting_id == setting_id, , drop = FALSE]
  slice_niches <- unique(current[, c("sample_key", "niche")])
  data.frame(
    setting_id = setting_id,
    centers = nrow(current),
    assigned_fraction = mean(current$niche != "Uncertain"),
    inflammatory_slices = sum(
      slice_niches$niche == niche_labels[["inflammatory_infarct"]]
    ),
    remodeling_slices = sum(
      slice_niches$niche ==
        niche_labels[["macrophage_fibroblast_remodeling"]]
    ),
    vascular_slices = sum(
      slice_niches$niche == niche_labels[["vascular_repair"]]
    ),
    remote_slices = sum(
      slice_niches$niche == niche_labels[["remote_resident"]]
    ),
    stringsAsFactors = FALSE
  )
}))
parameter_summary <- merge(grid, parameter_summary, by = "setting_id")
write_tsv(
  parameter_summary,
  file.path(root, "metadata", "phase_09_niche_parameter_summary.tsv")
)

# Baseline spatial null: labels are permuted within each physical slice while
# the observed center graph and label counts remain fixed.
permutation_rows <- list()
for (sample in samples) {
  current <- baseline[baseline$sample_key == sample, , drop = FALSE]
  current_key <- current$center_barcode
  current_edges <- edges[
    edges$sample_key == sample &
      abs(edges$radius_multiplier - 1.2) < 1e-8 &
      !edges$is_center &
      edges$center_barcode %in% current_key &
      edges$neighbor_barcode %in% current_key, ,
    drop = FALSE
  ]
  if (nrow(current_edges)) {
    pair <- t(apply(
      current_edges[, c("center_barcode", "neighbor_barcode")],
      1L, sort
    ))
    current_edges <- current_edges[!duplicated(
      paste(pair[, 1L], pair[, 2L], sep = "::")
    ), , drop = FALSE]
  }
  labels <- setNames(current$niche, current$center_barcode)
  if (!nrow(current_edges)) {
    permutation_rows[[sample]] <- data.frame(
      sample_key = sample,
      centers = nrow(current),
      undirected_center_edges = 0L,
      observed_same_label_fraction = NA_real_,
      null_mean = NA_real_,
      null_sd = NA_real_,
      enrichment = NA_real_,
      empirical_p = NA_real_,
      permutations = 0L,
      stringsAsFactors = FALSE
    )
    next
  }
  observed <- mean(
    labels[current_edges$center_barcode] ==
      labels[current_edges$neighbor_barcode]
  )
  set.seed(20260731L + match(sample, samples))
  null <- replicate(1000L, {
    permuted <- sample(labels, length(labels), replace = FALSE)
    names(permuted) <- names(labels)
    mean(
      permuted[current_edges$center_barcode] ==
        permuted[current_edges$neighbor_barcode]
    )
  })
  permutation_rows[[sample]] <- data.frame(
    sample_key = sample,
    centers = nrow(current),
    undirected_center_edges = nrow(current_edges),
    observed_same_label_fraction = observed,
    null_mean = mean(null),
    null_sd = sd(null),
    enrichment = observed - mean(null),
    empirical_p = (1 + sum(null >= observed)) / 1001,
    permutations = 1000L,
    stringsAsFactors = FALSE
  )
}
spatial_null <- do.call(rbind, permutation_rows)
write_tsv(
  spatial_null,
  file.path(root, "metadata", "phase_09_niche_spatial_null.tsv")
)

# Slice-level continuous scores and leave-one-slice-out effect stability.
continuous_columns <- c(
  inflammatory_axis = "score_inflammatory_infarct",
  remodeling_axis = "score_macrophage_fibroblast_remodeling",
  vascular_axis = "score_vascular_repair",
  remote_resident_axis = "score_remote_resident"
)
slice_continuous <- do.call(rbind, lapply(
  split(baseline, baseline$sample_key),
  function(x) data.frame(
    sample_key = x$sample_key[[1L]],
    condition = x$condition[[1L]],
    timepoint = x$timepoint[[1L]],
    analysis_role = x$analysis_role[[1L]],
    centers = nrow(x),
    inflammatory_axis = mean(x$score_inflammatory_infarct),
    remodeling_axis = mean(x$score_macrophage_fibroblast_remodeling),
    vascular_axis = mean(x$score_vascular_repair),
    remote_resident_axis = mean(x$score_remote_resident),
    stringsAsFactors = FALSE
  )
))
write_tsv(
  slice_continuous,
  file.path(root, "metadata", "phase_09_niche_continuous_slice_summary.tsv")
)

slice_niche <- as.data.frame(table(
  sample_key = baseline$sample_key,
  niche = baseline$niche
), stringsAsFactors = FALSE)
slice_niche <- slice_niche[slice_niche$Freq > 0L, ]
slice_total <- aggregate(Freq ~ sample_key, slice_niche, sum)
names(slice_total)[2L] <- "centers"
slice_niche <- merge(slice_niche, slice_total, by = "sample_key", all.x = TRUE)
slice_niche$fraction <- slice_niche$Freq / slice_niche$centers

contrast_definitions <- list(
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

metric_rows <- list()
for (axis in names(continuous_columns)) {
  metric_rows[[length(metric_rows) + 1L]] <- data.frame(
    sample_key = slice_continuous$sample_key,
    metric = axis,
    value = slice_continuous[[axis]],
    stringsAsFactors = FALSE
  )
}
for (label in unname(niche_labels)) {
  table <- data.frame(
    sample_key = primary_samples,
    metric = paste0(
      "fraction_",
      gsub("[^A-Za-z0-9]+", "_", label)
    ),
    value = 0,
    stringsAsFactors = FALSE
  )
  match_rows <- slice_niche$niche == label &
    slice_niche$sample_key %in% primary_samples
  matched_index <- match(
    slice_niche$sample_key[match_rows], table$sample_key
  )
  table$value[matched_index] <- slice_niche$fraction[match_rows]
  metric_rows[[length(metric_rows) + 1L]] <- table
}
metric_table <- do.call(rbind, metric_rows)
loso_rows <- list()
for (contrast_name in names(contrast_definitions)) {
  definition <- contrast_definitions[[contrast_name]]
  for (metric in unique(metric_table$metric)) {
    values <- metric_table[metric_table$metric == metric, , drop = FALSE]
    numerator <- setNames(
      values$value[match(definition$numerator, values$sample_key)],
      definition$numerator
    )
    denominator <- setNames(
      values$value[match(definition$denominator, values$sample_key)],
      definition$denominator
    )
    full_effect <- exact_effect(numerator, denominator)
    drop_samples <- c(names(numerator), names(denominator))
    loo_effects <- vapply(drop_samples, function(drop_sample) {
      numerator_loo <- numerator[names(numerator) != drop_sample]
      denominator_loo <- denominator[names(denominator) != drop_sample]
      if (!length(numerator_loo) || !length(denominator_loo)) return(NA_real_)
      exact_effect(numerator_loo, denominator_loo)
    }, numeric(1L))
    sign_match <- sign(loo_effects) == sign(full_effect)
    loso_rows[[length(loso_rows) + 1L]] <- data.frame(
      contrast = contrast_name,
      metric = metric,
      full_effect = full_effect,
      leave_one_out_runs = sum(is.finite(loo_effects)),
      min_leave_one_out_effect = min(loo_effects, na.rm = TRUE),
      max_leave_one_out_effect = max(loo_effects, na.rm = TRUE),
      sign_stability_fraction = mean(sign_match[is.finite(loo_effects)]),
      stringsAsFactors = FALSE
    )
  }
}
loso <- do.call(rbind, loso_rows)
write_tsv(
  loso,
  file.path(root, "metadata", "phase_09_niche_leave_one_slice_out.tsv")
)

# Independent unsupervised recovery. K=2..8 is evaluated; K=4 is not forced.
unsupervised_features <- features[
  features$radius_multiplier == 1.2 &
    features$center_total_myeloid_prop >= 0.20, ,
  drop = FALSE
]
unsupervised_key <- paste(
  unsupervised_features$sample_key,
  unsupervised_features$center_barcode,
  sep = "::"
)
unsupervised_columns <- c(
  "local_cell_Cardiomyocyte", "local_cell_Endothelial",
  "local_cell_Fibroblast", "local_cell_Macrophage",
  "local_cell_Monocyte", "local_cell_Neutrophil",
  "local_state_Resident_Folr2_macrophage",
  "local_state_Trem2_Spp1_macrophage",
  "local_state_MHCII_macrophage",
  "module_inflammation", "module_ECM", "module_vascular",
  "module_cardiomyocyte_stress"
)
unsupervised_matrix <- as.matrix(
  unsupervised_features[, unsupervised_columns, drop = FALSE]
)
unsupervised_matrix <- scale(unsupervised_matrix)
unsupervised_matrix[!is.finite(unsupervised_matrix)] <- 0
rule_labels <- baseline$niche[match(unsupervised_key, baseline_key)]
clara_rows <- list()
clara_objects <- list()
for (k in 2:8) {
  set.seed(20260731L + k)
  fit <- clara(
    unsupervised_matrix,
    k = k,
    metric = "euclidean",
    samples = 100L,
    sampsize = min(1000L, nrow(unsupervised_matrix)),
    pamLike = TRUE,
    rngR = TRUE
  )
  clara_objects[[as.character(k)]] <- fit
  assigned <- rule_labels != "Uncertain"
  clara_rows[[as.character(k)]] <- data.frame(
    k = k,
    average_silhouette = fit$silinfo$avg.width,
    ari_all_vs_rule_labels = adjusted_rand(rule_labels, fit$clustering),
    ari_assigned_only = adjusted_rand(
      rule_labels[assigned], fit$clustering[assigned]
    ),
    centers = nrow(unsupervised_matrix),
    stringsAsFactors = FALSE
  )
}
unsupervised_summary <- do.call(rbind, clara_rows)
best_k <- unsupervised_summary$k[
  which.max(unsupervised_summary$average_silhouette)
]
unsupervised_summary$selected_by_silhouette <-
  unsupervised_summary$k == best_k
write_tsv(
  unsupervised_summary,
  file.path(root, "metadata", "phase_09_niche_unsupervised_recovery.tsv")
)
best_fit <- clara_objects[[as.character(best_k)]]
unsupervised_assignment <- data.frame(
  sample_key = unsupervised_features$sample_key,
  center_barcode = unsupervised_features$center_barcode,
  rule_based_niche = rule_labels,
  selected_k = best_k,
  unsupervised_cluster = best_fit$clustering,
  stringsAsFactors = FALSE
)
write_tsv_gz(
  unsupervised_assignment,
  file.path(root, "metadata", "phase_09_niche_unsupervised_assignments.tsv.gz")
)

# Save baseline labels to new objects and render native scop maps.
niche_maps <- character(length(samples))
axis_maps <- character(length(samples))
object_audit_rows <- list()
for (i in seq_along(samples)) {
  sample <- samples[[i]]
  input_path <- file.path(
    root, "objects", "phase_09", "macrophage_mapping",
    paste0(sample, "_phase09_macrophage.rds")
  )
  object <- readRDS(input_path)
  current <- baseline[baseline$sample_key == sample, , drop = FALSE]
  metadata <- data.frame(
    P09N_center = colnames(object) %in% current$center_barcode,
    P09N_niche = "Not_macrophage_center",
    P09N_niche_score = NA_real_,
    P09N_niche_margin = NA_real_,
    P09N_inflammatory_axis = NA_real_,
    P09N_remodeling_axis = NA_real_,
    P09N_vascular_axis = NA_real_,
    P09N_remote_axis = NA_real_,
    row.names = colnames(object),
    stringsAsFactors = FALSE
  )
  metadata[current$center_barcode, "P09N_niche"] <- current$niche
  metadata[current$center_barcode, "P09N_niche_score"] <- current$niche_score
  metadata[current$center_barcode, "P09N_niche_margin"] <- current$niche_margin
  metadata[current$center_barcode, "P09N_inflammatory_axis"] <-
    current$score_inflammatory_infarct
  metadata[current$center_barcode, "P09N_remodeling_axis"] <-
    current$score_macrophage_fibroblast_remodeling
  metadata[current$center_barcode, "P09N_vascular_axis"] <-
    current$score_vascular_repair
  metadata[current$center_barcode, "P09N_remote_axis"] <-
    current$score_remote_resident
  object <- AddMetaData(object, metadata)
  object@misc$phase_09_niche <- list(
    baseline_setting_id = baseline_id,
    center_threshold = 0.20,
    center_threshold_definition =
      "RCTD Macrophage + Monocyte + Neutrophil + Dendritic",
    macrophage_state_mapping_eligibility =
      "P09E_RCTD_prop_Macrophage >= 0.10",
    radius_multiplier = 1.2,
    program_gate = 0.10,
    margin_gate = 0.15,
    grid_settings_tested = 108L,
    spatial_permutations = 1000L,
    spots_are_not_replicates = TRUE
  )
  output_path <- file.path(
    root, "objects", "phase_09", "niches",
    paste0(sample, "_phase09_niche.rds")
  )
  saveRDS(object, output_path, compress = FALSE)
  readback <- readRDS(output_path)
  if (sum(readback$P09N_center) != nrow(current)) {
    stop("Baseline niche object readback failed: ", sample)
  }
  image_name <- Images(readback)[[1L]]
  niche_maps[[i]] <- file.path(
    root, "results", "phase_09", "niche_robustness",
    paste0(sample, "_baseline_niche.png")
  )
  axis_maps[[i]] <- file.path(
    root, "results", "phase_09", "niche_robustness",
    paste0(sample, "_remodeling_axis.png")
  )
  niche_values <- as.character(readback$P09N_niche)
  p_niche <- SpatialSpotPlot(
    readback,
    values = setNames(niche_values, colnames(readback)),
    plot_type = "point",
    image = image_name,
    crop = FALSE,
    coord.cols = c("x", "y"),
    overlay_image = TRUE,
    image.alpha = 0.85,
    pt.alpha = 0.9,
    palcolor = niche_palette[unique(niche_values)],
    legend.title = "Corrected niche"
  ) + ggtitle(
    paste0(sample, ": corrected baseline niches"),
    subtitle = "Baseline is one prespecified member of a 108-setting grid"
  )
  ggsave(
    niche_maps[[i]], p_niche,
    width = 9, height = 8, dpi = 260, bg = "white"
  )
  axis_values <- as.numeric(readback$P09N_remodeling_axis)
  p_axis <- SpatialSpotPlot(
    readback,
    values = setNames(axis_values, colnames(readback)),
    plot_type = "point",
    image = image_name,
    crop = FALSE,
    coord.cols = c("x", "y"),
    overlay_image = TRUE,
    image.alpha = 0.85,
    pt.alpha = 0.9,
    palcolor = c("#2F6FB0", "#F7F7F7", "#C81D25"),
    legend.title = "Remodeling axis"
  ) + ggtitle(paste0(sample, ": continuous remodeling axis"))
  ggsave(
    axis_maps[[i]], p_axis,
    width = 9, height = 8, dpi = 260, bg = "white"
  )
  object_audit_rows[[sample]] <- data.frame(
    sample_key = sample,
    spots = ncol(readback),
    centers = sum(readback$P09N_center),
    images = length(Images(readback)),
    output_path = output_path,
    niche_map = niche_maps[[i]],
    axis_map = axis_maps[[i]],
    stringsAsFactors = FALSE
  )
  rm(object, readback)
  gc()
}
object_audit <- do.call(rbind, object_audit_rows)
write_tsv(
  object_audit,
  file.path(root, "metadata", "phase_09_niche_object_audit.tsv")
)
montage_status <- c(
  system2("montage", c(
    shQuote(niche_maps), "-tile", "3x4", "-geometry", "650x560+8+8",
    shQuote(file.path(
      root, "results", "phase_09", "niche_robustness", "visual_qc",
      "phase_09_niche_contact_sheet.jpg"
    ))
  )),
  system2("montage", c(
    shQuote(axis_maps), "-tile", "3x4", "-geometry", "650x560+8+8",
    shQuote(file.path(
      root, "results", "phase_09", "niche_robustness", "visual_qc",
      "phase_09_remodeling_axis_contact_sheet.jpg"
    ))
  ))
)
if (any(montage_status != 0L)) stop("Niche contact sheet failed")

validation <- data.frame(
  check = c(
    "parameter_grid_108_complete",
    "all_grid_settings_have_assignments",
    "baseline_setting_unique",
    "physical_slices_12",
    "cross_slice_edges_zero",
    "neighborhood_min_one",
    "total_myeloid_center_fraction_valid",
    "module_min_four_genes",
    "grid_stability_108",
    "spatial_null_12_slices",
    "spatial_null_1000_per_nonisolated_slice",
    "leave_one_slice_out_complete",
    "unsupervised_k2_to_k8_complete",
    "unsupervised_k_not_forced_to_4",
    "object_readback_12",
    "native_scop_niche_maps_12",
    "native_scop_axis_maps_12",
    "contact_sheets_2"
  ),
  value = c(
    nrow(grid) == 108L,
    length(unique(assignments$setting_id)) == 108L,
    length(baseline_id) == 1L,
    length(unique(features$sample_key)) == 12L,
    all(edges$sample_key == edges$neighbor_sample_key),
    min(features$neighborhood_size) >= 1L,
    all(
      is.finite(features$center_total_myeloid_prop) &
        features$center_total_myeloid_prop >=
          features$center_macrophage_prop &
        features$center_total_myeloid_prop <= 1 + 1e-8
    ),
    min(module_audit$genes_detected) >= 4L,
    nrow(stability) == 108L,
    nrow(spatial_null) == 12L,
    all(
      spatial_null$permutations == 1000L |
        spatial_null$undirected_center_edges == 0L
    ),
    nrow(loso) == 3L * (length(continuous_columns) + length(niche_labels)),
    nrow(unsupervised_summary) == 7L &&
      setequal(unsupervised_summary$k, 2:8),
    sum(unsupervised_summary$selected_by_silhouette) == 1L,
    nrow(object_audit) == 12L &&
      all(object_audit$images == 1L),
    sum(file.exists(niche_maps)) == 12L,
    sum(file.exists(axis_maps)) == 12L,
    all(file.exists(file.path(
      root, "results", "phase_09", "niche_robustness", "visual_qc",
      c(
        "phase_09_niche_contact_sheet.jpg",
        "phase_09_remodeling_axis_contact_sheet.jpg"
      )
    )))
  ),
  status = "",
  stringsAsFactors = FALSE
)
validation$status <- ifelse(as.logical(validation$value), "PASS", "FAIL")
write_tsv(
  validation,
  file.path(root, "metadata", "phase_09_niche_robustness_validation.tsv")
)
log_msg(
  "Niche robustness complete:",
  " settings=", nrow(grid),
  " baseline_centers=", nrow(baseline),
  " best_unsupervised_k=", best_k,
  " validation=", paste(table(validation$status), collapse = ",")
)
if (any(validation$status != "PASS")) {
  stop("Niche robustness validation failed")
}
