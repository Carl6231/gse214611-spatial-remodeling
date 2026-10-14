args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else {
  "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
}
.libPaths(c(file.path(root, "Rlib"), .libPaths()))

dir.create(file.path(root, "objects", "phase_09"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results", "phase_09", "reference_build"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "logs"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(root, "logs", "phase_09_build_references.log")
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
file_signature <- function(path) {
  info <- file.info(path)
  data.frame(
    path = normalizePath(path, mustWork = TRUE),
    bytes = as.numeric(info$size),
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )
}
select_pseudobulk_features <- function(
    counts, labels, required_labels, target_features, top_per_label,
    mandatory_markers, exclude_pattern = "^(mt-|Rpl|Rps)") {
  if (!identical(colnames(counts), names(labels))) {
    stop("Counts columns and label names are not identical")
  }
  if (!setequal(unique(labels), required_labels)) {
    stop("Reference labels differ from required labels")
  }
  pseudobulk <- vapply(required_labels, function(label) {
    Matrix::rowSums(counts[, labels == label, drop = FALSE])
  }, numeric(nrow(counts)))
  rownames(pseudobulk) <- rownames(counts)
  colnames(pseudobulk) <- required_labels
  totals <- colSums(pseudobulk)
  if (any(totals <= 0)) stop("A reference class has zero total counts")
  log_cpm <- log1p(sweep(pseudobulk, 2L, totals, "/") * 1e6)
  variance <- apply(log_cpm, 1L, stats::var)
  specificity <- vapply(seq_along(required_labels), function(i) {
    if (length(required_labels) == 1L) return(log_cpm[, i])
    log_cpm[, i] - rowMeans(log_cpm[, -i, drop = FALSE])
  }, numeric(nrow(counts)))
  rownames(specificity) <- rownames(counts)
  colnames(specificity) <- required_labels
  eligible <- !grepl(exclude_pattern, rownames(counts), ignore.case = TRUE, perl = TRUE)
  top_by_label <- lapply(required_labels, function(label) {
    score <- specificity[, label]
    candidates <- which(eligible & is.finite(score) & score > 0)
    rownames(counts)[candidates[order(score[candidates], decreasing = TRUE)][
      seq_len(min(top_per_label, length(candidates)))
    ]]
  })
  names(top_by_label) <- required_labels
  mandatory <- intersect(unique(mandatory_markers), rownames(counts))
  selected <- unique(c(mandatory, unlist(top_by_label, use.names = FALSE)))
  fill <- rownames(counts)[order(variance, decreasing = TRUE)]
  fill <- fill[eligible[match(fill, rownames(counts))]]
  selected <- unique(c(selected, fill))
  if (length(selected) < target_features) {
    stop("Unable to select target features")
  }
  selected <- selected[seq_len(target_features)]
  if (!all(mandatory %in% selected)) stop("Mandatory marker preservation failed")

  reason <- setNames(rep("variance_fill", nrow(counts)), rownames(counts))
  for (label in required_labels) {
    reason[top_by_label[[label]]] <- paste0("top_specific_", label)
  }
  reason[mandatory] <- "mandatory_marker"
  table <- data.frame(
    feature = rownames(counts),
    selected = rownames(counts) %in% selected,
    selected_rank = match(rownames(counts), selected),
    selection_reason = unname(reason[rownames(counts)]),
    mandatory_marker = rownames(counts) %in% mandatory,
    profile_variance = variance,
    max_specificity = apply(specificity, 1L, max),
    dominant_reference_profile = required_labels[
      max.col(log_cpm, ties.method = "first")
    ],
    stringsAsFactors = FALSE
  )
  for (label in required_labels) {
    safe <- gsub("[^A-Za-z0-9]+", "_", label)
    table[[paste0("logCPM_", safe)]] <- log_cpm[, label]
    table[[paste0("specificity_", safe)]] <- specificity[, label]
  }
  list(selected = selected, table = table, log_cpm = log_cpm)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

source_srt <- "/home/zhenhanli/myocardial/scop/srt.rds"
phase6_path <- file.path(
  root, "objects", "macrophage", "reference_macrophage_clustered.rds"
)
spatial_example_path <- file.path(
  root, "objects", "qc", "Sham_Day3_rep1_qc.rds"
)
for (path in c(source_srt, phase6_path, spatial_example_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}
source_signature_before <- file_signature(source_srt)

cluster_annotation <- data.frame(
  P6_cluster = as.character(0:11),
  audited_state = c(
    "Fibroblast_contaminant",
    "Resident_Folr2_macrophage",
    "Trem2_Spp1_macrophage",
    "MHCII_macrophage",
    "Inflammatory_monocyte",
    "Endothelial_contaminant",
    "Neutrophil",
    "Dendritic_cell",
    "Cycling_macrophage",
    "B_cell",
    "T_NK",
    "Stress_response_macrophage"
  ),
  expanded_celltype = c(
    "Fibroblast", "Macrophage", "Macrophage", "Macrophage",
    "Monocyte", "Endothelial", "Neutrophil", "Dendritic",
    "Macrophage", "B_cell", "T_NK", "Macrophage"
  ),
  retain_in_refined_macrophage = c(
    FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
    TRUE, FALSE, FALSE, TRUE
  ),
  primary_support = c(
    "Col1a1/Col1a2/Col3a1/Pdgfra/Postn/Col11a1/Cdh11; weak macrophage core",
    "Lyve1/Folr2/Timd4/Vsig4/Cd163/Mrc1 plus macrophage core",
    "Trem2/Spp1/Gpnmb/Lgals3/Fabp5/Ctsd plus macrophage core",
    "Cd74/H2-Aa/H2-Ab1/H2-Eb1 plus 10/10 macrophage-core markers",
    "Ccr2/Vcan/Ly6c2/Chil3/Plac8/Clec4e",
    "Pecam1/Kdr/Emcn/Egfl7/Robo4/Flt1/Cdh5/Esam; weak macrophage core",
    "S100a8/S100a9/Csf3r/Retnlg/Mpo/Elane/Camp/Ngp",
    "Flt3/Clec9a/Itgax/Itgae/Zbtb46/Xcr1 plus MHC-II",
    "Mki67/Top2a/Birc5/Pclaf/Cenpf/Ube2c/Foxm1 plus 10/10 macrophage-core markers",
    "Cd79a/Cd79b/Cd19/Ms4a1/Cd22/Iglc2",
    "Cd3d/Cd3e/Trbc1/Trbc2/Tcf7/Bcl11b/Nkg7/Klrd1",
    "macrophage core with Hspa1a/Hspa1b/Rgs1/Ccl3/Ccl4/Nr4a1 stress program"
  ),
  ambiguity_note = c(
    "Cardiomyocyte transcripts are treated as ambient myocardial RNA because the defining cluster markers are stromal and Ptprc/macrophage-core expression is weak.",
    "None",
    "None",
    "MHC-II expression alone is not used to call dendritic lineage; macrophage-core expression and absence of Clec9a/Flt3 support macrophage.",
    "Monocyte-to-macrophage continuum; retained as a separate broad class for sensitivity analysis.",
    "None",
    "Some myeloid-core expression is expected; granulocyte markers define the broad class.",
    "MHC-II overlaps macrophages, but Clec9a/Flt3/Zbtb46/Xcr1 support dendritic lineage.",
    "Cell-cycle state is separated from lineage; macrophage-core and Trem2/Spp1 markers support macrophage lineage.",
    "Cardiomyocyte transcripts are treated as ambient myocardial RNA.",
    "Cardiomyocyte transcripts are treated as ambient myocardial RNA.",
    "Small cluster (n=74); retained only as a stress-response macrophage state and reported cautiously."
  ),
  evidence_source = "Phase 6 DE markers + Phase 9 pseudobulk CPM/detection + native scop GroupHeatmap",
  stringsAsFactors = FALSE
)
write_tsv(
  cluster_annotation,
  file.path(root, "metadata", "phase_09_cluster_annotation.tsv")
)

log_msg("Loading Phase 6 clustered object")
mac <- readRDS(phase6_path)
if (!inherits(mac, "Seurat") ||
    !"RNA" %in% Assays(mac) ||
    !"counts" %in% Layers(mac[["RNA"]])) {
  stop("Phase 6 clustered object is incompatible")
}
cluster_by_cell <- setNames(as.character(mac$P6_cluster), colnames(mac))
state_by_cluster <- setNames(
  cluster_annotation$audited_state, cluster_annotation$P6_cluster
)
retain_cluster <- cluster_annotation$P6_cluster[
  cluster_annotation$retain_in_refined_macrophage
]
retained_cells <- colnames(mac)[as.character(mac$P6_cluster) %in% retain_cluster]
if (length(retained_cells) < 8000L) stop("Too few refined macrophage cells")

log_msg("Building refined macrophage-state reference")
mac_counts <- LayerData(mac, assay = "RNA", layer = "counts")
spatial_example <- readRDS(spatial_example_path)
shared_mac <- intersect(rownames(mac_counts), rownames(spatial_example))
rm(spatial_example)
if (length(shared_mac) < 15000L) stop("Too few macrophage-spatial shared genes")
mac_counts <- mac_counts[shared_mac, retained_cells, drop = FALSE]
mac_state <- unname(state_by_cluster[cluster_by_cell[retained_cells]])
names(mac_state) <- retained_cells
required_states <- cluster_annotation$audited_state[
  cluster_annotation$retain_in_refined_macrophage
]
mac_markers <- c(
  "Lyve1", "Folr2", "Timd4", "Vsig4", "Cd163", "Mrc1",
  "Trem2", "Spp1", "Gpnmb", "Lgals3", "Fabp5", "Ctsd", "Apoe",
  "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Ciita",
  "Mki67", "Top2a", "Birc5", "Pclaf", "Cenpf", "Ube2c", "Foxm1",
  "Hspa1a", "Hspa1b", "Rgs1", "Ccl3", "Ccl4", "Nr4a1",
  "Ptprc", "Lyz2", "Csf1r", "Adgre1", "Fcgr1", "C1qa", "C1qb", "C1qc"
)
mac_selection <- select_pseudobulk_features(
  counts = mac_counts,
  labels = mac_state,
  required_labels = required_states,
  target_features = 1800L,
  top_per_label = 300L,
  mandatory_markers = mac_markers
)
write_tsv(
  mac_selection$table,
  file.path(root, "metadata", "phase_09_refined_macrophage_feature_selection.tsv")
)
mac_meta <- mac[[]][retained_cells, c(
  "sample_id", "timepoint", "P6_cluster"
), drop = FALSE]
mac_meta$P09_macrophage_state <- factor(
  mac_state[rownames(mac_meta)], levels = required_states
)
mac_reference <- CreateSeuratObject(
  counts = mac_counts[mac_selection$selected, , drop = FALSE],
  assay = "RNA",
  meta.data = mac_meta,
  min.cells = 0,
  min.features = 0,
  project = "GSE214611_Phase9_refined_macrophage_reference"
)
mac_reference@misc$phase_09_reference <- list(
  source_path = phase6_path,
  retained_clusters = retain_cluster,
  excluded_clusters = setdiff(cluster_annotation$P6_cluster, retain_cluster),
  selected_features = mac_selection$selected,
  annotation_table = cluster_annotation
)
mac_reference_path <- file.path(
  root, "objects", "phase_09", "reference_macrophage_refined.rds"
)
saveRDS(mac_reference, mac_reference_path, compress = FALSE)
mac_readback <- readRDS(mac_reference_path)
if (!inherits(mac_readback, "Seurat") ||
    ncol(mac_readback) != length(retained_cells) ||
    nrow(mac_readback) != 1800L) {
  stop("Refined macrophage reference readback failed")
}
mac_reference_validation <- c(
  readable = inherits(mac_readback, "Seurat"),
  features_1800 = nrow(mac_readback) == 1800L,
  states_5 = length(unique(as.character(
    mac_readback$P09_macrophage_state
  ))) == 5L
)

state_composition <- as.data.frame(table(
  state = as.character(mac_reference$P09_macrophage_state),
  sample_id = as.character(mac_reference$sample_id),
  timepoint = as.character(mac_reference$timepoint)
), stringsAsFactors = FALSE)
state_composition <- state_composition[state_composition$Freq > 0L, ]
state_total <- aggregate(Freq ~ state, state_composition, sum)
names(state_total)[2L] <- "state_cells"
state_composition <- merge(state_composition, state_total, by = "state", all.x = TRUE)
state_composition$fraction_within_state <-
  state_composition$Freq / state_composition$state_cells
write_tsv(
  state_composition,
  file.path(root, "metadata", "phase_09_refined_macrophage_composition.tsv")
)
rm(mac_reference, mac_readback, mac_counts)
gc()

log_msg("Loading read-only 34 GB srt.rds for expanded all-cell reference")
srt <- readRDS(source_srt)
if (!inherits(srt, "Seurat") ||
    !"RNA" %in% Assays(srt) ||
    !all(c("celltype", "sample_id", "timepoint") %in% colnames(srt[[]]))) {
  stop("Full reference srt is incompatible")
}
if (!all(names(cluster_by_cell) %in% colnames(srt))) {
  stop("Phase 6 cell names are not all present in full srt")
}
if (!all(as.character(srt$celltype[names(cluster_by_cell)]) == "Macrophage")) {
  stop("Phase 6 cells no longer map to original Macrophage labels")
}

log_msg("Joining full RNA count layers in memory")
rna_joined <- JoinLayers(srt[["RNA"]])
if (!"counts" %in% Layers(rna_joined)) stop("Joined RNA counts are absent")
counts <- LayerData(rna_joined, layer = "counts", fast = FALSE)
if (!inherits(counts, "sparseMatrix") ||
    !identical(colnames(counts), colnames(srt))) {
  stop("Joined counts do not match the full reference")
}
spatial_example <- readRDS(spatial_example_path)
shared <- intersect(rownames(counts), rownames(spatial_example))
rm(spatial_example)
if (length(shared) < 15000L) stop("Too few all-cell-spatial shared genes")
counts_shared <- counts[shared, , drop = FALSE]

meta <- srt[[]][colnames(counts_shared), c(
  "celltype", "sample_id", "timepoint"
), drop = FALSE]
meta$P09_source_cluster <- NA_character_
meta$P09_audited_state <- NA_character_
mac_positions <- match(names(cluster_by_cell), rownames(meta))
meta$P09_source_cluster[mac_positions] <- unname(cluster_by_cell)
meta$P09_audited_state[mac_positions] <- unname(
  state_by_cluster[cluster_by_cell]
)
meta$P09_expanded_celltype <- as.character(meta$celltype)
expanded_by_cluster <- setNames(
  cluster_annotation$expanded_celltype, cluster_annotation$P6_cluster
)
meta$P09_expanded_celltype[mac_positions] <- unname(
  expanded_by_cluster[cluster_by_cell]
)

required_expanded <- c(
  "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage",
  "Monocyte", "Neutrophil", "Dendritic", "B_cell", "T_NK",
  "Pericyte", "Schwann_cell"
)
expanded_labels <- as.character(meta$P09_expanded_celltype)
names(expanded_labels) <- rownames(meta)
expanded_counts <- table(factor(expanded_labels, levels = required_expanded))
if (any(expanded_counts < 25L)) {
  stop("An expanded reference class has fewer than 25 cells: ",
       paste(names(expanded_counts)[expanded_counts < 25L], collapse = ", "))
}

expanded_markers <- c(
  "Tnnt2", "Myh6", "Myh7", "Actc1", "Ttn", "Ryr2",
  "Pecam1", "Cdh5", "Kdr", "Emcn", "Egfl7", "Robo4",
  "Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Postn",
  "Ptprc", "Lyz2", "Csf1r", "Adgre1", "Fcgr1", "C1qa", "C1qb",
  "Ccr2", "Vcan", "Ly6c2", "Chil3", "Plac8", "Clec4e",
  "S100a8", "S100a9", "Csf3r", "Retnlg", "Mpo", "Elane", "Camp", "Ngp",
  "Flt3", "Clec9a", "Itgax", "Itgae", "Zbtb46", "Xcr1",
  "Cd79a", "Cd79b", "Cd19", "Ms4a1", "Cd22", "Iglc2",
  "Cd3d", "Cd3e", "Trbc1", "Trbc2", "Tcf7", "Bcl11b", "Nkg7", "Klrd1",
  "Rgs5", "Pdgfrb", "Cspg4", "Kcnj8", "Acta2", "Tagln",
  "S100b", "Plp1", "Sox10", "Mpz"
)
log_msg(
  "Selecting expanded-reference features from ", length(shared),
  " shared genes and ", ncol(counts_shared), " cells"
)
expanded_selection <- select_pseudobulk_features(
  counts = counts_shared,
  labels = expanded_labels,
  required_labels = required_expanded,
  target_features = 3000L,
  top_per_label = 250L,
  mandatory_markers = expanded_markers
)
write_tsv(
  expanded_selection$table,
  file.path(root, "metadata", "phase_09_expanded_feature_selection.tsv")
)

reference_counts <- counts_shared[
  expanded_selection$selected, , drop = FALSE
]
nonzero <- Matrix::colSums(reference_counts) > 0
nonzero_by_type <- table(factor(
  expanded_labels[nonzero], levels = required_expanded
))
if (any(nonzero_by_type < 25L)) {
  stop("Selected features leave an expanded class below 25 nonzero cells")
}

log_msg(
  "Creating expanded compact reference: cells=", ncol(reference_counts),
  " features=", nrow(reference_counts),
  " classes=", length(required_expanded)
)
expanded_reference <- CreateSeuratObject(
  counts = reference_counts,
  assay = "RNA",
  meta.data = meta,
  min.cells = 0,
  min.features = 0,
  project = "GSE214611_Phase9_expanded_immune_reference"
)
expanded_reference$P09_expanded_celltype <- factor(
  as.character(expanded_reference$P09_expanded_celltype),
  levels = required_expanded
)
expanded_reference@misc$phase_09_reference <- list(
  source_path = source_srt,
  source_signature = source_signature_before,
  all_source_cells_retained = TRUE,
  source_cells = ncol(srt),
  selected_features = expanded_selection$selected,
  expanded_labels = required_expanded,
  cluster_annotation = cluster_annotation
)
expanded_reference_path <- file.path(
  root, "objects", "phase_09", "reference_expanded_immune.rds"
)
saveRDS(expanded_reference, expanded_reference_path, compress = FALSE)
expanded_readback <- readRDS(expanded_reference_path)
if (!inherits(expanded_readback, "Seurat") ||
    ncol(expanded_readback) != ncol(srt) ||
    nrow(expanded_readback) != 3000L ||
    !identical(colnames(expanded_readback), colnames(srt)) ||
    !setequal(
      unique(as.character(expanded_readback$P09_expanded_celltype)),
      required_expanded
    )) {
  stop("Expanded reference readback failed")
}

reference_summary <- data.frame(
  celltype = required_expanded,
  source_cells = as.integer(expanded_counts[required_expanded]),
  selected_feature_nonzero_cells = as.integer(
    nonzero_by_type[required_expanded]
  ),
  retained_fraction = as.integer(nonzero_by_type[required_expanded]) /
    as.integer(expanded_counts[required_expanded]),
  stringsAsFactors = FALSE
)
write_tsv(
  reference_summary,
  file.path(root, "metadata", "phase_09_expanded_reference_summary.tsv")
)

original_vs_expanded <- as.data.frame(table(
  original_celltype = as.character(meta$celltype),
  expanded_celltype = expanded_labels
), stringsAsFactors = FALSE)
original_vs_expanded <- original_vs_expanded[original_vs_expanded$Freq > 0L, ]
write_tsv(
  original_vs_expanded,
  file.path(root, "metadata", "phase_09_original_vs_expanded_reference.tsv")
)

source_signature_after <- file_signature(source_srt)
source_unchanged <- identical(source_signature_before, source_signature_after)
validation <- data.frame(
  check = c(
    "source_srt_unchanged",
    "expanded_reference_readback",
    "expanded_reference_all_source_cells",
    "expanded_reference_features_3000",
    "expanded_reference_classes_11",
    "expanded_class_min_cells_25",
    "refined_macrophage_reference_readback",
    "refined_macrophage_features_1800",
    "refined_macrophage_states_5",
    "nonmacrophage_clusters_excluded_from_refined_reference",
    "all_expanded_mandatory_markers_retained",
    "all_macrophage_mandatory_markers_retained"
  ),
  value = c(
    source_unchanged,
    inherits(expanded_readback, "Seurat"),
    ncol(expanded_readback) == ncol(srt),
    nrow(expanded_readback) == 3000L,
    length(unique(as.character(expanded_readback$P09_expanded_celltype))) == 11L,
    min(expanded_counts) >= 25L,
    mac_reference_validation[["readable"]],
    mac_reference_validation[["features_1800"]],
    mac_reference_validation[["states_5"]],
    all(!cluster_annotation$retain_in_refined_macrophage[
      cluster_annotation$expanded_celltype %in%
        c("Fibroblast", "Endothelial", "Monocyte", "Neutrophil",
          "Dendritic", "B_cell", "T_NK")
    ]),
    all(intersect(expanded_markers, shared) %in% expanded_selection$selected),
    all(intersect(mac_markers, shared_mac) %in% mac_selection$selected)
  ),
  status = "",
  stringsAsFactors = FALSE
)
validation$status <- ifelse(as.logical(validation$value), "PASS", "FAIL")
write_tsv(
  validation,
  file.path(root, "metadata", "phase_09_reference_build_validation.tsv")
)
capture.output(
  sessionInfo(),
  file = file.path(root, "results", "phase_09", "reference_build",
                   "session_info.txt")
)
log_msg(
  "Reference build complete:",
  "expanded_cells=", ncol(expanded_reference),
  "expanded_features=", nrow(expanded_reference),
  "refined_macrophage_cells=", length(retained_cells),
  "validation=", paste(table(validation$status), collapse = ",")
)
if (any(validation$status != "PASS")) {
  stop("Phase 9 reference build validation failed")
}
