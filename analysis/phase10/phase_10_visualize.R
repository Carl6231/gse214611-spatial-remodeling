root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
source(file.path(root, "scripts", "phase_10_common.R"))
root <- phase10_root()
dirs <- phase10_dirs(root)
phase10_log_start(root, "phase_10_visualize")

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scop)
})

set.seed(20260805)
figure_root <- file.path(dirs$results, "figures")
panel_dir <- file.path(figure_root, "panels")
final_dir <- file.path(figure_root, "final")
source_dir <- file.path(figure_root, "source_data")
qc_dir <- file.path(figure_root, "visual_qc")
invisible(lapply(c(panel_dir, final_dir, source_dir, qc_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

signature <- fread(require_file(file.path(dirs$metadata, "phase_10_day7_remodeling_signature.tsv")))
gsea <- fread(require_file(file.path(dirs$metadata, "phase_10_scop_gsea.tsv.gz")))
cell_scores <- fread(require_file(file.path(dirs$metadata, "phase_10_reference_celltype_signature_scores.tsv")))
cell_effects <- fread(require_file(file.path(dirs$metadata, "phase_10_reference_celltype_effects.tsv")))
macro_effects <- fread(require_file(file.path(dirs$metadata, "phase_10_macrophage_state_effects.tsv")))
quartile <- fread(require_file(file.path(dirs$metadata, "phase_10_day7_axis_decomposition_summary.tsv")))
quality <- fread(require_file(file.path(dirs$metadata, "phase_10_axis_quality_covariates.tsv")))
moran <- fread(require_file(file.path(dirs$metadata, "phase_10_signature_spatial_autocorrelation.tsv")))
external_annotation <- fread(require_file(file.path(dirs$metadata, "phase_10_external_annotation_projection.tsv.gz")))
external_lineage <- fread(require_file(file.path(dirs$metadata, "phase_10_external_lineage_projection.tsv")))
external_direction <- fread(require_file(file.path(dirs$metadata, "phase_10_external_signature_directions.tsv")))
spot_scores <- fread(require_file(file.path(dirs$metadata, "phase_10_spot_signature_scores.tsv.gz")))
map_manifest <- fread(require_file(file.path(dirs$metadata, "phase_10_spatial_map_manifest.tsv")))

direction_colors <- c(negative = "#2F6FB0", positive = "#C81D25")
time_colors <- c(Day0 = "#BDBDBD", `1HR` = "#F5C04A", `4HR` = "#E88E3D", Day1 = "#4DAF4A", Day3 = "#C95A49", Day7 = "#7E3F8F")
celltype_colors <- c(
  Cardiomyocyte = "#D73027", Endothelial = "#2C7BB6", Fibroblast = "#FDAE61",
  Macrophage = "#7B3294", Monocyte = "#C2A5CF", Neutrophil = "#E66101",
  Dendritic = "#5E3C99", B_cell = "#1B9E77", T_NK = "#66A61E",
  Pericyte = "#A6761D", Schwann_cell = "#666666"
)
theme_p10 <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 8, color = "#555555"),
        axis.text = element_text(color = "black"), legend.title = element_text(face = "bold"))

save_panel <- function(plot, stem, width = 7.2, height = 5.4) {
  png_path <- file.path(panel_dir, paste0(stem, ".png"))
  pdf_path <- file.path(panel_dir, paste0(stem, ".pdf"))
  if (file.exists(png_path) && file.exists(pdf_path) && file.info(png_path)$size > 10000 && file.info(pdf_path)$size > 1000) {
    return(c(png_path, pdf_path))
  }
  ggsave(png_path, plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  c(png_path, pdf_path)
}

pretty_term <- function(x) {
  x <- gsub("^(GOBP|HALLMARK|REACTOME)_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  tools::toTitleCase(tolower(x))
}

stamp("Creating 12 table-derived scientific panels")
gene_panel <- rbind(head(signature[direction == "positive"], 12L), head(signature[direction == "negative"], 12L))
gene_panel[, gene := factor(gene, levels = rev(gene))]
write_tsv(gene_panel, file.path(source_dir, "panel_A_recurrent_genes.tsv"))
pA <- ggplot(gene_panel, aes(median_rho, gene, color = direction)) +
  geom_vline(xintercept = 0, color = "#888888", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = median_rho, yend = gene), linewidth = 0.65) +
  geom_point(size = 2.2) + scale_color_manual(values = direction_colors) +
  labs(title = "A  Recurrent non-circular spatial genes", subtitle = "Same direction in all three Day7 MI slices; block-permutation gate", x = "Median within-slice Spearman rho", y = NULL) + theme_p10
save_panel(pA, "A_recurrent_spatial_genes")

gsea[, axis_direction := fifelse(NES >= 0, "positive", "negative")]
gsea_use <- gsea[is.finite(p.adjust) & is.finite(NES)]
gsea_panel <- gsea_use[order(p.adjust, -abs(NES)), head(.SD, 5L), by = .(collection, axis_direction)]
gsea_panel[, term_label := pretty_term(ID)]
gsea_panel[, term_label := factor(term_label, levels = rev(unique(term_label)))]
write_tsv(gsea_panel, file.path(source_dir, "panel_B_scop_gsea.tsv"))
pB <- ggplot(gsea_panel, aes(NES, term_label, color = axis_direction, size = -log10(p.adjust + 1e-300))) +
  geom_vline(xintercept = 0, color = "#888888", linewidth = 0.3) + geom_point(alpha = 0.9) +
  facet_wrap(~collection, scales = "free_y", ncol = 1) + scale_color_manual(values = direction_colors) +
  scale_size(range = c(1.5, 4.5)) +
  labs(title = "B  Pathways along the remodeling axis", subtitle = "SCOP::RunGSEA; native mouse MSigDB 2025.1", x = "Normalized enrichment score", y = NULL, size = "-log10 BH") + theme_p10 + theme(strip.text = element_text(face = "bold"))
save_panel(pB, "B_scop_gsea", 8.3, 8.5)

source_panel <- cell_scores[timepoint == "Day7" & cells >= 20]
source_panel[, P09_expanded_celltype := factor(P09_expanded_celltype, levels = names(sort(tapply(positive_mean_logcpm, P09_expanded_celltype, mean))))]
write_tsv(source_panel, file.path(source_dir, "panel_C_day7_celltype_expression.tsv"))
pC <- ggplot(source_panel, aes(positive_mean_logcpm, P09_expanded_celltype, color = P09_expanded_celltype)) +
  geom_point(size = 2.2, alpha = 0.9, position = position_jitter(height = 0.08)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.8, color = "black") +
  scale_color_manual(values = celltype_colors) +
  labs(title = "C  Day7 cell-class expression source", subtitle = "sample_id x cell type pseudobulk; black diamond = sample mean", x = "Positive-program mean logCPM", y = NULL) + theme_p10 + theme(legend.position = "none")
save_panel(pC, "C_day7_celltype_expression")

effect_panel <- cell_effects[comparison == "Day7_vs_Day3" & n_target >= 2 & n_reference >= 2]
effect_panel[, source_class := factor(source_class, levels = source_class[order(mean_difference)])]
write_tsv(effect_panel, file.path(source_dir, "panel_D_celltype_effects.tsv"))
pD <- ggplot(effect_panel, aes(mean_difference, source_class, color = source_class)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) + geom_point(size = 2.8) +
  scale_color_manual(values = celltype_colors) +
  labs(title = "D  Cell-class temporal change", subtitle = "Day7 minus Day3; sample-level pseudobulk, descriptive exact tests", x = "Mean signature-score difference", y = NULL) + theme_p10 + theme(legend.position = "none")
save_panel(pD, "D_celltype_temporal_effects")

macro_panel <- macro_effects[comparison == "Day7_vs_Day3" & n_target >= 2 & n_reference >= 2]
macro_panel[, source_class := factor(source_class, levels = source_class[order(mean_difference)])]
write_tsv(macro_panel, file.path(source_dir, "panel_E_macrophage_state_effects.tsv"))
pE <- ggplot(macro_panel, aes(mean_difference, source_class)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) + geom_point(size = 3, color = "#7B3294") +
  labs(title = "E  Refined macrophage-state change", subtitle = "Day7 minus Day3; only states with >=2 sample pseudobulks/group", x = "Mean signature-score difference", y = NULL) + theme_p10
save_panel(pE, "E_macrophage_state_effects")

decomp_features <- c(
  "P09E_RCTD_prop_Fibroblast", "P09E_RCTD_prop_Macrophage", "P09M_RCTD_prop_Trem2_Spp1_macrophage",
  "P09M_RCTD_prop_MHCII_macrophage", "P09M_RCTD_prop_Resident_Folr2_macrophage",
  "P09E_RCTD_prop_Cardiomyocyte", "P09E_RCTD_prop_Endothelial", "P09N_inflammatory_axis",
  "P09N_vascular_axis", "P09N_remote_axis"
)
quartile_panel <- quartile[feature %in% decomp_features]
quartile_panel[, feature_label := gsub("^(P09E_RCTD_prop_|P09M_RCTD_prop_|P09N_)", "", feature)]
quartile_panel[, feature_label := gsub("_", " ", feature_label, fixed = TRUE)]
quartile_panel[, feature_label := factor(feature_label, levels = feature_label[order(median_q4_minus_q1)])]
write_tsv(quartile_panel, file.path(source_dir, "panel_F_axis_quartile_decomposition.tsv"))
pF <- ggplot(quartile_panel, aes(median_q4_minus_q1, feature_label, color = median_q4_minus_q1 > 0)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) + geom_point(size = 3) +
  scale_color_manual(values = c(`FALSE` = "#2F6FB0", `TRUE` = "#C81D25"), guide = "none") +
  labs(title = "F  Biological decomposition of axis extremes", subtitle = "Median Q4-Q1 across three Day7 slices", x = "Median Q4 minus Q1", y = NULL) + theme_p10
save_panel(pF, "F_axis_quartile_decomposition")

external_ann_panel <- head(external_annotation[external_group == "LAD_Day7" & cells >= 20][order(-mean_score)], 14L)
external_ann_panel[, annotation := factor(annotation, levels = rev(annotation))]
write_tsv(external_ann_panel, file.path(source_dir, "panel_G_external_annotations.tsv"))
pG <- ggplot(external_ann_panel, aes(mean_score, annotation, size = cells)) +
  geom_point(color = "#C81D25", alpha = 0.85) +
  labs(title = "G  Independent LAD Day7 annotation projection", subtitle = "GSE280373; cell-level descriptive summary, no animal-level inference", x = "Mean projected score", y = NULL, size = "Cells") + theme_p10
save_panel(pG, "G_external_annotations")

external_groups_use <- c("Sham_Day3", "LAD_Day3", "LAD_Day7", "Cryo_Day3", "Cryo_Day7")
external_lineage_panel <- external_lineage[external_group %in% external_groups_use]
external_lineage_panel[, external_group := factor(external_group, levels = external_groups_use)]
write_tsv(external_lineage_panel, file.path(source_dir, "panel_H_external_lineages.tsv"))
pH <- ggplot(external_lineage_panel, aes(external_group, mean_score, color = lineage, group = lineage)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.3) + geom_line(linewidth = 0.7) + geom_point(size = 2.3) +
  scale_color_manual(values = c(Fibroblast = "#FDAE61", Macrophage = "#7B3294", Other = "#777777")) +
  labs(title = "H  External lineage projection", subtitle = "Descriptive only; groups are not biological replicates", x = NULL, y = "Mean projected score", color = "Lineage") + theme_p10 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(pH, "H_external_lineages")

quality_panel <- quality[condition == "MI" & timepoint == "Day7"]
write_tsv(quality_panel, file.path(source_dir, "panel_I_quality_covariates.tsv"))
pI <- ggplot(quality_panel, aes(covariate, spearman_rho, color = sample_key, group = sample_key)) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.3) + geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  labs(title = "I  Technical-covariate challenge", subtitle = "Within-slice Day7 Spearman correlations; confounding not assumed absent", x = NULL, y = "rho", color = "Slice") + theme_p10 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(pI, "I_quality_covariates")

moran_panel <- moran[condition == "MI"]
write_tsv(moran_panel, file.path(source_dir, "panel_J_spatial_autocorrelation.tsv"))
pJ <- ggplot(moran_panel, aes(sample_key, moran_I, color = feature, shape = feature)) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.3) + geom_point(size = 2.6) +
  labs(title = "J  Spatial autocorrelation", subtitle = "Within-slice 6-nearest-neighbor Moran I; 499 permutations", x = NULL, y = "Moran I", color = NULL, shape = NULL) +
  theme_p10 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(pJ, "J_spatial_autocorrelation")

external_scatter <- external_direction
label_genes <- unique(c(head(signature[direction == "positive", gene], 5L), head(signature[direction == "negative", gene], 5L)))
external_scatter[, label := fifelse(gene %in% label_genes, gene, "")]
write_tsv(external_scatter, file.path(source_dir, "panel_K_external_gene_direction.tsv"))
pK <- ggplot(external_scatter, aes(internal_median_rho, external_LAD_D7_vs_D3, color = signature_direction)) +
  geom_hline(yintercept = 0, color = "#999999", linewidth = 0.3) + geom_vline(xintercept = 0, color = "#999999", linewidth = 0.3) +
  geom_point(alpha = 0.65, size = 1.7) + geom_text(data = external_scatter[label != ""], aes(label = label), size = 2.5, check_overlap = TRUE, vjust = -0.5) +
  scale_color_manual(values = direction_colors) +
  labs(title = "K  External temporal direction challenge", subtitle = "Gene-level descriptive concordance; negative rho is retained", x = "Internal median spatial rho", y = "External LAD Day7 vs Day3 log2FC", color = "Signature") + theme_p10
save_panel(pK, "K_external_gene_direction")

slice_panel <- spot_scores[center == TRUE & is.finite(signature_score), .(
  centers = .N, mean_signature_score = mean(signature_score), median_signature_score = median(signature_score)
), by = .(sample_key, condition, timepoint)]
slice_panel[, group := factor(paste(condition, timepoint), levels = c("Sham Day3", "Needle-pass Day3", "MI 1HR", "MI 4HR", "MI Day3", "MI Day7"))]
write_tsv(slice_panel, file.path(source_dir, "panel_L_slice_signature_summary.tsv"))
pL <- ggplot(slice_panel, aes(group, mean_signature_score, color = timepoint)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.3) + geom_point(size = 3, position = position_jitter(width = 0.08)) +
  scale_color_manual(values = time_colors) +
  labs(title = "L  Slice-level signature summary", subtitle = "Each point is one physical slice; spots are not replicates", x = NULL, y = "Mean center-spot score", color = "Time") +
  theme_p10 + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(pL, "L_slice_signature_summary")

stamp("Creating native SCOP GroupHeatmap")
pb_object <- readRDS(require_file(file.path(dirs$objects, "phase_10_reference_pseudobulk_seurat.rds")))
pb_object <- NormalizeData(pb_object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
heatmap_features <- c(
  head(signature[direction == "positive" & gene %in% rownames(pb_object), gene], 15L),
  head(signature[direction == "negative" & gene %in% rownames(pb_object), gene], 15L)
)
if (length(heatmap_features) < 20L) stop("Too few signature genes for SCOP GroupHeatmap")
pb_object$timepoint <- factor(as.character(pb_object$timepoint), levels = c("Day0", "1HR", "4HR", "Day1", "Day3", "Day7"))
heatmap_result <- GroupHeatmap(
  pb_object, features = heatmap_features,
  group.by = c("P09_expanded_celltype", "timepoint"),
  layer = "data", assay = "RNA", exp_method = "zscore", lib_normalize = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = FALSE,
  heatmap_palcolor = c("#2F6FB0", "#F7F7F7", "#C81D25"),
  group_palcolor = list(celltype_colors, time_colors),
  use_raster = FALSE, verbose = TRUE
)
heatmap_png <- file.path(final_dir, "phase_10_scop_groupheatmap.png")
heatmap_pdf <- file.path(final_dir, "phase_10_scop_groupheatmap.pdf")
png(heatmap_png, width = 5400, height = 3300, res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(heatmap_result$g_tree)
dev.off()
cairo_pdf(heatmap_pdf, width = 18, height = 11)
grid::grid.newpage()
grid::grid.draw(heatmap_result$g_tree)
dev.off()

stamp("Assembling multiscale scientific-structure preview")
overview <- (pA | pB | pC) / (pD | pE | pF) / (pG | pH | pI) / (pJ | pK | pL) +
  plot_annotation(
    title = "Phase 10: multi-scale decomposition of the Day7 spatial remodeling continuum",
    subtitle = "Measured/computed panels only; causal direction is not inferred",
    theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11))
  )
overview_png <- file.path(final_dir, "phase_10_multiscale_summary.png")
overview_pdf <- file.path(final_dir, "phase_10_multiscale_summary.pdf")
ggsave(overview_png, overview, width = 23, height = 27, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(overview_pdf, overview, width = 23, height = 27, device = cairo_pdf, limitsize = FALSE)

stamp("Creating SCOP SpatialSpotPlot contact sheet")
contact_sheet <- file.path(final_dir, "phase_10_spatial_signature_contact_sheet.jpg")
montage_status <- system2(
  "/usr/bin/montage",
  args = c(map_manifest$map, "-tile", "4x3", "-geometry", "900x800+12+12", "-background", "white", contact_sheet),
  stdout = TRUE, stderr = TRUE
)
if (!file.exists(contact_sheet) || file.info(contact_sheet)$size < 100000) {
  stop("Spatial map contact sheet was not created: ", paste(montage_status, collapse = " | "))
}

figure_manifest <- rbindlist(list(
  data.table(panel = LETTERS[1:12], file = file.path(panel_dir, paste0(c(
    "A_recurrent_spatial_genes", "B_scop_gsea", "C_day7_celltype_expression", "D_celltype_temporal_effects",
    "E_macrophage_state_effects", "F_axis_quartile_decomposition", "G_external_annotations", "H_external_lineages",
    "I_quality_covariates", "J_spatial_autocorrelation", "K_external_gene_direction", "L_slice_signature_summary"
  ), ".png")), method = "table-derived ggplot2", inference_unit = c(
    "gene", "gene set", "sample_id", "sample_id", "sample_id", "physical slice",
    "external descriptive cells", "external descriptive cells", "within-slice correlation", "physical slice",
    "gene-level descriptive", "physical slice"
  )),
  data.table(panel = c("SCOP_heatmap", "SCOP_maps", "overview"),
             file = c(heatmap_png, contact_sheet, overview_png),
             method = c("scop::GroupHeatmap", "scop::SpatialSpotPlot", "patchwork scientific-structure preview"),
             inference_unit = c("sample_id pseudobulk", "physical slice", "mixed; see panel manifest"))
), fill = TRUE)
figure_manifest[, bytes := file.info(file)$size]
write_tsv(figure_manifest, file.path(dirs$metadata, "phase_10_figure_manifest.tsv"))

validation <- data.table(
  check = c("table_panels", "panel_pngs", "overview_png", "overview_pdf", "scop_groupheatmap_png", "scop_groupheatmap_pdf", "scop_spatial_maps", "contact_sheet", "source_data_tables"),
  value = c(12L, sum(file.exists(figure_manifest[panel %in% LETTERS[1:12], file])), file.info(overview_png)$size, file.info(overview_pdf)$size,
            file.info(heatmap_png)$size, file.info(heatmap_pdf)$size, nrow(map_manifest), file.info(contact_sheet)$size,
            length(list.files(source_dir, pattern = "\\.tsv$"))),
  status = c("PASS", ifelse(sum(file.exists(figure_manifest[panel %in% LETTERS[1:12], file])) == 12L, "PASS", "FAIL"),
             ifelse(file.info(overview_png)$size > 500000, "PASS", "FAIL"), ifelse(file.info(overview_pdf)$size > 50000, "PASS", "FAIL"),
             ifelse(file.info(heatmap_png)$size > 200000, "PASS", "FAIL"), ifelse(file.info(heatmap_pdf)$size > 10000, "PASS", "FAIL"),
             ifelse(nrow(map_manifest) == 12L && all(map_manifest$bytes > 10000), "PASS", "FAIL"),
             ifelse(file.info(contact_sheet)$size > 100000, "PASS", "FAIL"),
             ifelse(length(list.files(source_dir, pattern = "\\.tsv$")) == 12L, "PASS", "FAIL"))
)
write_tsv(validation, file.path(dirs$metadata, "phase_10_visualization_validation.tsv"))
if (!all(validation$status == "PASS")) stop("Phase 10 visualization validation failed")
stamp("Visualization module PASS")
phase10_log_stop("phase_10_visualize")
