#!/usr/bin/env Rscript
source("/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_redesign/script/figure_common_v2.R")

STAGE_V2 <- file.path(V2,"tmp","figure5_6_refined")
ASSET_V2 <- file.path(STAGE_V2,"assets")
save_asset_original <- save_spatial_asset_v2
save_spatial_asset_v2 <- function(plot,id,width_mm,height_mm) {
 p<-file.path(ASSET_V2,paste0(id,'.png'))
 if(file.exists(p))return(p)
 save_asset_original(plot,id,width_mm,height_mm)
}
dir.create(ASSET_V2, recursive = TRUE, showWarnings = FALSE)

P05 <- file.path(ROOT, "21_scop_expansion", "P05_SECRETED_AXIS_GATE_AMENDED", "tables")
X01 <- file.path(ROOT, "21_scop_expansion", "X01_OMICVERSE_CELLORACLE", "tables")
MER <- file.path(V1, "source_data", "meringue")
d7 <- c("MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3")
targets <- c("Bhlhe40", "Runx2", "Arid5b", "Mef2c", "Maf", "Spi1")
contractile_genes <- c("Fabp3", "Myl3", "Myl2", "Myh6")
fibrotic_genes <- c("Serpinh1", "Thbs2", "Col5a1", "Cthrc1")

message("Figure 5 V2: reading frozen inputs")
spot <- read_tsv_v2(file.path(ROOT, "metadata", "phase_10_spot_signature_scores.tsv.gz"))
signature <- read_tsv_v2(file.path(ROOT, "metadata", "phase_10_day7_remodeling_signature.tsv"))
gsea <- read_tsv_v2(file.path(ROOT, "metadata", "phase_10_scop_gsea.tsv.gz"))
cellscore <- read_tsv_v2(file.path(ROOT, "metadata", "phase_10_reference_celltype_signature_scores.tsv"))
decomp <- read_tsv_v2(file.path(ROOT, "metadata", "phase_10_day7_axis_decomposition_summary.tsv"))
axes <- read_tsv_v2(file.path(P05, "frozen_secreted_axis.tsv"))[order(selection_rank)]
axis_matrix <- read_tsv_v2(file.path(P05, "secreted_axis_candidate_matrix.tsv"))
sample_effect <- read_tsv_v2(file.path(X01, "sample_level_effects.tsv"))
run_summary <- read_tsv_v2(file.path(X01, "run_summary.tsv"))

if (!identical(axes$selection_rank, 1:10)) stop("Frozen secreted-axis rank is not exactly 1:10")
if (!identical(targets, run_summary$target[match(targets, run_summary$target)])) stop("CellOracle target set mismatch")
if (any(run_summary[target %in% targets, status] != "NEGATIVE_COMPLETE")) stop("CellOracle status drift detected")

auto <- rbindlist(lapply(d7, function(s) {
  f <- file.path(MER, paste0(sub("Day", "D", s), "_autocorrelation.tsv"))
  read_tsv_v2(f)
}), fill = TRUE)
auto[, sample_key := sub("^MI_D7_", "MI_Day7_", biological_sample_id)]
fu <- read_tsv_v2(file.path(ROOT, "21_scop_expansion", "P03_MERINGUE", "tables", "feature_universe.tsv"))
tech_genes <- unique(fu[mitochondrial == TRUE | ribosomal == TRUE | hemoglobin == TRUE, feature])
auto[, technical := feature %in% tech_genes |
       grepl("^(mt-|Rpl|Rps|Hb[ab]|Gm42418$|Malat1$)", feature, ignore.case = FALSE)]

object_path <- function(sample_key) {
  file.path(ROOT, "objects", "phase_09", "niches", paste0(sample_key, "_phase09_niche.rds"))
}

message("Figure 5 V2: rendering spatial source panels")
day7_scores <- spot[sample_key %in% d7 & is.finite(signature_score), signature_score]
day7_lim <- max(abs(as.numeric(quantile(day7_scores, c(0.01, 0.99), na.rm = TRUE))))
if (!is.finite(day7_lim) || day7_lim <= 0) stop("Invalid shared Day 7 continuum limit")
continuum_assets <- lapply(seq_along(d7), function(i) {
  sk <- d7[i]
  vals <- setNames(spot[sample_key == sk, signature_score], spot[sample_key == sk, barcode])
  spatial_continuous_asset_v2(
    object_path(sk), vals, paste0("continuum_", sk),
    limits = c(-day7_lim, day7_lim),
    palette = c(COL_CONTRACTILE, COL_ZERO, COL_FIBROTIC), midpoint = 0,
    width_mm = 58, height_mm = 43, point_size = 0.54
  )
})

get_gene_values <- function(sample_key, gene) {
  obj <- readRDS(object_path(sample_key))
  assay <- DefaultAssay(obj)
  mat <- suppressWarnings(tryCatch(LayerData(obj, assay = assay, layer = "data"), error = function(e) NULL))
  raw_layer <- is.null(mat) || !gene %in% rownames(mat)
  if (raw_layer) mat <- LayerData(obj, assay = assay, layer = "counts")
  if (!gene %in% rownames(mat)) stop("Gene missing: ", gene, " in ", sample_key)
  vals <- as.numeric(mat[gene, colnames(obj)])
  if (raw_layer) vals <- log1p(vals)
  names(vals) <- colnames(obj)
  rm(obj, mat)
  gc(verbose = FALSE)
  vals
}

gene_limits <- list()
gene_assets <- list()
for (gene in c("Sdf4", "Abca1")) {
  vals_by_sample <- lapply(d7, get_gene_values, gene = gene)
  lim <- as.numeric(quantile(unlist(vals_by_sample, use.names = FALSE), c(0.01, 0.99), na.rm = TRUE))
  if (!all(is.finite(lim)) || diff(lim) <= 0) stop("Invalid gene-specific spatial limit for ", gene)
  gene_limits[[gene]] <- lim
  for (i in seq_along(d7)) {
    key <- paste(gene, d7[i], sep = "|")
    gene_assets[[key]] <- spatial_continuous_asset_v2(
      object_path(d7[i]), vals_by_sample[[i]], paste0("meringue_", gene, "_", d7[i]),
      limits = lim, palette = colorRampPalette(c("#EFF5FA","#74B5C1","#234F70"))(9),
      width_mm = 51, height_mm = 35, point_size = 0.50
    )
  }
}

message("Figure 5 V2: assembling a-h editorial blocks")
hero_plate <- spatial_plate_v2(
  continuum_assets, nrow = 1, ncol = 3,
  col_labels = unname(sample_labels_v2[d7]), header_height = 0.09
)
hero_legend <- plot_spacer() /
  continuous_bar_v2(c(-day7_lim, day7_lim), c(COL_CONTRACTILE, COL_ZERO, COL_FIBROTIC),
                    "Continuum score", midpoint = 0, vertical = TRUE) +
  plot_layout(heights = c(0.09, 1))
panel_a <- evidence_header_v2(
  (hero_plate | hero_legend) + plot_layout(widths = c(1, 0.045)),
  "a", "Day 7 fibrotic-to-contractile continuum"
)

role_group_v2 <- function(key) {
  fifelse(grepl("^Sham", key), "Sham",
    fifelse(grepl("^MI_Day3", key), "MI Day 3",
      fifelse(grepl("^MI_Day7", key), "MI Day 7",
        fifelse(grepl("^MI_(1HR|4HR)", key), "Descriptive", "Needle-pass"))))
}
physical <- spot[center == TRUE & is.finite(signature_score),
  .(median_score = median(signature_score), n_center_spots = .N),
  by = .(sample_key, condition, timepoint)]
if (!setequal(physical$sample_key, sample_order_v2) || nrow(physical) != 12L) stop("Physical-tissue summary is not 12 unique tissues")
physical[, sample_key := factor(sample_key, levels = sample_order_v2)]
physical[, role := role_group_v2(as.character(sample_key))]
score_span <- diff(range(physical$median_score))
strip_y <- max(physical$median_score) + 0.11 * score_span
physical_labels <- sample_labels_v2
physical_labels[c("MI_1HR_rep1", "MI_4HR_rep1")] <- paste0(physical_labels[c("MI_1HR_rep1", "MI_4HR_rep1")], "\ndescriptive")
physical_plot <- ggplot(physical, aes(sample_key, median_score)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "#7B8188") +
  geom_tile(aes(y = strip_y, fill = role), width = 0.86, height = 0.035 * score_span) +
  geom_point(aes(colour = role), size = 2.0) +
  scale_fill_manual(values = role_cols_v2, guide = "none") +
  scale_colour_manual(values = role_cols_v2, guide = "none") +
  scale_x_discrete(labels = physical_labels[sample_order_v2]) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = NULL, y = "Median center-spot score") +
  theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, lineheight = 0.88),
        plot.margin = margin(0.5, 2, 1, 2, unit = "mm"))
panel_b <- evidence_header_v2(physical_plot, "b", "Continuum by physical tissue",
                              "One summary point per physical tissue; 1 h and 4 h are descriptive")

gene_order <- c(contractile_genes, fibrotic_genes)
gene_panel <- signature[gene %in% gene_order, .(gene, median_rho)]
if (nrow(gene_panel) != 8L || uniqueN(gene_panel$gene) != 8L) stop("Frozen eight-gene panel mismatch")
gene_panel[, pole := ifelse(gene %in% contractile_genes, "Contractile", "Fibrotic")]
gene_panel[, gene := factor(gene, levels = rev(gene_order))]
gene_plot <- ggplot(gene_panel, aes(median_rho, gene, colour = pole)) +
  geom_vline(xintercept = 0, linewidth = 0.30, colour = "#7B8188") +
  geom_segment(aes(x = 0, xend = median_rho, yend = gene), linewidth = 0.55) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c(Contractile = COL_CONTRACTILE, Fibrotic = COL_FIBROTIC), guide = "none") +
  scale_x_continuous(limits = c(-0.56, 0.56), breaks = c(-0.5, 0, 0.5)) +
  labs(x = "Median within-slice ρ", y = NULL) +
  theme_editorial_v2() +
  theme(axis.text.y = element_text(face = "italic"), plot.margin = margin(1, 2, 1, 2, unit = "mm"))
panel_c <- evidence_header_v2(gene_plot, "c", "Genes outside axis construction", "No P values")

path_spec <- data.table(
  collection = rep(c("GO_BP", "Hallmark", "Reactome"), each = 3),
  ID = c(
    "GOBP_AEROBIC_RESPIRATION",
    "GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",
    "GOBP_CELLULAR_RESPIRATION",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ADIPOGENESIS",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_COMPLEX_I_BIOGENESIS"
  ),
  term = c(
    "Aerobic Respiration",
    "Generation of Precursor Metabolites and Energy",
    "Cellular Respiration",
    "Oxidative Phosphorylation",
    "Epithelial Mesenchymal Transition",
    "Adipogenesis",
    "Aerobic Respiration and Respiratory Electron Transport",
    "Respiratory Electron Transport",
    "Complex I Biogenesis"
  ),
  term_label = c(
    "Aerobic\nRespiration",
    "Generation of Precursor\nMetabolites and Energy",
    "Cellular\nRespiration",
    "Oxidative\nPhosphorylation",
    "Epithelial Mesenchymal\nTransition",
    "Adipogenesis",
    "Aerobic Respiration\nand Respiratory\nElectron Transport",
    "Respiratory\nElectron Transport",
    "Complex I Biogenesis"
  )
)
path_spec[, display_order := .I]
path_df <- merge(path_spec, gsea[, .(collection, ID, NES)], by = c("collection", "ID"), all.x = TRUE, sort = FALSE)
if (nrow(path_df) != 9L || any(!is.finite(path_df$NES))) stop("Frozen nine-pathway panel mismatch")
setorder(path_df, display_order)
path_df[, collection_label := factor(collection, levels = c("GO_BP", "Hallmark", "Reactome"),
                                     labels = c("GO BP", "Hallmark", "Reactome"))]
path_df[, col_id := match(collection, c("GO_BP", "Hallmark", "Reactome")) - 1L]
path_df[, row_id := seq_len(.N), by = collection]
path_df[, y := c(0.70, 0.42, 0.14)[row_id]]
path_df[, text_x := col_id + 0.02]
path_df[, zero_x := col_id + 0.865]
path_df[, end_x := col_id + 0.75 + (NES + 4.15) / (3.15 + 4.15) * 0.23]
collection_headers <- data.table(
  x = c(0.5, 1.5, 2.5), y = 0.95,
  label = c("GO BP", "Hallmark", "Reactome")
)
zero_lines <- data.table(x = c(0.865, 1.865, 2.865), xend = c(0.865, 1.865, 2.865),
                         y = 0.06, yend = 0.79)
path_plot <- ggplot(path_df) +
  geom_vline(xintercept = c(1, 2), linewidth = 0.18, colour = "#E2E5E8") +
  geom_segment(data = zero_lines, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.28, colour = "#7B8188") +
  geom_text(aes(x = text_x, y = y, label = term_label), hjust = 0, lineheight = 0.82,
            family = FONT_V2, size = 6.5 / GG_PT, colour = COL_AXIS) +
  geom_segment(aes(x = zero_x, xend = end_x, y = y, yend = y, colour = NES), linewidth = 0.50) +
  geom_point(aes(x = end_x, y = y, colour = NES), size = 1.9) +
  geom_text(data = collection_headers, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, family = FONT_V2, fontface = "bold",
            size = 7.2 / GG_PT, colour = COL_AXIS) +
  scale_colour_gradient2(low = COL_CONTRACTILE, mid = COL_ZERO, high = COL_FIBROTIC,
                         midpoint = 0, limits = c(-4.15, 4.15), guide = "none") +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 1), expand = FALSE, clip = "on") +
  theme_void(base_family = FONT_V2) +
  theme(plot.margin = margin(0, 1, 0, 1, unit = "mm"))
panel_d <- evidence_header_v2(
  path_plot, "d", "Pathway polarity",
  "Shared NES scale: -4 to +3 | nine frozen pathways | no significance stars"
)

cell_order <- c("Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage", "Monocyte",
                "Neutrophil", "Dendritic", "B_cell", "T_NK", "Pericyte", "Schwann_cell")
time_order <- c("Day0", "1HR", "4HR", "Day1", "Day3", "Day7")
cell_summary <- cellscore[, .(
  median_signature = median(remodeling_signature_score, na.rm = TRUE),
  n_samples = uniqueN(sample_id)
), by = .(P09_expanded_celltype, timepoint)]
cell_summary <- cell_summary[P09_expanded_celltype %in% cell_order & timepoint %in% time_order]
cell_summary[, P09_expanded_celltype := factor(P09_expanded_celltype, levels = rev(cell_order))]
cell_summary[, timepoint := factor(timepoint, levels = time_order,
                                   labels = c("Day 0", "1 h", "4 h", "Day 1", "Day 3", "Day 7"))]
cell_lim <- max(abs(as.numeric(quantile(cell_summary$median_signature, c(0.02, 0.98), na.rm = TRUE))))
cell_plot <- ggplot(cell_summary, aes(timepoint, P09_expanded_celltype)) +
  geom_point(aes(size = n_samples, fill = median_signature), shape = 21, stroke = 0.25, colour = "white") +
  geom_point(data = cell_summary[n_samples == 1L], aes(size = n_samples), shape = 1,
             stroke = 0.55, colour = COL_AXIS, show.legend = FALSE) +
  scale_fill_gradient2(low = COL_CONTRACTILE, mid = COL_ZERO, high = COL_FIBROTIC,
                       midpoint = 0, limits = c(-cell_lim, cell_lim), oob = scales::squish,
                       name = "Median signature") +
  scale_size_continuous(range = c(1.6, 3.4), breaks = sort(unique(cell_summary$n_samples)),
                        name = "Samples") +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colorbar(order = 1, barheight = unit(14, "mm"), barwidth = unit(2.2, "mm")),
         size = "none") +
  theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right",
        plot.margin = margin(0.5, 1, 0.5, 1, unit = "mm"))

decomp_top <- copy(decomp[order(-abs(median_q4_minus_q1))][1:4])
decomp_top[, label := gsub("^(P09E_RCTD_prop_|P09M_RCTD_prop_|P09N_)", "", feature)]
decomp_top[, label := gsub("_", " ", label)]
decomp_top[, label := factor(label, levels = label[order(median_q4_minus_q1)])]
decomp_plot <- ggplot(decomp_top, aes(median_q4_minus_q1, label, colour = median_q4_minus_q1)) +
  geom_vline(xintercept = 0, linewidth = 0.30, colour = "#7B8188") +
  geom_segment(aes(x = 0, xend = median_q4_minus_q1, yend = label), linewidth = 0.50) +
  geom_point(size = 2.0) +
  scale_colour_gradient2(low = COL_CONTRACTILE, mid = COL_ZERO, high = COL_FIBROTIC,
                         midpoint = 0, guide = "none") +
  labs(x = "Fibrotic Q4 - contractile Q1", y = NULL) +
  theme_editorial_v2() +
  theme(axis.text.y = element_text(size = 6.5), plot.margin = margin(0.5, 1, 0.5, 1, unit = "mm"))
panel_e <- evidence_header_v2(
  (cell_plot | decomp_plot) + plot_layout(widths = c(0.62, 0.38)),
  "e", "Cell-source localization and Q4-Q1 decomposition",
  "Color: sample median; area: sample count; open ring: n = 1"
)

d7_top <- auto[sample_key %in% d7 & !technical & rank <= 100]
rep_genes <- d7_top[, .N, by = feature][N == 3, feature]
ranker <- d7_top[feature %in% rep_genes,
  .(median_rank = median(rank), worst_rank = max(rank), median_stat = median(statistic)), by = feature]
setorder(ranker, median_rank, worst_rank, -median_stat, feature)
rep_genes <- head(ranker$feature, 16)
if (length(rep_genes) != 16L) stop("MERINGUE 3/3 top-100 frozen panel does not contain 16 genes")
mer_matrix <- auto[feature %in% rep_genes & sample_key %in% d7]
mer_matrix[, feature := factor(feature, levels = rep_genes)]
mer_matrix[, sample_key := factor(sample_key, levels = rev(d7), labels = rev(unname(sample_labels_v2[d7])))]
mer_plot <- ggplot(mer_matrix, aes(feature, sample_key, fill = statistic)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  geom_point(aes(size = pmin(rank, 100)), shape = 21, fill = NA, stroke = 0.28, colour = COL_AXIS) +
  scale_size_continuous(range = c(2.1, 0.55), limits = c(1, 100), guide = "none") +
  scale_fill_gradientn(colours=c("#EFF5FA","#74B5C1","#234F70"), name="Statistic") +
  labs(x = NULL, y = NULL) +
  theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, face = "italic", size = 6.5),
        legend.position = "right", legend.key.height = unit(5, "mm"),
        plot.margin = margin(0, 1, 0, 1, unit = "mm"))

atlas_header <- wrap_plots(
  c(list(plot_spacer()), lapply(unname(sample_labels_v2[d7]), strip_text_v2)),
  ncol = 4, widths = c(0.075, 1, 1, 1)
)
atlas_sdf4 <- wrap_plots(
  c(list(strip_text_v2("Sdf4", hjust = 1)), lapply(unname(gene_assets[paste("Sdf4", d7, sep = "|")]), img_panel_compact_v2)),
  ncol = 4, widths = c(0.075, 1, 1, 1)
)
atlas_abca1 <- wrap_plots(
  c(list(strip_text_v2("Abca1", hjust = 1)), lapply(unname(gene_assets[paste("Abca1", d7, sep = "|")]), img_panel_compact_v2)),
  ncol = 4, widths = c(0.075, 1, 1, 1)
)
atlas_core <- wrap_plots(list(atlas_header, atlas_sdf4, atlas_abca1), ncol = 1, heights = c(0.09, 1, 1))
atlas_legends <- wrap_plots(list(
  strip_text_v2("log1p expr.", size = 6.5),
  continuous_bar_v2(gene_limits$Sdf4, colorRampPalette(c("#EFF5FA","#74B5C1","#234F70"))(9), "") +
    theme(axis.title.y.right = element_blank()),
  continuous_bar_v2(gene_limits$Abca1, colorRampPalette(c("#EFF5FA","#74B5C1","#234F70"))(9), "") +
    theme(axis.title.y.right = element_blank())
), ncol = 1, heights = c(0.09, 1, 1))
atlas <- wrap_plots(list(atlas_core, atlas_legends), ncol = 2, widths = c(1, 0.047))
panel_f <- evidence_header_v2(
  mer_plot / atlas + plot_layout(heights = c(0.22, 0.78)),
  "f", "MERINGUE reproducibility", "Descriptive statistic/rank; no P values"
)

axis_context <- axis_matrix[axis_id %in% axes$axis_id]
axis_aux <- axis_context[, .(
  neighborhood = any(neighborhood_support %in% TRUE),
  gradient = if (any(gradient_support %in% TRUE)) "TRUE" else
    if (all(gradient_support == "NOT_EVALUABLE" | is.na(gradient_support))) "NE" else "FALSE"
), by = axis_id]
axis_plot_data <- merge(
  axes[, .(axis_id, selection_rank, ligand_secact_human, receptor,
           n_supported_stages, n_discovery_supporting_tissues,
           minimum_component_coverage, spatialcellchat_median_score)],
  axis_aux, by = "axis_id", all.x = TRUE, sort = FALSE
)[order(selection_rank)]
if (nrow(axis_plot_data) != 10L) stop("Frozen axis display contract failed")
axis_plot_data[, label := paste0(selection_rank, "  ", ligand_secact_human, "→", receptor)]
axis_plot_data[, label := factor(label, levels = rev(label))]
axis_plot_data[, tissue_size := 1.8 + 2.0 * n_discovery_supporting_tissues / 6]
if (any(!is.finite(axis_plot_data$spatialcellchat_median_score))) stop("Non-finite frozen CellChat score")
axis_plot_data[, cellchat_size := if (diff(range(spatialcellchat_median_score)) > 0)
  1.7 + 2.1 * scales::rescale(spatialcellchat_median_score) else 2.7]

axis_half_plot_v2 <- function(ranks) {
  d <- copy(axis_plot_data[selection_rank %in% ranks])
  d[, label := factor(as.character(label), levels = rev(as.character(label)))]
  ggplot(d, aes(y = label)) +
    geom_hline(yintercept = seq_len(nrow(d)), linewidth = 0.18, colour = "#E8EAED") +
    geom_text(aes(x = 1, label = paste0(n_supported_stages, "/2")), family = FONT_V2,
              size = 6.5 / GG_PT, colour = COL_AXIS) +
    geom_point(aes(x = 2, size = tissue_size), shape = 21, stroke = 0.35,
               colour = COL_AXIS, fill = "#9FB9D4") +
    geom_text(aes(x = 2, label = paste0(n_discovery_supporting_tissues, "/6")), family = FONT_V2,
              size = 6.5 / GG_PT, colour = COL_AXIS) +
    geom_tile(aes(x = 3, fill = minimum_component_coverage), width = 0.72, height = 0.72,
              colour = "white", linewidth = 0.2) +
    geom_text(aes(x = 3, label = sprintf("%.0f%%", 100 * minimum_component_coverage)), family = FONT_V2,
              size = 6.5 / GG_PT, colour = COL_AXIS) +
    geom_point(data = d[neighborhood == TRUE], aes(x = 4), shape = 21, size = 2.4,
               stroke = 0.45, colour = COL_AXIS, fill = "#5B8F73") +
    geom_point(data = d[neighborhood != TRUE | is.na(neighborhood)], aes(x = 4), shape = 21,
               size = 2.4, stroke = 0.45, colour = COL_AXIS, fill = "white") +
    geom_point(data = d[gradient == "TRUE"], aes(x = 5), shape = 21, size = 2.4,
               stroke = 0.45, colour = COL_AXIS, fill = "#79549A") +
    geom_point(data = d[gradient == "FALSE"], aes(x = 5), shape = 21, size = 2.4,
               stroke = 0.45, colour = COL_AXIS, fill = "white") +
    geom_text(data = d[gradient == "NE"], aes(x = 5, label = "NE"), family = FONT_V2,
              size = 6.5 / GG_PT, colour = COL_UNCERTAIN) +
    geom_point(aes(x = 6, size = cellchat_size), shape = 21, stroke = 0.35,
               colour = COL_AXIS, fill = "#4F78A7") +
    scale_size_identity() +
    scale_fill_gradient(low = "#EDF3F8", high = "#2E6F9E", limits = c(0, 1), guide = "none") +
    scale_x_continuous(
      breaks = 1:6, labels = c("Stage", "Tissue", "Coverage", "Neighbor", "Gradient", "CellChat"),
      position = "top", limits = c(0.55, 6.45), expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) + theme_editorial_v2() +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 6.5, angle = 42, hjust = 0, vjust = 0.5),
          axis.text.y = element_text(size = 6.5),
          plot.margin = margin(0.5, 1, 0.5, 1, unit = "mm"))
}
axis_plot <- wrap_plots(list(axis_half_plot_v2(1:5), axis_half_plot_v2(6:10)), ncol = 2)
panel_g <- evidence_header_v2(
  axis_plot, "g", "Frozen candidate secreted axes"
)

oracle <- sample_effect[target %in% targets, .(
  median_displacement = median(median_displacement, na.rm = TRUE),
  transition_mass = median(median_transition_mass_to_Macrophage_Spp1, na.rm = TRUE)
), by = target]
oracle <- merge(oracle, run_summary[, .(target, status)], by = "target", all.x = TRUE, sort = FALSE)
oracle[, target := factor(target, levels = rev(targets))]
oracle_long <- melt(oracle, id.vars = c("target", "status"),
                    measure.vars = c("median_displacement", "transition_mass"),
                    variable.name = "metric", value.name = "raw_value")
oracle_long[, metric := factor(metric, levels = c("median_displacement", "transition_mass"),
                               labels = c("Median\ndisplacement", "Transition\nmass"))]
oracle_long[, within_metric_z := if (sd(raw_value, na.rm = TRUE) > 0)
  as.numeric(scale(raw_value)) else 0, by = metric]
status_points <- unique(oracle[, .(target, status)])
metric_levels <- c("Median\ndisplacement", "Transition\nmass", "Status")
oracle_long[, metric := factor(as.character(metric), levels = metric_levels)]
status_points[, metric := factor("Status", levels = metric_levels)]

oracle_plot <- ggplot(oracle_long, aes(metric, target, fill = within_metric_z)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", raw_value)), family = FONT_V2,
            size = 6.5 / GG_PT, colour = COL_AXIS) +
  geom_point(data = status_points, aes(x = metric, y = target), inherit.aes = FALSE,
             shape = 21, size = 2.2, stroke = 0.45, colour = "#68717D", fill = "#EEF0F2") +
  scale_fill_gradient2(low = COL_CONTRACTILE, mid = COL_ZERO, high = COL_FIBROTIC,
                       midpoint = 0, guide = "none") +
  scale_x_discrete(position = "top", limits = metric_levels, expand = expansion(add = c(0.15, 0.15))) +
  labs(x = NULL, y = NULL, caption = "○  NEGATIVE_COMPLETE") +
  coord_cartesian(clip = "off") +
  theme_editorial_v2() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 6.5, lineheight = 0.88),
        axis.text.y = element_text(size = 6.5),
        plot.caption = element_text(size = 6.5, hjust = 0.5, colour = COL_AXIS,
                                    margin = margin(t = 1, unit = "mm")),
        plot.margin = margin(0.5, 1, 0, 1, unit = "mm"))
panel_h <- evidence_header_v2(
  oracle_plot, "h", "CellOracle computational boundary",
  "6/6 NEGATIVE_COMPLETE\nunder the preregistered module-consistency rule", header_mm = 8.0
)


# Separate the expression-axis question from recurrence and computational candidates.
physical_plot <- physical_plot + labs(y="Median signature") +
 scale_x_discrete(labels=setNames(c('Sh1','Sh2','D3r1','D3r2','D3r3','D7r1','D7r2','D7r3','1 h','4 h','N3','N4'),sample_order_v2))
path_plot <- path_plot + geom_text(aes(x=end_x,y=y-.072,label=sprintf('%.2f',NES)),
 family=FONT_V2,size=5.8/GG_PT,colour=COL_AXIS)
draw_block <- function(p,x,y,w,h,letter,title,sub=NULL) {
 pushViewport(viewport(x=unit(x,'mm'),y=unit(225-y,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','top')))
 grid.text(letter,x=unit(0,'mm'),y=unit(1,'npc'),just=c('left','top'),gp=gpar(fontfamily=FONT_V2,fontsize=9,fontface='bold'))
 grid.text(title,x=unit(4,'mm'),y=unit(1,'npc'),just=c('left','top'),gp=gpar(fontfamily=FONT_V2,fontsize=7.4,fontface='bold'))
 hh<-if(is.null(sub))6 else 10
 if(!is.null(sub))grid.text(sub,x=unit(4,'mm'),y=unit(1,'npc')-unit(4,'mm'),just=c('left','top'),gp=gpar(fontfamily=FONT_V2,fontsize=6.2,col='#555D66'))
 pushViewport(viewport(x=.5,y=0,width=1,height=unit(h-hh,'mm'),just=c('center','bottom')))
 grid.draw(if(inherits(p,'patchwork'))patchworkGrob(p)else ggplotGrob(p));upViewport(2)
}
map_direct<-function(path,x,y,w,h){im<-png::readPNG(path);ratio<-dim(im)[2]/dim(im)[1];ww<-min(w,h*ratio);hh<-ww/ratio;grid.raster(im,x=unit(x+w/2,'mm'),y=unit(225-y-h/2,'mm'),width=unit(ww,'mm'),height=unit(hh,'mm'),interpolate=FALSE)}
fig5<-function(){
 grid.newpage()
 draw_block(ggplot()+theme_void(),2,2,179,55,'a','Day 7 expression geography','Shared symmetric score range across the three tissues')
 for(i in 1:3){x<-8+(i-1)*53;grid.text(paste('MI D7 r',i,sep=''),x=unit(x+25,'mm'),y=unit(210,'mm'),gp=gpar(fontfamily=FONT_V2,fontsize=7,fontface='bold'));map_direct(continuum_assets[[i]],x,18,50,35)}
 draw_block(continuous_bar_v2(c(-day7_lim,day7_lim),c(COL_CONTRACTILE,COL_ZERO,COL_FIBROTIC),'Score',midpoint=0),168,13,12,40,'','')
 draw_block(physical_plot,2,59,179,30,'b','Tissue-level signature summaries','One point per physical tissue; 1 h and 4 h are descriptive')
 draw_block(gene_plot,2,94,53,60,'c','Nonconstruction genes','Descriptive correlations')
 draw_block(path_plot,58,94,123,60,'d','Pathway polarity','Numbers show saved normalized enrichment scores')
 draw_block(cell_plot|decomp_plot,2,161,179,61,'e','Cell-source distribution and Q4–Q1 decomposition','Colour: sample median; area: sample count; open ring: n = 1')
}
fig6<-function(){
 grid.newpage()
 draw_block(mer_plot,2,2,179,41,'a','Sixteen recurrent spatial genes','Day 7 top 100 in all three tissues; statistics and ranks are descriptive')
 draw_block(ggplot()+theme_void(),2,46,179,86,'b','Spatial expression in the three Day 7 tissues','Frozen Sdf4 and Abca1 selection; a shared scale within each gene')
 for(j in 1:2)for(i in 1:3){g<-c('Sdf4','Abca1')[j];x<-15+(i-1)*49;y<-60+(j-1)*34
  if(j==1)grid.text(paste('r',i,sep=''),x=unit(x+23,'mm'),y=unit(225-y+2,'mm'),gp=gpar(fontfamily=FONT_V2,fontsize=7,fontface='bold'))
  if(i==1)grid.text(g,x=unit(9,'mm'),y=unit(225-y-16,'mm'),rot=90,gp=gpar(fontfamily=FONT_V2,fontsize=7,fontface='italic'))
  map_direct(gene_assets[[paste(g,d7[i],sep='|')]],x,y,46,32)
 }
 for(j in 1:2){g<-c('Sdf4','Abca1')[j];draw_block(continuous_bar_v2(gene_limits[[g]],colorRampPalette(c('#EFF5FA','#74B5C1','#234F70'))(9),'log1p'),164,56+(j-1)*34,16,37,'','')}
 draw_block(axis_plot,2,137,179,40,'c','Ten frozen ligand–receptor candidates','Original communication labels; no receptor-specific causal claim')
 draw_block(oracle_plot,2,182,179,40,'d','CellOracle simulations under the fixed model','All six targets failed the prespecified module-consistency criterion')
}
for(n in 5:6){stem<-file.path(STAGE_V2,paste0('Figure_',n,'_final'));fun<-get(paste0('fig',n))
 cairo_pdf(paste0(stem,'.pdf'),width=183/25.4,height=225/25.4);fun();dev.off()
 svglite::svglite(paste0(stem,'.svg'),width=183/25.4,height=225/25.4);fun();dev.off()
 rsvg::rsvg_png(paste0(stem,'.svg'),paste0(stem,'.png'),width=2161,height=2657)
}
saveRDS(list(physical=physical,ranker=ranker,axes=axis_plot_data,oracle=oracle),file.path(STAGE_V2,'plotted_values.rds'))
