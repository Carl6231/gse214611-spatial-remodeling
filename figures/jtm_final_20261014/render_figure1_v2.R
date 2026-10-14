#!/usr/bin/env Rscript

source("/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_polished/script/figure_common_v2.R")

BUILD_ROOT <- Sys.getenv("REDESIGN_V2_BUILD_ROOT", unset = "")
if (!nzchar(BUILD_ROOT) || !dir.exists(BUILD_ROOT)) {
  stop("REDESIGN_V2_BUILD_ROOT must name an existing empty build root")
}
STAGE_V2 <- file.path(BUILD_ROOT, "Figure_1_refined")
if (file.exists(STAGE_V2)) stop("Figure 1 staging directory already exists: ", STAGE_V2)
ASSET_V2 <- file.path(STAGE_V2, "assets")

if (!exists("SpatialSpotPlot", envir = asNamespace("scop"), mode = "function", inherits = FALSE)) {
  stop("scop::SpatialSpotPlot is unavailable in the live runtime")
}

F1_SRC <- file.path(ROOT, "results", "figures", "Figure_1", "source_data")
role_from_key_v2 <- function(key) {
  out <- rep(NA_character_, length(key))
  out[grepl("^Sham_", key)] <- "Sham"
  out[grepl("^MI_Day3_", key)] <- "MI Day 3"
  out[grepl("^MI_Day7_", key)] <- "MI Day 7"
  out[grepl("^MI_[14]HR_", key)] <- "Descriptive"
  out[grepl("^NeedlePass_", key)] <- "Needle-pass"
  if (anyNA(out)) stop("Unmapped sample role: ", paste(unique(key[is.na(out)]), collapse = ", "))
  out
}

role_display_cols <- c(
  "Formal comparison" = "#333333",
  "Descriptive only" = COL_DESC,
  "Auxiliary" = COL_NEEDLE
)
celltype_cols <- c(
  Cardiomyocyte = "#5B9BD5", Endothelial = "#22A7C7",
  Fibroblast = "#35B89A", Macrophage = "#D97FB2",
  Pericyte = "#B5A2D8", Schwann_cell = "#81CCEA"
)

save_spot_asset <- function(object_path, id, display_values, display_levels, colours,
                            width_mm, height_mm, point_size = 0.45) {
  if (!file.exists(object_path)) stop("Missing spatial object: ", object_path)
  obj <- readRDS(object_path)
  if (length(display_values) == 1L) display_values <- rep(display_values, ncol(obj))
  if (length(display_values) != ncol(obj)) stop("Spatial display vector length mismatch: ", id)
  obj$.v2_display <- factor(display_values, levels = display_levels)
  present <- display_levels[display_levels %in% unique(as.character(obj$.v2_display))]
  p <- suppressMessages(suppressWarnings(scop::SpatialSpotPlot(
    obj,
    group.by = ".v2_display",
    assay = "Spatial",
    layer = "counts",
    image = Seurat::Images(obj)[[1L]],
    overlay_image = TRUE,
    image.alpha = 0.18,
    crop = TRUE,
    show_axes = FALSE,
    pt.size = point_size,
    pt.alpha = 0.92,
    stroke = 0,
    palcolor = unname(colours[present]),
    legend.position = "none",
    verbose = FALSE
  )))
  p <- suppressMessages(p + scale_fill_manual(values = colours,
                                              limits = display_levels, drop = FALSE))
  stopifnot(identical(unname(ggplot_build(p)$plot$scales$get_scales("fill")$map(present)),
                      unname(colours[present])))
  p <- suppressMessages(suppressWarnings(p + coord_fixed())) +
    theme_void(base_family = FONT_V2) +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  out <- save_spatial_asset_v2(p, id, width_mm, height_mm)
  rm(obj, p)
  gc(verbose = FALSE)
  out
}

labeled_asset_grid <- function(paths, labels, ncol) {
  if (length(paths) != length(labels)) stop("Labeled-grid path/label mismatch")
  cells <- Map(function(path, label) {
    wrap_plots(list(strip_text_v2(label, size = 7.2), img_panel_v2(path)), ncol = 1,
               heights = c(0.11, 0.89))
  }, paths, labels)
  wrap_plots(cells, ncol = ncol)
}

registry <- read_tsv_v2(file.path(F1_SRC, "Figure_1a_sample_design.tsv"))
mapping <- read_tsv_v2(file.path(F1_SRC, "Figure_1b_GSM_to_tissue.tsv"))
selection <- read_tsv_v2(file.path(F1_SRC, "Figure_1e_representative_selection.tsv"))
spot_counts <- read_tsv_v2(file.path(F1_SRC, "Figure_1f_spot_accounting.tsv"))
qc <- read_tsv_v2(file.path(F1_SRC, "Figure_1g_spot_QC.tsv.gz"))
reference <- read_tsv_v2(file.path(F1_SRC, "Figure_1h_reference_celltypes.tsv"))
common <- read_tsv_v2(file.path(F1_SRC, "Figure_1i_common_genes.tsv"))
marker <- read_tsv_v2(file.path(F1_SRC, "Figure_1i_celltype_marker_summary.tsv"))
flow <- read_tsv_v2(file.path(F1_SRC, "Figure_1j_spot_flow.tsv.gz"))

if (uniqueN(registry$GSM) != 10L || uniqueN(registry$sample_key) != 12L || nrow(registry) != 12L) {
  stop("Figure 1 capture/tissue contract changed")
}
role_counts <- registry[, .N, by = analysis_role_display]
expected_roles <- c("Formal comparison" = 8L, "Descriptive only" = 2L, "Auxiliary" = 2L)
actual_roles <- setNames(role_counts$N, role_counts$analysis_role_display)
if (!identical(as.integer(actual_roles[names(expected_roles)]), as.integer(expected_roles))) {
  stop("Figure 1 analysis-role counts changed")
}
if (nrow(qc) != 25548L || uniqueN(qc$spot_key) != 25548L ||
    sum(spot_counts$retained_spots) != 25548L || uniqueN(spot_counts$sample_key) != 12L) {
  stop("Figure 1 retained-spot contract changed")
}
flow_counts <- table(flow$stage_niche)
if (nrow(flow) != 25548L || sum(flow$stage_macrophage == "Macrophage-rich") != 5937L ||
    unname(flow_counts[["Assigned niche"]]) != 4143L ||
    unname(flow_counts[["Niche Uncertain"]]) != 1794L) {
  stop("Figure 1 analytical-flow contract changed")
}
if (sum(reference$reference_cells) != 176526L ||
    reference[celltype == "Schwann_cell", reference_cells] != 28L || nrow(reference) != 6L) {
  stop("Figure 1 reference-cell contract changed")
}
if (uniqueN(common$sample_key) != 12L || any(common$exact_common_genes != 22490L) ||
    nrow(marker) != 72L || uniqueN(marker$sample_key) != 12L || uniqueN(marker$celltype) != 6L) {
  stop("Figure 1 reference-compatibility contract changed")
}

# a | Compact study architecture.
group_levels <- c("Day 3 Sham", "1 h MI", "4 h MI", "Day 3 MI", "Day 7 MI", "Day 3 needle-pass")
group_labels <- c("Day 3\nSham", "MI\n1 h", "MI\n4 h", "MI\nDay 3", "MI\nDay 7", "Day 3\nNeedle-pass")
group_cols <- c(COL_SHAM, COL_DESC, COL_DESC, COL_D3, COL_D7, COL_NEEDLE)
group_summary <- unique(registry[, .(sample_key, analysis_group_display)])[,
  .(n_tissues = .N), by = analysis_group_display]
group_summary <- group_summary[match(group_levels, analysis_group_display)]
if (anyNA(group_summary$n_tissues) || !identical(group_summary$n_tissues, c(2L, 1L, 1L, 3L, 3L, 2L))) {
  stop("Figure 1 study-group counts changed")
}
group_summary[, `:=`(
  x = seq_len(.N), label = group_labels, colour = group_cols,
  box_label = paste0(group_labels, "\nn=", n_tissues)
)]

architecture_plot <- ggplot() +
  geom_rect(data = group_summary,
            aes(xmin = x - 0.41, xmax = x + 0.41, ymin = 0.70, ymax = 1.44),
            fill = "#F5F6F7", colour = group_summary$colour, linewidth = 0.48) +
  geom_text(data = group_summary, aes(x, 1.08, label = box_label), family = FONT_V2,
            size = 7.1 / GG_PT, lineheight = 0.72, colour = COL_AXIS) +
  annotate("segment", x = 1.42, xend = 4.30, y = 0.20, yend = 0.20,
           linewidth = 0.35, colour = "#8A9199", arrow = arrow(length = unit(1.7, "mm"))) +
  annotate("text", x = 0.55, y = 0.28, label = "10 GEO\ncaptures", hjust = 0,
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 3.0, y = 0.47,
           label = "8 one-to-one · GSM6613080 / GSM6613081 branch only",
           family = FONT_V2, size = 7.1 / GG_PT, colour = "#555D66") +
  annotate("text", x = 5.95, y = 0.28,
           label = "12 physical tissues\n8 formal · 2 descriptive · 2 auxiliary", hjust = 1,
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold", lineheight = 0.9) +
  coord_cartesian(xlim = c(0.45, 6.15), ylim = c(0.12, 1.47), expand = FALSE, clip = "off") +
  theme_void(base_family = FONT_V2) + theme(plot.margin = margin(0, 1, 0, 1, unit = "mm"))
panel_a <- evidence_header_v2(architecture_plot, "a", "Study design")

# b | Two-by-two mixed-capture hero plate.
mixed_raw <- c(
  GSM6613080 = file.path(ROOT, "objects", "raw", "GSM6613080_raw.rds"),
  GSM6613081 = file.path(ROOT, "objects", "raw", "GSM6613081_raw.rds")
)
mixed_labeled <- c(
  GSM6613080 = file.path(ROOT, "objects", "mixed_labeled", "GSM6613080_mixed_labeled.rds"),
  GSM6613081 = file.path(ROOT, "objects", "mixed_labeled", "GSM6613081_mixed_labeled.rds")
)
before_assets <- Map(function(path, gsm) {
  obj <- readRDS(path)
  values <- rep("Unassigned", ncol(obj))
  rm(obj)
  save_spot_asset(path, paste0("F1b_before_", gsm), values, "Unassigned",
                  c(Unassigned = COL_UNCERTAIN), 74, 47, 0.47)
}, mixed_raw, names(mixed_raw))
resolved_assets <- Map(function(path, gsm) {
  obj <- readRDS(path)
  values <- as.character(obj$region_condition)
  rm(obj)
  save_spot_asset(path, paste0("F1b_resolved_", gsm), values, c("Sham", "Needle-pass"),
                  c(Sham = COL_SHAM, "Needle-pass" = COL_NEEDLE), 74, 47, 0.49)
}, mixed_labeled, names(mixed_labeled))
hero_columns <- lapply(names(mixed_raw), function(gsm) {
  pair <- wrap_plots(list(
    wrap_plots(list(strip_text_v2("Before assignment",6.8),img_panel_compact_v2(before_assets[[gsm]])),ncol=1,heights=c(.12,.88)),
    wrap_plots(list(strip_text_v2("Author-resolved",6.8),img_panel_compact_v2(resolved_assets[[gsm]])),ncol=1,heights=c(.12,.88))),ncol=2)
  wrap_plots(list(strip_text_v2(gsm),pair),ncol=1,heights=c(.13,.87))
})
hero_plate <- wrap_plots(hero_columns,ncol=2)
mixed_key <- ggplot(data.frame(
  x = 1:2, y = 1, label = c("Sham", "Needle-pass"),
  fill = c(COL_SHAM, COL_NEEDLE)
), aes(x, y)) +
  geom_point(aes(fill = fill), shape = 21, size = 2.8, stroke = 0.35, colour = COL_AXIS) +
  geom_text(aes(label = label), nudge_x = 0.13, hjust = 0, family = FONT_V2,
            size = 7.1 / GG_PT, colour = COL_AXIS) +
  scale_fill_identity() + coord_cartesian(xlim = c(0.82, 2.48), ylim = c(0.8, 1.2), expand = FALSE) +
  theme_void(base_family = FONT_V2)
hero_block <- wrap_plots(list(hero_plate, mixed_key), ncol = 1, heights = c(0.93, 0.07))
panel_b <- evidence_header_v2(
  hero_block, "b", "Mixed-capture resolution",
  "Author-provided barcode assignment; morphology does not define identity"
)

# c | Frozen six-tissue alignment plate.
expected_selection <- c(
  "Sham_Day3_rep1", "Sham_Day3_rep2", "MI_1HR_rep1",
  "MI_4HR_rep1", "MI_Day3_rep3", "MI_Day7_rep3"
)
selection <- selection[match(expected_selection, sample_key)]
if (anyNA(selection$sample_key) || !identical(selection$sample_key, expected_selection)) {
  stop("Figure 1 representative-tissue selection changed")
}
alignment_assets <- lapply(seq_len(nrow(selection)), function(i) {
  save_spot_asset(selection$qc_object[[i]], paste0("F1c_alignment_", selection$sample_key[[i]]),
                  "Retained spot", "Retained spot", c("Retained spot" = COL_DESC), 52, 36, 0.32)
})
alignment_labels <- c("Sham r1", "Sham r2", "MI 1 h",
                      "MI 4 h", "MI D3 r3", "MI D7 r3")
alignment_plate <- labeled_asset_grid(alignment_assets, alignment_labels, 6)
panel_c <- evidence_header_v2(
  alignment_plate, "c", "Tissue alignment",
  "Six fixed tissues; complete atlas in Figure S1"
)

# d | Retained-spot lollipop.
spot_counts[, sample_key := factor(sample_key, levels = rev(sample_order_v2))]
spot_counts[, role := factor(role_from_key_v2(as.character(sample_key)), levels = names(role_cols_v2))]
spot_count_plot <- ggplot(spot_counts, aes(retained_spots, sample_key)) +
  geom_segment(aes(x = 0, xend = retained_spots, yend = sample_key),
               linewidth = 0.42, colour = "#C7CCD1") +
  geom_point(aes(fill = role), shape = 21, size = 2.6, stroke = 0.35, colour = COL_AXIS) +
  geom_text(aes(label = scales::comma(retained_spots)), hjust = -0.18,
            family = FONT_V2, size = 7.1 / GG_PT, colour = COL_AXIS) +
  scale_fill_manual(values = role_cols_v2, guide = "none") +
  scale_y_discrete(labels = sample_labels_v2[rev(sample_order_v2)]) +
  scale_x_continuous(limits = c(0, 3300), breaks = c(0, 1000, 2000, 3000),
                     labels = scales::comma, expand = c(0, 0)) +
  labs(x = "Retained spots", y = NULL) + theme_editorial_v2() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(1, 2, 1, 1, unit = "mm"))
panel_d <- evidence_header_v2(spot_count_plot, "d", "Retained spots",
                              "All 25,548 tissue spots retained")

# e | Median and IQR of the frozen descriptive spot distributions.
qc_long <- melt(qc,
  id.vars = c("sample_key", "spot_key"),
  measure.vars = c("log10_nCount", "log10_nFeature", "percent.mt"),
  variable.name = "metric", value.name = "value"
)
qc_summary <- qc_long[is.finite(value), .(
  median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75)
), by = .(sample_key, metric)]
qc_summary[, role := factor(role_from_key_v2(sample_key), levels = names(role_cols_v2))]
qc_summary[, `:=`(
  sample_key = factor(sample_key, levels = sample_order_v2),
  metric = factor(metric,
    levels = c("log10_nCount", "log10_nFeature", "percent.mt"),
    labels = c("log10(UMI + 1)", "log10(features + 1)", "Mitochondrial\nreads (%)"))
)]
qc_plot <- ggplot(qc_summary, aes(sample_key, median, colour = role)) +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.55) +
  geom_point(size = 1.8) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = role_cols_v2, guide = "none") +
  scale_x_discrete(labels = setNames(c("Sh1", "Sh2", "D3r1", "D3r2", "D3r3",
                                      "D7r1", "D7r2", "D7r3", "1h", "4h", "N3", "N4"),
                                    sample_order_v2)) +
  labs(x = NULL, y = "Median [IQR]") + theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 6.2),
        panel.spacing.x = unit(3.0, "mm"), plot.margin = margin(1, 1, 1, 1, unit = "mm"))
panel_e <- evidence_header_v2(qc_plot, "e", "Descriptive QC summaries",
                              "Within-tissue median and IQR; no universal cutoff")

# f | Six-cell-type reference composition.
reference[, display := gsub("_", " ", celltype)]
reference[, display := factor(display, levels = display[order(log10_reference_cells)])]
reference[, colour := celltype_cols[celltype]]
reference_plot <- ggplot(reference, aes(log10_reference_cells, display)) +
  geom_segment(aes(x = 1, xend = log10_reference_cells, yend = display),
               linewidth = 0.45, colour = "#C7CCD1") +
  geom_point(aes(fill = colour), shape = 21, size = 2.7, stroke = 0.35, colour = COL_AXIS) +
  geom_text(aes(label = scales::comma(reference_cells),
                hjust = ifelse(reference_cells > 20000, 1.18, -0.16)),
            family = FONT_V2, size = 7.1 / GG_PT) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(1, 5.55), breaks = 1:5,
                     labels = function(x) paste0("10", c("¹", "²", "³", "⁴", "⁵")),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Reference nuclei (log10 view)", y = NULL) + theme_editorial_v2() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(1, 2, 1, 1, unit = "mm"))
panel_f <- evidence_header_v2(reference_plot, "f", "snRNA-seq reference",
                              "176,526 nuclei · Schwann n=28 (caution)")

# g | Common genes plus 6 × 12 marker-detection heatmap.
common_one <- unique(common[, .(exact_common_genes, fraction_spatial_genes_common,
                                fraction_reference_genes_common)])
if (nrow(common_one) != 1L) stop("Figure 1 common-gene summary is not invariant across tissues")
common_plot_data <- data.table(
  label = factor(c("Spatial genes", "Reference genes"), levels = c("Reference genes", "Spatial genes")),
  fraction = c(common_one$fraction_spatial_genes_common, common_one$fraction_reference_genes_common)
)
common_plot <- ggplot(common_plot_data, aes(fraction, label)) +
  geom_col(width = 0.52, fill = "#6F8FA8") +
  geom_text(aes(label = scales::percent(fraction, accuracy = 0.1)), hjust = -0.16,
            family = FONT_V2, size = 7.1 / GG_PT) +
  annotate("text", x = 0.5, y = 2.65, label = "22,490 exact\ncommon genes",
           family = FONT_V2, size = 7.2 / GG_PT, fontface = "bold", lineheight = 0.92) +
  scale_x_continuous(limits = c(0, 1.08), labels = scales::percent,
                     breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.55, 2.85), clip = "off") +
  labs(x = "Genes represented", y = NULL) + theme_editorial_v2() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(2, 2, 1, 1, unit = "mm"))

celltype_order <- c("Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage", "Pericyte", "Schwann_cell")
marker[, `:=`(
  sample_key = factor(sample_key, levels = sample_order_v2),
  celltype = factor(celltype, levels = rev(celltype_order))
)]
marker_heat <- ggplot(marker, aes(sample_key, celltype, fill = median_positive_fraction)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  scale_fill_gradient(low = "#F7F9FB", high = "#2D6598", limits = c(0, 1),
                      oob = scales::squish, name = "Detection\nfraction") +
  scale_x_discrete(labels = sample_labels_v2[sample_order_v2]) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = NULL, y = NULL) + theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 6.2),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", plot.margin = margin(1, 1, 1, 1, unit = "mm"))
compatibility_block <- wrap_plots(list(common_plot, marker_heat), ncol = 2, widths = c(0.27, 0.73))
panel_g <- evidence_header_v2(compatibility_block, "g", "Reference compatibility",
                              "Six prespecified broad-cell-type marker sets across 12 physical tissues")

# h | Low-ink analytical spot accounting.
accounting_plot <- ggplot() +
  annotate("rect", xmin = 0.2, xmax = 2.15, ymin = 0.23, ymax = 0.77,
           fill = "#F0F2F4", colour = "#7C858E", linewidth = 0.38) +
  annotate("rect", xmin = 3.2, xmax = 5.15, ymin = 0.23, ymax = 0.77,
           fill = "#EEE5EC", colour = "#9B6688", linewidth = 0.38) +
  annotate("rect", xmin = 6.25, xmax = 7.95, ymin = 0.53, ymax = 0.96,
           fill = "#E2F0ED", colour = COL_DESC, linewidth = 0.38) +
  annotate("rect", xmin = 6.25, xmax = 7.95, ymin = 0.04, ymax = 0.47,
           fill = "#ECEDEF", colour = COL_UNCERTAIN, linewidth = 0.38) +
  annotate("segment", x = 2.25, xend = 3.05, y = 0.5, yend = 0.5,
           colour = "#7C858E", linewidth = 0.38, arrow = arrow(length = unit(1.8, "mm"))) +
  annotate("segment", x = 5.25, xend = 6.08, y = 0.5, yend = 0.74,
           colour = "#7C858E", linewidth = 0.38, arrow = arrow(length = unit(1.8, "mm"))) +
  annotate("segment", x = 5.25, xend = 6.08, y = 0.5, yend = 0.25,
           colour = "#7C858E", linewidth = 0.38, arrow = arrow(length = unit(1.8, "mm"))) +
  annotate("text", x = 1.175, y = 0.5, label = "25,548 retained\ntissue spots",
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold", lineheight = 0.92) +
  annotate("text", x = 4.175, y = 0.5, label = "5,937 macrophage-rich\n23.2% of retained",
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold", lineheight = 0.92) +
  annotate("text", x = 7.10, y = 0.745, label = "4,143 assigned\n69.8% of selected",
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold", lineheight = 0.92) +
  annotate("text", x = 7.10, y = 0.255, label = "1,794 Uncertain\n30.2% of selected",
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold", lineheight = 0.92) +
  coord_cartesian(xlim = c(0, 8.15), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void(base_family = FONT_V2) + theme(plot.margin = margin(0, 1, 0, 1, unit = "mm"))
panel_h <- evidence_header_v2(accounting_plot, "h", "Spot accounting",
                              "Counts are descriptive pipeline units, not biological replicates", header_mm = 9)

figure1 <- wrap_plots(list(
  panel_a,
  panel_b,
  panel_c,
  wrap_plots(list(panel_d, panel_e), ncol = 2, widths = c(0.32, 0.68)),
  wrap_plots(list(panel_f, panel_g), ncol = 2, widths = c(0.28, 0.72)),
  panel_h
), ncol = 1, heights = c(0.9, 1.85, 1.15, 1.85, 1.8, 1.05)) &
  theme(plot.margin = margin(0.6, 0.8, 0.6, 0.8, unit = "mm"))

paths <- export_page_v2(figure1)
named_paths <- c(
  pdf = file.path(STAGE_V2, "figures", "main", "Figure_1.pdf"),
  svg = file.path(STAGE_V2, "figures", "main", "Figure_1.svg"),
  qa = file.path(STAGE_V2, "qa", "Figure_1_render.png")
)
if (any(file.exists(named_paths))) stop("Refusing to overwrite named Figure 1 staging outputs")
if (!all(file.rename(unname(paths), unname(named_paths)))) stop("Could not rename staged Figure 1 outputs")
if (any(file.info(named_paths)$size <= 1000)) stop("Named Figure 1 output is incomplete")
cat("FIGURE_1_STAGED_OK\n")
cat(paste(names(named_paths), named_paths, sep = "=", collapse = "\n"), "\n")

saveRDS(list(panel_a,panel_b,panel_c,panel_d,panel_e,panel_f,panel_g,panel_h),file.path(STAGE_V2,"panels.rds"))
