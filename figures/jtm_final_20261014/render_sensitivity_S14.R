#!/usr/bin/env Rscript

source("/home/zhenhanli/myocardial/scop/GSE214611_spatial/22_JTM_completion_20261014/figures_polished/script/figure_common_v2.R")

BUILD_ROOT <- Sys.getenv("REDESIGN_V2_BUILD_ROOT", unset = "")
if (!nzchar(BUILD_ROOT) || !dir.exists(BUILD_ROOT)) {
  stop("REDESIGN_V2_BUILD_ROOT must name an existing build root")
}
STAGE_V2 <- file.path(BUILD_ROOT, "Figure_S14")
if (file.exists(STAGE_V2)) stop("Figure 4 staging directory already exists: ", STAGE_V2)
ASSET_V2 <- file.path(STAGE_V2, "assets")

if (!exists("SpatialSpotPlot", envir = asNamespace("scop"), mode = "function", inherits = FALSE)) {
  stop("scop::SpatialSpotPlot is unavailable in the live runtime")
}

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

F4_SRC <- file.path(ROOT, "results", "figures", "Figure_4", "source_data")
ASSIGN_PATH <- file.path(ROOT, "metadata", "phase_07_neighborhood_assignments.tsv.gz")
SUMMARY_PATH <- file.path(ROOT, "metadata", "phase_07_niche_summary.tsv")
RECURRENCE_PATH <- file.path(ROOT, "metadata", "phase_07_niche_recurrence.tsv")
COHERENCE_PATH <- file.path(ROOT, "metadata", "phase_07_spatial_coherence.tsv")

niche_order <- c(
  "inflammatory-infarct niche",
  "macrophage-fibroblast remodeling niche",
  "vascular-repair niche",
  "remote resident niche",
  "Uncertain"
)
niche_labels <- setNames(c(
  "Inflammatory", "Mac–fibro remodeling", "Vascular repair",
  "Remote resident", "Uncertain"
), niche_order)
niche_cols <- c(
  "inflammatory-infarct niche" = "#EE6E79",
  "macrophage-fibroblast remodeling niche" = "#A48ACB",
  "vascular-repair niche" = "#39B7A5",
  "remote resident niche" = "#428BD0",
  "Uncertain" = "#6F6F6F",
  "Not_macrophage_rich" = "#D0D0D0"
)
niche_shapes <- c(
  "Not_macrophage_rich" = 1,
  "inflammatory-infarct niche" = 21,
  "macrophage-fibroblast remodeling niche" = 22,
  "vascular-repair niche" = 24,
  "remote resident niche" = 23,
  "Uncertain" = 25
)

a_source <- read_tsv_v2(file.path(F4_SRC, "Figure_4a_local_neighborhood.tsv"))
b_source <- read_tsv_v2(file.path(F4_SRC, "Figure_4b_evidence_breadth.tsv"))
c_source <- read_tsv_v2(file.path(F4_SRC, "Figure_4c_joint_feature_heatmap.tsv"))
assignments <- read_tsv_v2(ASSIGN_PATH)
niche_summary <- read_tsv_v2(SUMMARY_PATH)
recurrence <- read_tsv_v2(RECURRENCE_PATH)
coherence <- read_tsv_v2(COHERENCE_PATH)

if (nrow(assignments) != 5937L || uniqueN(assignments$center_barcode) < 1L ||
    anyDuplicated(assignments[, paste(sample_key, center_barcode, sep = "::")])) {
  stop("Figure 4 macrophage-rich-center contract changed")
}
if (sum(assignments$P7_niche == "Uncertain") != 1794L ||
    sum(assignments$P7_niche != "Uncertain") != 4143L ||
    !setequal(unique(assignments$P7_niche), niche_order) ||
    uniqueN(assignments$sample_key) != 12L) {
  stop("Figure 4 assignment contract changed")
}
if (nrow(a_source) != 7L || sum(a_source$is_center) != 1L || any(a_source$cross_slice_edge)) {
  stop("Figure 4 saved neighborhood is not one center plus six within-slice neighbors")
}
expected_denominators <- c(
  "inflammatory-infarct niche" = 1418L,
  "macrophage-fibroblast remodeling niche" = 1900L,
  "vascular-repair niche" = 685L,
  "remote resident niche" = 140L,
  "Uncertain" = 1794L
)
observed_denominators <- table(assignments$P7_niche)
if (!identical(as.integer(observed_denominators[names(expected_denominators)]),
               as.integer(expected_denominators))) {
  stop("Figure 4 niche denominators changed")
}
expected_recurrence <- c(
  "inflammatory-infarct niche" = 10L,
  "macrophage-fibroblast remodeling niche" = 7L,
  "vascular-repair niche" = 10L,
  "remote resident niche" = 9L,
  "Uncertain" = 10L
)
rec_lookup <- setNames(recurrence$primary_MI_slices, recurrence$niche)
if (nrow(recurrence) != 5L ||
    !identical(as.integer(rec_lookup[names(expected_recurrence)]), as.integer(expected_recurrence))) {
  stop("Figure 4 frozen non-auxiliary recurrence counts changed")
}
if (nrow(coherence) != 12L || uniqueN(coherence$sample_key) != 12L ||
    coherence[sample_key == "Sham_Day3_rep1", enrichment] >= 0) {
  stop("Figure 4 coherence contract changed or the frozen negative slice was lost")
}

mi6 <- c(
  "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
  "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3"
)
niche_objects <- setNames(
  file.path(ROOT, "objects", "niches", paste0(mi6, "_phase07.rds")), mi6
)
if (any(!file.exists(niche_objects))) stop("One or more formal-MI niche objects are missing")

# a | One saved within-slice neighborhood, rendered on its native H&E image.
day3_object <- readRDS(niche_objects[["MI_Day3_rep3"]])
a_source[, role := factor(role, levels = c("Direct neighbor", "Center"))]
neighborhood_cols <- c("Direct neighbor" = "#39B7A5", "Center" = "#EE6E79")
p_a_native <- suppressMessages(suppressWarnings(scop::SpatialSpotPlot(
  day3_object,
  plot.data = a_source,
  spot.by = "neighbor_barcode",
  color.by = "role",
  image = Seurat::Images(day3_object)[[1L]],
  overlay_image = TRUE,
  image.alpha = 0.18,
  crop = TRUE,
  show_axes = FALSE,
  pt.size = 3.2,
  pt.alpha = 0.92,
  stroke = 0.08,
  palcolor = unname(neighborhood_cols),
  legend.position = "none",
  verbose = FALSE
)))
p_a_native <- suppressMessages(p_a_native +
  scale_fill_manual(values = neighborhood_cols, drop = FALSE))
a_build <- ggplot_build(p_a_native)
stopifnot(identical(unname(a_build$plot$scales$get_scales("fill")$map(names(neighborhood_cols))),
                    unname(neighborhood_cols)))
a_layer <- which(vapply(a_build$data, function(x) {
  nrow(x) == 7L && all(c("x", "y") %in% names(x))
}, logical(1)))
if (length(a_layer) != 1L) stop("Could not recover the seven native spatial coordinates for panel a")
a_plot <- copy(a_source)
a_plot[, `:=`(x = a_build$data[[a_layer]]$x, y = a_build$data[[a_layer]]$y)]
a_center <- a_plot[is_center == TRUE]
a_edges <- a_plot[is_center == FALSE]
a_edges[, `:=`(center_x = a_center$x[[1L]], center_y = a_center$y[[1L]])]
p_a <- p_a_native +
  geom_segment(data = a_edges, aes(center_x, center_y, xend = x, yend = y),
               inherit.aes = FALSE, linewidth = 0.82, colour = "white", alpha = 0.90) +
  geom_segment(data = a_edges, aes(center_x, center_y, xend = x, yend = y),
               inherit.aes = FALSE, linewidth = 0.38, colour = "#68717A", alpha = 0.92) +
  geom_point(data = a_plot, aes(x, y, fill = role), inherit.aes = FALSE,
             shape = 21, size = 3.15, stroke = 0.28, colour = COL_AXIS) +
  theme_void(base_family = FONT_V2) +
  theme(plot.margin = margin(0, 0, 0, 0))
p_a <- p_a + guides(fill = guide_legend(title = NULL, ncol = 1,
                                        override.aes = list(size = 1.8))) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.1),
        legend.key.height = unit(2.5, "mm"), legend.key.width = unit(2.5, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0, "mm"))
a_asset <- save_spatial_asset_v2(p_a, "F4a_saved_neighborhood", 38, 38)
panel_a <- evidence_header_v2(p_a, "a", "Neighborhood unit",
                              "Center + 6 direct neighbors")
rm(day3_object, p_a_native, p_a, a_build)
gc(verbose = FALSE)

# b | Evidence breadth across all 5,937 frozen centers.
b_counts <- b_source[, .N, by = .(P7_niche, P7_evidence_domains)]
b_counts[, fraction := N / sum(N), by = P7_niche]
b_counts[, niche := factor(P7_niche, levels = niche_order, labels = unname(niche_labels))]
b_counts[, breadth := factor(P7_evidence_domains, levels = 0:3)]
b_labels <- setNames(
  paste0(c("Inflam.", "Mac–fibro", "Vascular", "Remote", "Uncertain"),
         "  n=", scales::comma(expected_denominators)),
  unname(niche_labels)
)
p_b <- ggplot(b_counts, aes(fraction, niche, fill = breadth)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("0" = "#D7DADD", "1" = "#B8D8D4", "2" = "#5AA69D", "3" = "#176F66"),
                    name = NULL, drop = FALSE) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_discrete(labels = b_labels) +
  labs(x = "Centers", y = NULL) + theme_editorial_v2() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = unit(2.5, "mm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0, 1, 1, 1, unit = "mm")) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
panel_b <- evidence_header_v2(p_b, "b", "Evidence breadth",
                              "Domains 0–3 per frozen center")

# c | Physical-slice-equal joint feature structure.
c_summary <- unique(c_source[, .(
  niche, feature, display_median_slice_mean, n_slices, display_summary
)])
c_records <- c_source[, .N, by = .(niche, feature)]
c_records <- merge(c_records, c_summary[, .(niche, feature, n_slices)],
                   by = c("niche", "feature"), all = TRUE)
if (nrow(c_summary) != 60L || any(!is.finite(c_summary$display_median_slice_mean)) ||
    any(c_records$N != c_records$n_slices) ||
    any(c_summary$display_summary != "median_of_per_slice_niche_means")) {
  stop("Figure 4 joint-feature summary is no longer a complete 12 x 5 slice-equal matrix")
}
feature_order <- c(
  "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage",
  "Resident Folr2", "Trem2 Spp1 remodeling", "Inflammatory Ccr2", "Antigen presentation",
  "Inflammation module", "ECM module", "Vascular module", "Cardiomyocyte-stress module"
)
c_limit <- max(abs(c_summary$display_median_slice_mean))
feature_domains <- list(
  "Cell composition" = feature_order[1:4],
  "Macrophage\nprograms" = feature_order[5:8],
  "Multi-gene\nmodules" = feature_order[9:12]
)
feature_display <- c(
  Cardiomyocyte = "CM", Endothelial = "Endo", Fibroblast = "Fibro", Macrophage = "Macro",
  `Resident Folr2` = "Folr2", `Trem2 Spp1 remodeling` = "Trem2/SPP1",
  `Inflammatory Ccr2` = "Ccr2", `Antigen presentation` = "Ag present.",
  `Inflammation module` = "Inflam.", `ECM module` = "ECM",
  `Vascular module` = "Vascular", `Cardiomyocyte-stress module` = "CM stress"
)
make_c_domain <- function(domain_name) {
  features <- feature_domains[[domain_name]]
  d <- c_summary[feature %in% features]
  d[, `:=`(x = match(niche, niche_order), y = 5L - match(feature, features))]
  ggplot(d, aes(x, y, fill = display_median_slice_mean)) +
    geom_tile(colour = "white", linewidth = 0.30) +
    geom_text(data = data.frame(x = 1:5, y = 4.63, label = c("I", "M", "V", "R", "U")),
              aes(x, y, label = label), inherit.aes = FALSE, family = FONT_V2,
              size = 5.3 / GG_PT, fontface = "bold", colour = COL_AXIS) +
    scale_fill_gradient2(low = "#397FC2", mid = "#FAFCFD", high = "#ED6A75",
                         midpoint = 0, limits = c(-c_limit, c_limit),
                         oob = scales::squish, guide = "none") +
    scale_x_continuous(breaks = NULL, limits = c(0.5, 5.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = 4:1, labels = unname(feature_display[features]),
                       limits = c(0.5, 4.95), expand = c(0, 0)) +
    labs(title = domain_name, x = NULL, y = NULL) + theme_editorial_v2() +
    theme(plot.title = element_text(size = 7.1, face = "bold", hjust = 0.5,
                                    margin = margin(b = 0.5, unit = "mm")),
          axis.text.y = element_text(size = 6, lineheight = 0.80, margin = margin(r = 1)),
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.margin = margin(0, 0.6, 0, 0, unit = "mm"))
}
c_domains <- lapply(names(feature_domains), make_c_domain)
c_bar <- continuous_bar_v2(c(-c_limit, c_limit),
                           c("#397FC2", "#FAFCFD", "#ED6A75"),
                           "Robust z", midpoint = 0, vertical = TRUE)
p_c <- wrap_plots(c(c_domains, list(c_bar)), ncol = 4, widths = c(0.27, 0.34, 0.25, 0.14))
panel_c <- evidence_header_v2(p_c, "c", "Joint niche feature structure",
                              "Slice-equal median · I/M/V/R/U follow frozen niche order")

# d | Six formal-MI maps from frozen Phase 7 objects.
niche_map_asset <- function(path, id) {
  obj <- readRDS(path)
  values <- as.character(obj$P7_niche)
  values[is.na(values)] <- "Not_macrophage_rich"
  obj$.v2_niche <- factor(values, levels = names(niche_shapes))
  p <- suppressMessages(suppressWarnings(scop::SpatialSpotPlot(
    obj,
    group.by = ".v2_niche",
    assay = "Spatial",
    layer = "counts",
    image = Seurat::Images(obj)[[1L]],
    plot_type = "point",
    overlay_image = TRUE,
    image.alpha = 0.18,
    crop = TRUE,
    show_axes = FALSE,
    pt.size = 0.44,
    pt.alpha = 0.92,
    stroke = 0.015,
    palcolor = unname(niche_cols[names(niche_shapes)]),
    legend.position = "none",
    legend.title = NULL,
    theme_use = "theme_blank",
    verbose = FALSE
  )))
  # SCOP point layers map categories to fill, not to palette position.
  p <- suppressMessages(p + scale_fill_manual(values = niche_cols,
                                              limits = names(niche_cols), drop = FALSE))
  mapped <- ggplot_build(p)$plot$scales$get_scales("fill")$map(names(niche_cols))
  stopifnot(identical(unname(mapped), unname(niche_cols)))
  if (!is.null(p$data) && all(c("x", "y", ".value") %in% names(p$data))) {
    p <- p + geom_point(data = p$data[p$data$.value != "Not_macrophage_rich", , drop = FALSE],
      aes(x, y, shape = .value), inherit.aes = FALSE, fill = NA, colour = COL_AXIS,
      size = 0.36, stroke = 0.055, alpha = 0.82, show.legend = FALSE) +
      scale_shape_manual(values = niche_shapes, drop = FALSE)
  }
  p <- suppressMessages(suppressWarnings(p + coord_fixed())) +
    theme_void(base_family = FONT_V2) + theme(plot.margin = margin(0, 0, 0, 0))
  p <- p & theme(legend.position = "none")
  out <- save_spatial_asset_v2(p, id, 53, 43)
  rm(obj, p)
  gc(verbose = FALSE)
  out
}
d_assets <- lapply(mi6, function(key) {
  niche_map_asset(niche_objects[[key]], paste0("F4d_", key))
})
d_plate <- spatial_plate_v2(
  d_assets, nrow = 2, ncol = 3,
  col_labels = c("r1", "r2", "r3"),
  row_labels = c("Day 3 MI", "Day 7 MI"),
  row_label_width = 0.075, header_height = 0.07
)
legend_df <- data.frame(
  x = rep(1:3, 2), y = rep(2:1, each = 3),
  key = names(niche_shapes),
  label = c("Not macrophage-rich", unname(niche_labels))
)
d_legend <- ggplot(legend_df, aes(x, y)) +
  geom_point(aes(shape = key, fill = key), size = 2.5, stroke = 0.38, colour = COL_AXIS) +
  geom_text(aes(label = label), nudge_x = 0.13, hjust = 0,
            family = FONT_V2, size = 7.1 / GG_PT, colour = COL_AXIS) +
  scale_shape_manual(values = niche_shapes, guide = "none") +
  scale_fill_manual(values = niche_cols, guide = "none") +
  coord_cartesian(xlim = c(0.88, 3.83), ylim = c(0.65, 2.35), expand = FALSE) +
  theme_void(base_family = FONT_V2)
d_block <- wrap_plots(list(d_plate, d_legend), ncol = 1, heights = c(0.90, 0.10))
panel_d <- evidence_header_v2(d_block, "d", "Recurrent niches across formal MI tissues",
                              "All neighborhoods are slice-local; six independent physical tissues")

# e | Prevalence matrix plus frozen n/10 non-auxiliary recurrence annotation.
grid_e <- CJ(sample_key = sample_order_v2, niche = niche_order, unique = TRUE)
grid_e <- merge(
  grid_e,
  niche_summary[, .(sample_key, niche, fraction = fraction_of_macrophage_rich_centers)],
  by = c("sample_key", "niche"), all.x = TRUE, sort = FALSE
)
grid_e[is.na(fraction), fraction := 0]
grid_e[, `:=`(
  x = match(sample_key, sample_order_v2),
  y = match(niche, niche_order)
)]
sample_role <- data.table(
  x = seq_along(sample_order_v2),
  role = role_from_key_v2(sample_order_v2)
)
rec_ann <- data.table(
  x = 13.0,
  y = seq_along(niche_order),
  label = paste0(expected_recurrence[niche_order], "/10")
)
p_e <- ggplot(grid_e, aes(x, y, fill = fraction)) +
  geom_tile(width = 0.92, height = 0.84, colour = "white", linewidth = 0.25) +
  geom_segment(data = sample_role,
               aes(x = x - 0.43, xend = x + 0.43, y = 5.56, yend = 5.56, colour = role),
               inherit.aes = FALSE, linewidth = 3.1, lineend = "butt") +
  geom_text(data = rec_ann, aes(x, y, label = label), inherit.aes = FALSE,
            family = FONT_V2, size = 5.3 / GG_PT, fontface = "bold", colour = COL_AXIS) +
  annotate("text", x = 13.0, y = 5.82, label = "n/10",
           family = FONT_V2, size = 7.1 / GG_PT, fontface = "bold") +
  scale_fill_gradient(low = "#F4F6F7", high = "#334E68", limits = c(0, 1),
                      oob = scales::squish, name = "Center\nfraction") +
  scale_colour_manual(values = role_cols_v2, guide = "none") +
  scale_x_continuous(breaks = seq_along(sample_order_v2),
                     labels = sample_labels_v2[sample_order_v2], limits = c(0.5, 13.75),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(niche_order),
                     labels = unname(niche_labels[niche_order]), limits = c(0.5, 6.45),
                     expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 52, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", plot.margin = margin(0, 1, 1, 1, unit = "mm"))
panel_e <- evidence_header_v2(p_e, "e", "Cross-tissue prevalence and recurrence",
                              "Matrix shows all 12 tissues; n/10 excludes the two auxiliary Needle-pass tissues")

# f | Frozen assignment-confidence composition.
certainty <- assignments[, .N, by = .(sample_key, P7_confidence)]
certainty[, fraction := N / sum(N), by = sample_key]
confidence_order <- c("High", "Moderate", "Low")
if (!all(unique(certainty$P7_confidence) %in% confidence_order)) {
  stop("Unexpected Figure 4 assignment-confidence category")
}
confidence_cols <- c(High = "#237A57", Moderate = "#78C7DB", Low = "#F4777F")
certainty[, `:=`(
  sample_key = factor(sample_key, levels = sample_order_v2),
  confidence = factor(P7_confidence, levels = confidence_order)
)]
p_f <- ggplot(certainty, aes(sample_key, fraction, fill = confidence)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = confidence_cols, drop = FALSE, name = NULL) +
  scale_x_discrete(labels = sample_labels_v2[sample_order_v2]) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Center fraction") + theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1),
        legend.position = "bottom", legend.direction = "horizontal",
        plot.margin = margin(0, 1, 1, 1, unit = "mm")) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
panel_f <- evidence_header_v2(p_f, "f", "Assignment certainty",
                              "All confidence strata retained")

# g | Observed-minus-permuted coherence, with the frozen negative value intact.
coherence[, role := factor(role_from_key_v2(sample_key), levels = names(role_cols_v2))]
coherence[, sample_key := factor(sample_key, levels = sample_order_v2)]
p_g <- ggplot(coherence, aes(sample_key, enrichment, colour = role)) +
  geom_hline(yintercept = 0, linewidth = 0.34, colour = "#7C858E") +
  geom_segment(aes(xend = sample_key, y = 0, yend = enrichment), linewidth = 0.60) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = role_cols_v2, guide = "none") +
  scale_x_discrete(labels = sample_labels_v2[sample_order_v2]) +
  labs(x = NULL, y = "Coherence Δ") + theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1),
        plot.margin = margin(0, 1, 1, 5, unit = "mm"))
panel_g <- evidence_header_v2(p_g, "g", "Within-slice coherence",
                              "One descriptive value per physical tissue; negative values retained")

# h | Four frozen modules summarized once per physical tissue.
module_fields <- c(
  "module_inflammation", "module_ECM", "module_vascular", "module_cardiomyocyte_stress"
)
modules <- melt(assignments,
  id.vars = "sample_key", measure.vars = module_fields,
  variable.name = "module", value.name = "score"
)
modules <- modules[is.finite(score), .(score = median(score)), by = .(sample_key, module)]
modules[, role := factor(role_from_key_v2(sample_key), levels = names(role_cols_v2))]
modules[, `:=`(
  sample_key = factor(sample_key, levels = sample_order_v2),
  module = factor(module, levels = module_fields,
                  labels = c("Inflammation", "ECM", "Vascular", "CM stress"))
)]
p_h <- ggplot(modules, aes(sample_key, score, colour = role)) +
  geom_hline(yintercept = 0, linewidth = 0.30, colour = "#A2A8AE") +
  geom_point(size = 1.8) +
  facet_wrap(~module, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = role_cols_v2, guide = "none") +
  scale_x_discrete(labels = setNames(c("Sh1", "Sh2", "D3r1", "D3r2", "D3r3",
                                      "D7r1", "D7r2", "D7r3", "1h", "4h", "N3", "N4"),
                                    sample_order_v2)) +
  labs(x = NULL, y = "Median center score") + theme_editorial_v2() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 5.5),
        panel.spacing.x = unit(2.8, "mm"), plot.margin = margin(0, 1, 1, 1, unit = "mm"))
panel_h <- evidence_header_v2(p_h, "h", "Frozen niche-module behavior",
                              "One median per physical tissue; module-specific y scales")

figure4 <- wrap_plots(list(
  wrap_plots(list(panel_a, panel_b, panel_c), ncol = 3, widths = c(0.23, 0.25, 0.52)),
  panel_d,
  panel_e,
  wrap_plots(list(panel_f, panel_g), ncol = 2, widths = c(0.44, 0.56)),
  panel_h
), ncol = 1, heights = c(1.10, 3.10, 1.25, 0.95, 1.10)) &
  theme(plot.margin = margin(0.6, 0.8, 0.6, 0.8, unit = "mm"))

paths <- export_page_v2(figure4)
named_paths <- c(
  pdf = file.path(STAGE_V2, "figures", "main", "Figure_S14.pdf"),
  svg = file.path(STAGE_V2, "figures", "main", "Figure_S14.svg"),
  qa = file.path(STAGE_V2, "qa", "Figure_S14_render.png")
)
if (any(file.exists(named_paths))) stop("Refusing to overwrite named Figure 4 staging outputs")
if (!all(file.rename(unname(paths), unname(named_paths)))) stop("Could not rename staged Figure 4 outputs")
if (any(file.info(named_paths)$size <= 1000)) stop("Named Figure 4 output is incomplete")
cat("FIGURE_4_STAGED_OK\n")
cat(paste(names(named_paths), named_paths, sep = "=", collapse = "\n"), "\n")
