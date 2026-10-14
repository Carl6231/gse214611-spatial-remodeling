suppressPackageStartupMessages(suppressWarnings({
  library(Seurat)
  library(SeuratObject)
  library(scop)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(grid)
}))

ROOT <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
V1 <- file.path(ROOT, "21_scop_expansion", "PUBLICATION_FIGURE_REDESIGN")
V2 <- file.path(ROOT, "22_JTM_completion_20261014", "figures_polished")
FIG_V2 <- file.path(V2, "figures", "main")
QA_V2 <- file.path(V2, "qa")
STAGE_V2 <- file.path(V2, "tmp", "figure5_build")
ASSET_V2 <- file.path(STAGE_V2, "assets")

FONT_V2 <- "Nimbus Sans"
font_path_v2 <- systemfonts::match_fonts(FONT_V2)$path
if (!length(font_path_v2) || !nzchar(font_path_v2) || !file.exists(font_path_v2)) {
  stop("Nimbus Sans is not available through systemfonts")
}

MM_TO_IN <- 1 / 25.4
PT_TO_MM <- 25.4 / 72
GG_PT <- 72.27 / 25.4
COL_AXIS <- "#20252B"
COL_SHAM <- "#4A90D9"
COL_D3 <- "#F4777F"
COL_D7 <- "#9A78C6"
COL_DESC <- "#36B6A4"
COL_NEEDLE <- "#78C7DB"
COL_UNCERTAIN <- "#AEB7C2"
COL_CONTRACTILE <- "#397FC2"
COL_ZERO <- "#FAFCFD"
COL_FIBROTIC <- "#ED6A75"

sample_order_v2 <- c(
  "Sham_Day3_rep1", "Sham_Day3_rep2",
  "MI_Day3_rep1", "MI_Day3_rep2", "MI_Day3_rep3",
  "MI_Day7_rep1", "MI_Day7_rep2", "MI_Day7_rep3",
  "MI_1HR_rep1", "MI_4HR_rep1",
  "NeedlePass_Day3_rep3", "NeedlePass_Day3_rep4"
)
sample_labels_v2 <- setNames(c(
  "Sham r1", "Sham r2", "MI D3 r1", "MI D3 r2", "MI D3 r3",
  "MI D7 r1", "MI D7 r2", "MI D7 r3", "MI 1 h", "MI 4 h",
  "Needle r3", "Needle r4"
), sample_order_v2)
role_cols_v2 <- c(
  Sham = COL_SHAM, `MI Day 3` = COL_D3, `MI Day 7` = COL_D7,
  Descriptive = COL_DESC, `Needle-pass` = COL_NEEDLE
)

read_tsv_v2 <- function(path) {
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  data.table::fread(path, na.strings = c("NA", ""), data.table = TRUE)
}

theme_editorial_v2 <- function(base_size = 7.1) {
  theme_classic(base_family = FONT_V2, base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = COL_AXIS),
      axis.ticks = element_line(linewidth = 0.30, colour = COL_AXIS),
      axis.text = element_text(size = 7.1, colour = COL_AXIS),
      axis.title = element_text(size = 7.3, colour = COL_AXIS),
      strip.text = element_text(size = 7.2, face = "bold", colour = COL_AXIS),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.title = element_text(size = 7.1, face = "bold", colour = COL_AXIS),
      legend.text = element_text(size = 7.1, colour = COL_AXIS),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(3.0, "mm"),
      panel.spacing = unit(3, "mm"),
      axis.title.x = element_text(margin=margin(t=3)),
      axis.title.y = element_text(margin=margin(r=3)),
      plot.margin = margin(1.2, 1.2, 1.2, 1.2, unit = "mm")
    )
}

strip_text_v2 <- function(label, size = 7.2, face = "bold", hjust = 0.5) {
  patchwork::wrap_elements(full = textGrob(
    label,
    x = unit(hjust, "npc"), y = unit(0.5, "npc"),
    just = c(if (hjust < 0.5) "left" else if (hjust > 0.5) "right" else "center", "center"),
    gp = gpar(fontfamily = FONT_V2, fontsize = size, fontface = face, col = COL_AXIS)
  ))
}

evidence_header_v2 <- function(plot, letter, title, subtitle = NULL, header_mm = if (is.null(subtitle)) 5 else 9) {
  header <- patchwork::wrap_elements(full = grobTree(
    textGrob(letter, x = unit(0, "npc"), y = unit(0.62, "npc"),
             just = c("left", "center"),
             gp = gpar(fontfamily = FONT_V2, fontsize = 10, fontface = "bold", col = COL_AXIS)),
    textGrob(title, x = unit(4.5, "mm"), y = unit(0.62, "npc"),
             just = c("left", "center"),
             gp = gpar(fontfamily = FONT_V2, fontsize = 7.7, fontface = "bold", col = COL_AXIS)),
    if (!is.null(subtitle)) textGrob(subtitle, x = unit(4.5, "mm"), y = unit(0.08, "npc"),
             just = c("left", "bottom"),
             gp = gpar(fontfamily = FONT_V2, fontsize = 7.1, col = "#555D66")) else nullGrob()
  ))
  wrap_plots(list(header, wrap_plots(list(plot))), ncol = 1,
             heights = unit.c(unit(header_mm, "mm"), unit(1, "null")))
}

img_panel_v2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Missing rendered spatial panel: ", path)
  im <- png::readPNG(path)
  ratio <- dim(im)[1] / dim(im)[2]
  ggplot() +
    annotation_raster(im, xmin = 0, xmax = 1, ymin = 0, ymax = ratio) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, ratio), expand = FALSE, clip = "on") +
    theme_void(base_family = FONT_V2) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

img_panel_compact_v2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Missing rendered spatial panel: ", path)
  im <- png::readPNG(path)
  ratio <- dim(im)[1] / dim(im)[2]
  patchwork::wrap_elements(full = grobTree(
    rectGrob(gp = gpar(fill = "white", col = NA)),
    rasterGrob(im, width = unit(1 / ratio, "snpc"), height = unit(1, "snpc"), interpolate = TRUE)
  ))
}

spatial_plate_v2 <- function(paths, nrow, ncol, col_labels, row_labels = NULL,
                             row_label_width = 0.085, header_height = 0.10) {
  if (length(paths) != nrow * ncol) stop("spatial_plate_v2 path count does not match layout")
  if (length(col_labels) != ncol) stop("spatial_plate_v2 column-label count mismatch")
  if (!is.null(row_labels) && length(row_labels) != nrow) stop("spatial_plate_v2 row-label count mismatch")

  header <- wrap_plots(lapply(col_labels, strip_text_v2), ncol = ncol)
  rows <- lapply(seq_len(nrow), function(i) {
    idx <- ((i - 1L) * ncol + 1L):(i * ncol)
    panels <- wrap_plots(lapply(paths[idx], img_panel_v2), ncol = ncol)
    if (is.null(row_labels)) return(wrap_plots(list(panels)))
    wrap_plots(list(strip_text_v2(row_labels[i], size = 7.2, hjust = 1), panels),
               ncol = 2, widths = c(row_label_width, 1))
  })
  body <- wrap_plots(rows, ncol = 1)
  if (is.null(row_labels)) {
    wrap_plots(list(header, body), ncol = 1, heights = c(header_height, 1))
  } else {
    header_full <- wrap_plots(list(plot_spacer(), header), ncol = 2,
                              widths = c(row_label_width, 1))
    wrap_plots(list(header_full, body), ncol = 1, heights = c(header_height, 1))
  }
}

status_key_v2 <- function(label = "NEGATIVE_COMPLETE") {
  ggplot() +
    annotate("point", x = 0.08, y = 0.5, shape = 21, size = 2.3,
             stroke = 0.45, colour = "#68717D", fill = "#EEF0F2") +
    annotate("text", x = 0.14, y = 0.5, label = label, hjust = 0,
             family = FONT_V2, size = 7.1 / GG_PT, colour = COL_AXIS) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = FONT_V2) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

save_spatial_asset_v2 <- function(plot, id, width_mm, height_mm) {
  dir.create(ASSET_V2, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(ASSET_V2, paste0(id, ".png"))
  if (file.exists(path)) stop("Refusing to overwrite staged spatial asset: ", path)
  ragg::agg_png(path, width = width_mm * MM_TO_IN, height = height_mm * MM_TO_IN,
                units = "in", res = 300, background = "white")
  print(plot)
  dev.off()
  if (!file.exists(path) || file.info(path)$size <= 1000) stop("Spatial asset export failed: ", id)
  path
}

spatial_continuous_asset_v2 <- function(object_path, values, id, limits,
                                        palette, midpoint = NULL,
                                        width_mm = 56, height_mm = 42,
                                        point_size = 0.52) {
  obj <- readRDS(object_path)
  values <- values[colnames(obj)]
  if (length(values) != ncol(obj) || !all(is.finite(values) | is.na(values))) {
    stop("Spatial value contract failed for ", basename(object_path), " / ", id)
  }
  p <- suppressMessages(suppressWarnings(scop::SpatialSpotPlot(
      obj, values = values, image = Images(obj)[[1]], plot_type = "point",
      crop = TRUE, overlay_image = TRUE, image.alpha = 0.18,
      pt.alpha = 0.92, pt.size = point_size, stroke = 0,
      palcolor = palette, legend.position = "none", legend.title = NULL,
      theme_use = "theme_blank"
    )))
  p <- suppressMessages(suppressWarnings(if (length(palette) == 3L && !is.null(midpoint)) {
    p + scale_colour_gradient2(low = palette[1], mid = palette[2], high = palette[3],
                               midpoint = midpoint, limits = limits, oob = scales::squish)
  } else {
    p + scale_colour_gradientn(colours = palette, limits = limits, oob = scales::squish)
  }))
  p <- suppressMessages(suppressWarnings(p + coord_fixed())) +
    theme_void(base_family = FONT_V2) +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  out <- save_spatial_asset_v2(p, id, width_mm, height_mm)
  rm(obj, p)
  gc(verbose = FALSE)
  out
}

continuous_bar_v2 <- function(limits, colours, title, midpoint = NULL, vertical = TRUE) {
  if (vertical) {
    d <- data.frame(x = 1, y = seq(limits[1], limits[2], length.out = 256))
    p <- ggplot(d, aes(x, y, fill = y)) + geom_raster()
    if (length(colours) == 3L && !is.null(midpoint)) {
      p <- p + scale_fill_gradient2(low = colours[1], mid = colours[2], high = colours[3], midpoint = midpoint,
                                    limits = limits, oob = scales::squish, guide = "none")
    } else {
      p <- p + scale_fill_gradientn(colours = colours, limits = limits, oob = scales::squish, guide = "none")
    }
    p + scale_x_continuous(NULL, breaks = NULL) +
      scale_y_continuous(title, position = "right", breaks = pretty(limits, n = 3)) +
      coord_cartesian(expand = FALSE) + theme_editorial_v2() +
      theme(axis.line.y.left = element_blank(), axis.line.x = element_blank(),
            axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
            axis.title.y.right = element_text(angle = 90, size = 7.1),
            plot.margin = margin(2, 1, 2, 1, unit = "mm"))
  } else {
    d <- data.frame(x = seq(limits[1], limits[2], length.out = 256), y = 1)
    p <- ggplot(d, aes(x, y, fill = x)) + geom_raster()
    if (length(colours) == 3L && !is.null(midpoint)) {
      p <- p + scale_fill_gradient2(low = colours[1], mid = colours[2], high = colours[3], midpoint = midpoint,
                                    limits = limits, oob = scales::squish, guide = "none")
    } else {
      p <- p + scale_fill_gradientn(colours = colours, limits = limits, oob = scales::squish, guide = "none")
    }
    p + scale_y_continuous(NULL, breaks = NULL) +
      scale_x_continuous(title, breaks = pretty(limits, n = 3)) +
      coord_cartesian(expand = FALSE) + theme_editorial_v2() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.x = element_text(size = 7.1), plot.margin = margin(0, 1, 1, 1, unit = "mm"))
  }
}

export_page_v2 <- function(plot) {
  dir.create(file.path(STAGE_V2, "figures", "main"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(STAGE_V2, "qa"), recursive = TRUE, showWarnings = FALSE)
  dir.create(ASSET_V2, recursive = TRUE, showWarnings = FALSE)
  pdf <- file.path(STAGE_V2, "figures", "main", "Figure_5.pdf")
  svg <- file.path(STAGE_V2, "figures", "main", "Figure_5.svg")
  qa <- file.path(STAGE_V2, "qa", "Figure_5_render.png")
  if (any(file.exists(c(pdf, svg, qa)))) stop("Refusing to overwrite staged Figure 5 outputs")

  cairo_pdf(pdf, width = 183 * MM_TO_IN, height = 225 * MM_TO_IN,
            family = "sans", onefile = TRUE)
  print(plot)
  dev.off()
  svglite::svglite(svg, width = 183 * MM_TO_IN, height = 225 * MM_TO_IN,
                   bg = "white", system_fonts = list(sans = FONT_V2))
  print(plot)
  dev.off()
  rsvg::rsvg_png(svg, qa, width = 2161, height = 2657)
  paths <- c(pdf = pdf, svg = svg, qa = qa)
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 1000)) stop("Figure 5 staged export is incomplete")
  paths
}
