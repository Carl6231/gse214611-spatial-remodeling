#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"

options(stringsAsFactors = FALSE, warn = 1)

paths <- list(
  workflow = file.path(root, "docs", "WORKFLOW.md"),
  status = file.path(root, "docs", "STATUS.md"),
  decisions = file.path(root, "docs", "DECISIONS.md"),
  design = file.path(root, "metadata", "sample_design.tsv"),
  objects = file.path(root, "objects", "raw"),
  audits = file.path(root, "results", "object_audit"),
  summary = file.path(root, "metadata", "spatial_object_summary_raw.tsv"),
  report = file.path(root, "reports", "phase_01_report.md"),
  log = file.path(root, "logs", "phase_01.log"),
  session = file.path(root, "results", "object_audit", "session_info.txt")
)

for (directory in c(paths$objects, paths$audits, dirname(paths$report), dirname(paths$log))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

sink(paths$log, append = FALSE, split = TRUE)
on.exit(sink(), add = TRUE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  flush.console()
}

format_time <- function(value) format(value, "%Y-%m-%d %H:%M:%S %z")

atomic_write <- function(lines, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("could not atomically write ", path)
  }
}

atomic_write_table <- function(value, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  write.table(value, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("could not atomically write ", path)
  }
}

replace_field <- function(lines, range, field, value) {
  prefix <- paste0("- ", field, ":")
  index <- range[startsWith(lines[range], prefix)]
  if (length(index) != 1L) stop("STATUS field is missing or ambiguous: ", field)
  lines[[index]] <- paste0(prefix, if (nzchar(value)) paste0(" ", value) else "")
  lines
}

update_status <- function(phase_status, start_time, finish_time = "", validation,
                          unresolved, next_phase, outputs) {
  lines <- readLines(paths$status, warn = FALSE, encoding = "UTF-8")
  phase0_start <- match("## Phase 0", lines)
  phase1_start <- match("## Phase 1", lines)
  phase2_start <- match("## Phase 2", lines)
  if (anyNA(c(phase0_start, phase1_start, phase2_start))) {
    stop("STATUS phase headings are incomplete")
  }
  top_range <- seq_len(phase0_start - 1L)
  phase1_range <- phase1_start:(phase2_start - 1L)
  scripts <- "scripts/phase_01_build_objects.R; scripts/launch_phase_01.sh"

  top_values <- list(
    current_phase = "Phase 1",
    phase_status = phase_status,
    start_time = start_time,
    finish_time = finish_time,
    scripts = scripts,
    outputs = outputs,
    validation = validation,
    unresolved_issues = unresolved,
    next_phase = next_phase
  )
  for (field in names(top_values)) {
    lines <- replace_field(lines, top_range, field, top_values[[field]])
  }

  phase_values <- list(
    status = phase_status,
    start_time = start_time,
    finish_time = finish_time,
    scripts = scripts,
    outputs = outputs,
    validation = validation,
    unresolved_issues = unresolved,
    next_phase = next_phase
  )
  for (field in names(phase_values)) {
    lines <- replace_field(lines, phase1_range, field, phase_values[[field]])
  }
  atomic_write(lines, paths$status)
}

append_decision <- function() {
  decision_text <- paste(
    "|", format(Sys.Date(), "%Y-%m-%d"),
    "| Phase 1 按当前 Seurat 5.5.0 函数定义使用 Read10X_h5、Read10X_Image 和 CreateSeuratObject 逐样本构建对象；递归识别标准 spatial/ 与同层布局，不重排原始文件。",
    "| 当前 Load10X_Spatial 固定查找 spatial/，而 8 个 TAR 的空间文件与 H5 同层；手动传入真实 image.dir 可保持原始目录只读。",
    "| Phase 1 | 否 |"
  )
  decisions <- readLines(paths$decisions, warn = FALSE, encoding = "UTF-8")
  if (!any(grepl("Phase 1 按当前 Seurat 5.5.0", decisions, fixed = TRUE))) {
    atomic_write(c(decisions, decision_text), paths$decisions)
  }
}

write_report <- function(phase_status, start_time, finish_time, summary = NULL,
                         unresolved = "none") {
  completed <- if (is.null(summary)) 0L else sum(summary$overall_status == "PASS")
  table_lines <- if (!is.null(summary) && nrow(summary)) {
    c(
      "",
      "| GSM | features | spots | barcode match | image bounds | readback | status |",
      "|---|---:|---:|---|---|---|---|",
      vapply(seq_len(nrow(summary)), function(i) {
        paste0(
          "| ", summary$GSM[[i]], " | ", summary$features[[i]], " | ",
          summary$spots[[i]], " | ", summary$barcode_setequal[[i]], " | ",
          summary$coordinates_within_image[[i]], " | ",
          summary$rds_readback_status[[i]], " | ", summary$overall_status[[i]], " |"
        )
      }, character(1))
    )
  } else {
    c("", "- 尚无成功对象记录。")
  }

  lines <- c(
    "# Phase 1 report",
    "",
    paste0("- Phase status: **", phase_status, "**"),
    paste0("- Start time: ", start_time),
    paste0("- Finish time: ", finish_time),
    "",
    "## 实际完成的任务",
    "",
    "- 重新读取 WORKFLOW.md、STATUS.md 和 DECISIONS.md，并确认 Phase 0 PASS。",
    "- 使用服务器当前 Seurat 5.5.0 的已核实函数逐 GSM 加载 filtered matrix、坐标和组织图像。",
    "- 每个 GSM 独立构建原始 Spatial assay Seurat 对象，并添加 GSM、sample_id、condition、timepoint、replicate、mixed_capture。",
    "- 检查 count barcodes 与 image/coordinate barcodes、图像边界、对象回读和 metadata 常量。",
    "- 为每个样本导出 H&E + spots 对齐图；未做 QC、过滤、归一化、整合或去卷积。",
    "",
    "## 输出文件",
    "",
    "- `objects/raw/GSMxxxx_raw.rds`",
    "- `metadata/spatial_object_summary_raw.tsv`",
    "- `results/object_audit/GSMxxxx/alignment_lowres.png`",
    "- `results/object_audit/GSMxxxx/object_audit.tsv`",
    "- `results/object_audit/GSMxxxx/object_structure.txt`",
    "- `results/object_audit/session_info.txt`",
    "- `logs/phase_01.log`",
    "- `reports/phase_01_report.md`",
    "",
    "## 验证结果",
    "",
    paste0("- 成功且通过对象自检：", completed, "/10。"),
    table_lines,
    "",
    "## 失败或未验证内容",
    "",
    if (identical(unresolved, "none")) "- 无。" else paste0("- ", unresolved),
    "",
    "## 风险",
    "",
    "- 本阶段仅验证对象结构、barcode/坐标一致性、图像边界和导出叠加图；没有进行生物学 QC。",
    "- GSM6613080 和 GSM6613081 仍为 mixed_capture，未将全部 spots 标记为 Sham，也未拆分组织区域。",
    "- spots 仍不是生物学重复。",
    "",
    "## 下一阶段入口条件",
    "",
    if (phase_status == "PASS") {
      "- Phase 2 为 NOT_STARTED，必须在单独执行中重新读取协议后才能进入。"
    } else {
      "- Phase 1 需要修复或人工检查；Phase 2 不得启动。"
    },
    "- 本次执行在 Phase 1 后停止。"
  )
  atomic_write(lines, paths$report)
}

required <- c(paths$workflow, paths$status, paths$decisions, paths$design)
missing <- required[!file.exists(required)]
if (length(missing)) stop("required file missing: ", paste(missing, collapse = ", "))

workflow_before <- readLines(paths$workflow, warn = FALSE, encoding = "UTF-8")
status_before <- readLines(paths$status, warn = FALSE, encoding = "UTF-8")
decisions_before <- readLines(paths$decisions, warn = FALSE, encoding = "UTF-8")
log_msg("control files re-read:", length(workflow_before), length(status_before), length(decisions_before))

phase0_start <- match("## Phase 0", status_before)
phase1_start <- match("## Phase 1", status_before)
phase2_start <- match("## Phase 2", status_before)
stopifnot(!anyNA(c(phase0_start, phase1_start, phase2_start)))
phase0_block <- status_before[phase0_start:(phase1_start - 1L)]
phase1_block <- status_before[phase1_start:(phase2_start - 1L)]
if (!any(grepl("- status: PASS", phase0_block, fixed = TRUE))) {
  stop("Phase 0 is not PASS")
}
if (!any(grepl("- unresolved_issues: none", phase0_block, fixed = TRUE))) {
  stop("Phase 0 has unresolved issues")
}
if (!any(grepl("- status: NOT_STARTED", phase1_block, fixed = TRUE))) {
  stop("Phase 1 is not NOT_STARTED")
}

start_time_value <- Sys.time()
start_time <- format_time(start_time_value)
running_outputs <- "pending"
update_status(
  phase_status = "RUNNING",
  start_time = start_time,
  validation = "Phase 1 object construction is running",
  unresolved = "pending",
  next_phase = "Phase 1",
  outputs = running_outputs
)
log_msg("Phase 1 status set to RUNNING")

result <- tryCatch({
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
  })

  log_msg("R:", R.version.string)
  log_msg("Seurat:", as.character(packageVersion("Seurat")))
  log_msg("SeuratObject:", as.character(packageVersion("SeuratObject")))
  for (fn in c("Read10X_h5", "Read10X_Image", "CreateSeuratObject", "SpatialDimPlot")) {
    if (!exists(fn, envir = asNamespace("Seurat"), mode = "function")) {
      stop("required Seurat function not found: ", fn)
    }
  }

  design <- read.delim(paths$design, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  expected_gsms <- sprintf("GSM%d", 6613080:6613089)
  if (!identical(design$GSM, expected_gsms)) stop("sample design GSM order/scope mismatch")

  writeLines(capture.output(sessionInfo()), paths$session)
  summary_rows <- vector("list", nrow(design))

  for (i in seq_len(nrow(design))) {
    gsm <- design$GSM[[i]]
    log_msg(gsm, "construction started")
    unpacked <- file.path(root, "spatial_raw", gsm, "unpacked")
    if (!dir.exists(unpacked)) stop(gsm, ": unpacked directory missing")

    h5 <- list.files(
      unpacked, pattern = "^filtered_feature_bc_matrix\\.h5$",
      recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    positions <- list.files(
      unpacked, pattern = "^tissue_positions(_list)?\\.csv$",
      recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    scalefactors <- list.files(
      unpacked, pattern = "^scalefactors_json\\.json$",
      recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    hires <- list.files(
      unpacked, pattern = "^tissue_hires_image\\.png$",
      recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    lowres <- list.files(
      unpacked, pattern = "^tissue_lowres_image\\.png$",
      recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    if (any(lengths(list(h5, positions, scalefactors, hires, lowres)) != 1L)) {
      stop(gsm, ": expected exactly one H5, positions, scalefactors, hires and lowres file")
    }
    spatial_dir <- dirname(positions[[1]])
    if (!all(dirname(c(scalefactors, hires, lowres)) == spatial_dir)) {
      stop(gsm, ": spatial payloads are not in one directory")
    }

    counts <- Read10X_h5(h5[[1]], use.names = TRUE, unique.features = TRUE)
    if (is.list(counts)) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else if (length(counts) == 1L) {
        counts <- counts[[1L]]
      } else {
        stop(gsm, ": H5 contains multiple feature types without a Gene Expression entry")
      }
    }
    if (anyDuplicated(colnames(counts))) stop(gsm, ": duplicated count barcodes")
    if (anyDuplicated(rownames(counts))) stop(gsm, ": duplicated feature names after import")

    image <- Read10X_Image(
      image.dir = spatial_dir,
      image.name = "tissue_lowres_image.png",
      assay = "Spatial",
      slice = gsm,
      filter.matrix = TRUE,
      image.type = "VisiumV2"
    )
    count_barcodes <- colnames(counts)
    coordinate_barcodes <- Cells(image)
    barcode_setequal <- setequal(count_barcodes, coordinate_barcodes)
    missing_coordinates <- setdiff(count_barcodes, coordinate_barcodes)
    extra_coordinates <- setdiff(coordinate_barcodes, count_barcodes)
    if (!barcode_setequal) {
      stop(
        gsm, ": barcode mismatch; missing coordinates=", length(missing_coordinates),
        ", extra coordinates=", length(extra_coordinates)
      )
    }

    object <- CreateSeuratObject(
      counts = counts,
      assay = "Spatial",
      project = gsm,
      min.cells = 0,
      min.features = 0
    )
    image <- image[Cells(object)]
    object[[gsm]] <- image

    mixed_capture <- tolower(as.character(design$mixed_capture[[i]])) == "true"
    object$GSM <- gsm
    object$sample_id <- design$sample_id[[i]]
    object$condition <- design$condition[[i]]
    object$timepoint <- design$timepoint[[i]]
    object$replicate <- as.integer(design$replicate[[i]])
    object$mixed_capture <- mixed_capture
    object@misc$phase_01_provenance <- list(
      source_h5 = h5[[1]],
      source_spatial_dir = spatial_dir,
      source_positions = positions[[1]],
      Seurat_version = as.character(packageVersion("Seurat")),
      construction_time = format_time(Sys.time()),
      qc_performed = FALSE,
      normalization_performed = FALSE,
      integration_performed = FALSE,
      deconvolution_performed = FALSE
    )

    required_metadata <- c("GSM", "sample_id", "condition", "timepoint", "replicate", "mixed_capture")
    metadata_ok <- all(required_metadata %in% colnames(object[[]])) &&
      all(object$GSM == gsm) &&
      all(object$sample_id == design$sample_id[[i]]) &&
      all(object$condition == design$condition[[i]]) &&
      all(object$timepoint == design$timepoint[[i]]) &&
      all(object$replicate == as.integer(design$replicate[[i]])) &&
      all(object$mixed_capture == mixed_capture)
    counts_layer_ok <- "counts" %in% Layers(object[["Spatial"]])
    image_names <- Images(object)
    image_ok <- identical(image_names, gsm)

    coordinates <- GetTissueCoordinates(object[[gsm]], scale = "lowres")
    image_array <- GetImage(object[[gsm]], mode = "raw")
    coordinate_count_ok <- nrow(coordinates) == ncol(object)
    bounds_ok <- all(
      is.finite(coordinates$x), is.finite(coordinates$y),
      coordinates$x >= 0, coordinates$x <= dim(image_array)[2],
      coordinates$y >= 0, coordinates$y <= dim(image_array)[1]
    )
    if (!all(c(metadata_ok, counts_layer_ok, image_ok, coordinate_count_ok, bounds_ok))) {
      stop(gsm, ": object structure or coordinate/image validation failed")
    }

    audit_dir <- file.path(paths$audits, gsm)
    dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
    object_path <- file.path(paths$objects, paste0(gsm, "_raw.rds"))
    object_temp <- paste0(object_path, ".tmp.", Sys.getpid())
    saveRDS(object, object_temp)
    if (!file.rename(object_temp, object_path)) {
      unlink(object_temp)
      stop(gsm, ": could not atomically save object")
    }

    readback <- readRDS(object_path)
    readback_ok <- inherits(readback, "Seurat") &&
      nrow(readback) == nrow(object) &&
      ncol(readback) == ncol(object) &&
      identical(Images(readback), gsm) &&
      all(required_metadata %in% colnames(readback[[]]))
    if (!readback_ok) stop(gsm, ": saved RDS failed readback validation")

    plot <- SpatialDimPlot(
      object = object,
      group.by = "GSM",
      images = gsm,
      cols = "#C81D25",
      crop = FALSE,
      label = FALSE,
      pt.size.factor = 1.3,
      alpha = c(0.80, 0.80),
      image.alpha = 1,
      image.scale = "lowres"
    ) +
      ggtitle(paste0(gsm, " — H&E + filtered spots")) +
      labs(subtitle = paste0(ncol(object), " spots; barcode match and image bounds PASS")) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        legend.position = "none"
      )
    plot_path <- file.path(audit_dir, "alignment_lowres.png")
    ggsave(plot_path, plot = plot, width = 7.2, height = 7.2, dpi = 300, bg = "white")
    if (!file.exists(plot_path) || file.info(plot_path)$size <= 0) {
      stop(gsm, ": alignment plot was not written")
    }

    sha256 <- system2("sha256sum", shQuote(object_path), stdout = TRUE)
    sha256 <- sub(" .*", "", sha256[[1]])
    object_audit <- data.frame(
      GSM = gsm,
      source_h5 = h5[[1]],
      source_spatial_dir = spatial_dir,
      features = nrow(object),
      spots = ncol(object),
      coordinate_spots = nrow(coordinates),
      count_barcodes = length(count_barcodes),
      coordinate_barcodes = length(coordinate_barcodes),
      barcode_setequal = barcode_setequal,
      missing_coordinate_barcodes = length(missing_coordinates),
      extra_coordinate_barcodes = length(extra_coordinates),
      image_name = gsm,
      image_class = class(object[[gsm]])[[1]],
      image_height = dim(image_array)[1],
      image_width = dim(image_array)[2],
      coordinate_x_min = min(coordinates$x),
      coordinate_x_max = max(coordinates$x),
      coordinate_y_min = min(coordinates$y),
      coordinate_y_max = max(coordinates$y),
      coordinates_within_image = bounds_ok,
      metadata_status = metadata_ok,
      counts_layer_status = counts_layer_ok,
      rds_readback_status = readback_ok,
      rds_size_bytes = file.info(object_path)$size,
      rds_sha256 = sha256,
      alignment_plot = plot_path,
      overall_status = "PASS",
      stringsAsFactors = FALSE
    )
    atomic_write_table(object_audit, file.path(audit_dir, "object_audit.tsv"))
    structure_lines <- c(
      capture.output(print(object)),
      "",
      paste0("Assays: ", paste(Assays(object), collapse = ",")),
      paste0("Spatial Layers: ", paste(Layers(object[["Spatial"]]), collapse = ",")),
      paste0("Images: ", paste(Images(object), collapse = ",")),
      paste0("Metadata: ", paste(colnames(object[[]]), collapse = ",")),
      paste0("Barcode setequal: ", barcode_setequal),
      paste0("Coordinates within lowres image: ", bounds_ok),
      paste0("RDS readback: ", readback_ok)
    )
    atomic_write(structure_lines, file.path(audit_dir, "object_structure.txt"))

    summary_rows[[i]] <- object_audit
    current_summary <- do.call(rbind, summary_rows[seq_len(i)])
    atomic_write_table(current_summary, paths$summary)
    log_msg(
      gsm, "PASS features=", nrow(object), " spots=", ncol(object),
      " image=", class(object[[gsm]])[[1]], " rds_bytes=", file.info(object_path)$size
    )
    rm(counts, image, object, readback, plot, image_array, coordinates)
    invisible(gc())
  }

  summary <- do.call(rbind, summary_rows)
  if (nrow(summary) != 10L || !all(summary$overall_status == "PASS")) {
    stop("Phase 1 summary did not contain 10 PASS samples")
  }
  object_files <- list.files(paths$objects, pattern = "^GSM661308[0-9]_raw\\.rds$", full.names = TRUE)
  plot_files <- list.files(paths$audits, pattern = "^alignment_lowres\\.png$", recursive = TRUE, full.names = TRUE)
  if (length(object_files) != 10L) stop("expected 10 raw RDS files, found ", length(object_files))
  if (length(plot_files) != 10L) stop("expected 10 alignment plots, found ", length(plot_files))
  list(summary = summary, object_files = object_files, plot_files = plot_files)
}, error = function(error) error)

finish_time <- format_time(Sys.time())
output_text <- paste(
  "objects/raw/GSMxxxx_raw.rds",
  "metadata/spatial_object_summary_raw.tsv",
  "results/object_audit/GSMxxxx/",
  "results/object_audit/session_info.txt",
  "logs/phase_01.log",
  "reports/phase_01_report.md",
  sep = "; "
)

if (inherits(result, "error")) {
  unresolved <- conditionMessage(result)
  log_msg("Phase 1 FAIL:", unresolved)
  partial_summary <- if (file.exists(paths$summary)) {
    read.delim(paths$summary, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
  write_report("FAIL", start_time, finish_time, partial_summary, unresolved)
  update_status(
    phase_status = "FAIL",
    start_time = start_time,
    finish_time = finish_time,
    validation = paste0("Phase 1 failed: ", unresolved),
    unresolved = unresolved,
    next_phase = "Phase 1 remediation or human review; Phase 2 must not start",
    outputs = output_text
  )
  log_msg("STOP: Phase 2 was not executed")
  quit(save = "no", status = 1L)
}

summary <- result$summary
validation <- paste0(
  "objects=", length(result$object_files), "/10; alignment_plots=", length(result$plot_files),
  "/10; barcode_setequal=", sum(summary$barcode_setequal), "/10",
  "; coordinate_bounds=", sum(summary$coordinates_within_image), "/10",
  "; rds_readback=", sum(summary$rds_readback_status), "/10"
)
write_report("PASS", start_time, finish_time, summary, "none")
append_decision()
update_status(
  phase_status = "PASS",
  start_time = start_time,
  finish_time = finish_time,
  validation = validation,
  unresolved = "none",
  next_phase = "Phase 2 (NOT_STARTED; requires a separate explicit execution)",
  outputs = output_text
)

workflow_after <- readLines(paths$workflow, warn = FALSE, encoding = "UTF-8")
status_after <- readLines(paths$status, warn = FALSE, encoding = "UTF-8")
phase1_start <- match("## Phase 1", status_after)
phase2_start <- match("## Phase 2", status_after)
phase3_start <- match("## Phase 3", status_after)
checks <- c(
  workflow_reread = length(workflow_after) > 0L,
  status_reread = length(status_after) > 0L,
  phase1_pass = any(grepl("- status: PASS", status_after[phase1_start:(phase2_start - 1L)], fixed = TRUE)),
  phase2_not_started = any(grepl("- status: NOT_STARTED", status_after[phase2_start:(phase3_start - 1L)], fixed = TRUE))
)
log_msg("final consistency:", paste(names(checks), checks, sep = "=", collapse = ";"))
log_msg("Phase 1 finished with status PASS")
log_msg("STOP: Phase 2 was not executed")
if (!all(checks)) quit(save = "no", status = 1L)
