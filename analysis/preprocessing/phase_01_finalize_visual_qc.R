#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"

paths <- list(
  workflow = file.path(root, "docs", "WORKFLOW.md"),
  status = file.path(root, "docs", "STATUS.md"),
  summary = file.path(root, "metadata", "spatial_object_summary_raw.tsv"),
  visual = file.path(root, "results", "object_audit", "visual_qc.tsv"),
  report = file.path(root, "reports", "phase_01_report.md"),
  log = file.path(root, "logs", "phase_01_visual_qc.log")
)

sink(paths$log, append = FALSE, split = TRUE)
on.exit(sink(), add = TRUE)

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")

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

required <- unlist(paths[c("workflow", "status", "summary", "report")])
missing <- required[!file.exists(required)]
if (length(missing)) stop("missing required file(s): ", paste(missing, collapse = ", "))

summary <- read.delim(paths$summary, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
expected <- sprintf("GSM%d", 6613080:6613089)
stopifnot(
  identical(summary$GSM, expected),
  all(summary$overall_status == "PASS"),
  all(summary$barcode_setequal),
  all(summary$coordinates_within_image),
  all(summary$rds_readback_status)
)

plot_paths <- file.path(root, "results", "object_audit", expected, "alignment_lowres.png")
stopifnot(all(file.exists(plot_paths)), all(file.info(plot_paths)$size > 0))

notes <- rep(
  "filtered spots follow tissue morphology; no global rotation, flip, offset or clipping observed",
  length(expected)
)
notes[expected == "GSM6613080"] <- paste(
  "multiple tissue pieces visible in one mixed capture;",
  "filtered spots follow tissue morphology; mixed_capture retained"
)
notes[expected == "GSM6613081"] <- paste(
  "multiple tissue pieces visible in one mixed capture;",
  "filtered spots follow tissue morphology; mixed_capture retained"
)
notes[expected == "GSM6613087"] <- paste(
  "low-resolution H&E has a dark cyan cast;",
  "spot grid still follows tissue boundaries with no global transform error"
)

visual <- data.frame(
  GSM = expected,
  plot_path = plot_paths,
  spots_follow_tissue = TRUE,
  no_global_rotation_or_flip = TRUE,
  no_global_offset = TRUE,
  no_edge_clipping = TRUE,
  visual_status = "PASS",
  notes = notes,
  stringsAsFactors = FALSE
)
atomic_write_table(visual, paths$visual)

finish_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
report <- readLines(paths$report, warn = FALSE, encoding = "UTF-8")
marker <- "<!-- PHASE_01_VISUAL_QC -->"
marker_index <- match(marker, report)
if (!is.na(marker_index)) report <- report[seq_len(marker_index - 1L)]
finish_index <- which(startsWith(report, "- Finish time:"))
if (length(finish_index)) report[[finish_index[[1L]]]] <- paste0("- Finish time: ", finish_time)

visual_lines <- c(
  "",
  marker,
  "",
  "## 实际导出图像视觉自检",
  "",
  "- 已查看 10 张 `alignment_lowres.png` 联系图，并单独放大检查 GSM6613087。",
  "- 10/10：filtered spots 跟随组织轮廓，未见全局旋转、翻转、平移或边缘裁切错误。",
  "- GSM6613080/81 显示同一 capture 中的多个组织块；继续保持 mixed_capture，不在 Phase 1 拆分或赋予 Sham 区域标签。",
  "- GSM6613087 的低分辨率 H&E 底图偏青暗，但 spots 与组织边界一致；记录为图像色调风险，不构成配准失败。",
  "- 逐样本视觉结论：`results/object_audit/visual_qc.tsv`。",
  "",
  "## 最终 Phase 1 结论",
  "",
  "- 对象结构、barcode/坐标、图像边界、RDS 回读和实际 H&E+spots 叠加图均为 10/10 PASS。",
  "- Phase 2 保持 NOT_STARTED；本次执行在 Phase 1 后停止。"
)
atomic_write(c(report, visual_lines), paths$report)

status <- readLines(paths$status, warn = FALSE, encoding = "UTF-8")
phase0_start <- match("## Phase 0", status)
phase1_start <- match("## Phase 1", status)
phase2_start <- match("## Phase 2", status)
phase3_start <- match("## Phase 3", status)
stopifnot(!anyNA(c(phase0_start, phase1_start, phase2_start, phase3_start)))
top_range <- seq_len(phase0_start - 1L)
phase1_range <- phase1_start:(phase2_start - 1L)

replace_one <- function(lines, range, prefix, value) {
  index <- range[startsWith(lines[range], prefix)]
  if (length(index) != 1L) stop("STATUS field missing or ambiguous: ", prefix)
  lines[[index]] <- paste0(prefix, " ", value)
  lines
}

scripts_value <- "scripts/phase_01_build_objects.R; scripts/launch_phase_01.sh; scripts/phase_01_finalize_visual_qc.R"
outputs_value <- paste(
  "objects/raw/GSMxxxx_raw.rds",
  "metadata/spatial_object_summary_raw.tsv",
  "results/object_audit/GSMxxxx/",
  "results/object_audit/visual_qc.tsv",
  "results/object_audit/session_info.txt",
  "logs/phase_01.log",
  "logs/phase_01_visual_qc.log",
  "reports/phase_01_report.md",
  sep = "; "
)
validation_value <- paste(
  "objects=10/10; alignment_plots=10/10; barcode_setequal=10/10;",
  "coordinate_bounds=10/10; rds_readback=10/10; visual_alignment=10/10"
)

for (range in list(top_range, phase1_range)) {
  status <- replace_one(status, range, "- finish_time:", finish_time)
  status <- replace_one(status, range, "- scripts:", scripts_value)
  status <- replace_one(status, range, "- outputs:", outputs_value)
  status <- replace_one(status, range, "- validation:", validation_value)
  status <- replace_one(status, range, "- unresolved_issues:", "none")
}
atomic_write(status, paths$status)

workflow_after <- readLines(paths$workflow, warn = FALSE, encoding = "UTF-8")
status_after <- readLines(paths$status, warn = FALSE, encoding = "UTF-8")
checks <- c(
  workflow_reread = length(workflow_after) > 0L,
  status_reread = length(status_after) > 0L,
  phase1_pass = any(grepl("- status: PASS", status_after[phase1_range], fixed = TRUE)),
  phase2_not_started = any(grepl(
    "- status: NOT_STARTED",
    status_after[phase2_start:(phase3_start - 1L)],
    fixed = TRUE
  )),
  visual_qc_complete = nrow(visual) == 10L && all(visual$visual_status == "PASS")
)
log_msg("visual QC written:", paths$visual)
log_msg("report/status refreshed at", finish_time)
log_msg("final consistency:", paste(names(checks), checks, sep = "=", collapse = ";"))
log_msg("STOP: Phase 2 was not executed")
if (!all(checks)) quit(save = "no", status = 1L)
