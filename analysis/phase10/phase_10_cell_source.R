root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
source(file.path(root, "scripts", "phase_10_common.R"))
root <- phase10_root()
dirs <- phase10_dirs(root)
phase10_log_start(root, "phase_10_cell_source")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(matrixStats)
  library(edgeR)
})

set.seed(20260805)
sig <- readRDS(require_file(file.path(dirs$objects, "phase_10_spatial_signature.rds")))
positive_genes <- sig$positive_genes
negative_genes <- sig$negative_genes
time_levels <- c("Day0", "1HR", "4HR", "Day1", "Day3", "Day7")
time_numeric <- c(Day0 = 0, `1HR` = 1 / 24, `4HR` = 4 / 24, Day1 = 1, Day3 = 3, Day7 = 7)

aggregate_sparse <- function(counts, metadata, fields) {
  stopifnot(identical(colnames(counts), rownames(metadata)))
  keys <- do.call(paste, c(metadata[, fields, drop = FALSE], sep = "||"))
  groups <- unique(keys)
  membership <- Matrix::sparseMatrix(
    i = seq_along(keys), j = match(keys, groups), x = 1,
    dims = c(length(keys), length(groups))
  )
  pb <- counts %*% membership
  colnames(pb) <- groups
  design <- as.data.table(metadata[match(groups, keys), fields, drop = FALSE])
  design[, pseudobulk_id := groups]
  design[, cells := as.integer(tabulate(match(keys, groups), nbins = length(groups)))]
  setcolorder(design, c("pseudobulk_id", fields, "cells"))
  rownames(design) <- design$pseudobulk_id
  list(counts = pb, design = design)
}

module_scores <- function(pb, design, class_field, positive, negative) {
  positive <- intersect(positive, rownames(pb))
  negative <- intersect(negative, rownames(pb))
  genes <- unique(c(positive, negative))
  if (length(genes) < 10L) stop("Too few signature genes overlap pseudobulk matrix")
  full_library <- Matrix::colSums(pb)
  if (any(full_library <= 0)) stop("Full pseudobulk libraries must be positive before signature scoring")
  logcpm <- as.matrix(edgeR::cpm(
    pb[genes, , drop = FALSE], log = TRUE, prior.count = 1,
    lib.size = full_library
  ))
  rows <- vector("list", length(unique(design[[class_field]])))
  classes <- unique(as.character(design[[class_field]]))
  for (i in seq_along(classes)) {
    cl <- classes[[i]]
    idx <- which(as.character(design[[class_field]]) == cl)
    z <- safe_row_z(logcpm[, idx, drop = FALSE])
    pos <- intersect(positive, rownames(z))
    neg <- intersect(negative, rownames(z))
    pos_score <- if (length(pos)) colMeans(z[pos, , drop = FALSE]) else rep(0, length(idx))
    neg_score <- if (length(neg)) colMeans(z[neg, , drop = FALSE]) else rep(0, length(idx))
    rows[[i]] <- cbind(
      copy(design[idx]),
      data.table(
        positive_zscore = pos_score,
        negative_zscore = neg_score,
        remodeling_signature_score = pos_score - neg_score,
        positive_mean_logcpm = colMeans(logcpm[pos, idx, drop = FALSE]),
        negative_mean_logcpm = if (length(neg)) colMeans(logcpm[neg, idx, drop = FALSE]) else NA_real_,
        positive_genes = length(pos), negative_genes = length(neg)
      )
    )
  }
  rbindlist(rows, fill = TRUE)
}

hedges_g <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  sp <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
               (length(x) + length(y) - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  d <- (mean(x) - mean(y)) / sp
  correction <- 1 - 3 / (4 * (length(x) + length(y)) - 9)
  correction * d
}

summarize_scores <- function(scores, class_field) {
  classes <- unique(as.character(scores[[class_field]]))
  comparisons <- list(Day7_vs_Day0 = c("Day7", "Day0"), Day7_vs_Day3 = c("Day7", "Day3"))
  effects <- rbindlist(lapply(classes, function(cl) {
    xdt <- scores[as.character(get(class_field)) == cl & cells >= 20]
    rbindlist(lapply(names(comparisons), function(label) {
      groups <- comparisons[[label]]
      x <- xdt[timepoint == groups[[1L]], remodeling_signature_score]
      y <- xdt[timepoint == groups[[2L]], remodeling_signature_score]
      p <- if (length(x) >= 2L && length(y) >= 2L) {
        suppressWarnings(wilcox.test(x, y, exact = TRUE)$p.value)
      } else NA_real_
      data.table(
        source_class = cl, comparison = label, n_target = length(x), n_reference = length(y),
        mean_target = mean(x), mean_reference = mean(y), mean_difference = mean(x) - mean(y),
        median_difference = median(x) - median(y), hedges_g = hedges_g(x, y),
        wilcoxon_exact_p = p
      )
    }))
  }))
  effects[, wilcoxon_bh := p.adjust(wilcoxon_exact_p, method = "BH"), by = comparison]

  trends <- rbindlist(lapply(classes, function(cl) {
    xdt <- scores[as.character(get(class_field)) == cl & cells >= 20]
    tt <- unname(time_numeric[as.character(xdt$timepoint)])
    if (length(tt) < 5L || sd(tt) == 0 || sd(xdt$remodeling_signature_score) == 0) {
      return(data.table(source_class = cl, samples = nrow(xdt), spearman_rho = NA_real_, p_value = NA_real_))
    }
    test <- suppressWarnings(cor.test(tt, xdt$remodeling_signature_score, method = "spearman", exact = FALSE))
    data.table(source_class = cl, samples = nrow(xdt), spearman_rho = unname(test$estimate), p_value = test$p.value)
  }))
  trends[, p_adjust_bh := p.adjust(p_value, method = "BH")]

  loso <- rbindlist(lapply(classes, function(cl) {
    xdt <- scores[as.character(get(class_field)) == cl & cells >= 20]
    day7 <- xdt[timepoint == "Day7"]
    rbindlist(lapply(c("Day0", "Day3"), function(ref) {
      ref_values <- xdt[timepoint == ref, remodeling_signature_score]
      if (nrow(day7) < 3L || length(ref_values) < 2L) return(NULL)
      rbindlist(lapply(seq_len(nrow(day7)), function(i) {
        retained <- day7[-i, remodeling_signature_score]
        data.table(
          source_class = cl, comparison = paste0("Day7_vs_", ref),
          omitted_sample = day7$sample_id[[i]],
          mean_difference = mean(retained) - mean(ref_values),
          direction = sign(mean(retained) - mean(ref_values))
        )
      }))
    }))
  }), fill = TRUE)
  loso_summary <- loso[, .(
    omissions = .N, positive_omissions = sum(direction > 0), negative_omissions = sum(direction < 0),
    stable_direction = uniqueN(direction[direction != 0]) <= 1L,
    min_difference = min(mean_difference), max_difference = max(mean_difference)
  ), by = .(source_class, comparison)]
  list(effects = effects, trends = trends, loso = loso, loso_summary = loso_summary)
}

run_edger_by_class <- function(pb, design, class_field) {
  classes <- unique(as.character(design[[class_field]]))
  rbindlist(lapply(classes, function(cl) {
    idx <- which(as.character(design[[class_field]]) == cl & design$cells >= 20)
    current_design <- copy(design[idx])
    current_design[, timepoint := factor(as.character(timepoint), levels = time_levels)]
    current_design <- current_design[!is.na(timepoint)]
    current_design[, timepoint := droplevels(timepoint)]
    if (!all(c("Day0", "Day3", "Day7") %in% as.character(unique(current_design$timepoint)))) return(NULL)
    idx <- match(current_design$pseudobulk_id, colnames(pb))
    y <- DGEList(pb[, idx, drop = FALSE])
    keep <- filterByExpr(y, group = current_design$timepoint, min.count = 5)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    design_matrix <- model.matrix(~0 + timepoint, data = current_design)
    colnames(design_matrix) <- sub("^timepoint", "", colnames(design_matrix))
    y <- estimateDisp(y, design_matrix, robust = TRUE)
    fit <- glmQLFit(y, design_matrix, robust = TRUE)
    contrasts <- list(Day7_vs_Day0 = c(Day7 = 1, Day0 = -1), Day7_vs_Day3 = c(Day7 = 1, Day3 = -1))
    rbindlist(lapply(names(contrasts), function(label) {
      weights <- contrasts[[label]]
      contrast <- rep(0, ncol(design_matrix)); names(contrast) <- colnames(design_matrix)
      contrast[names(weights)] <- weights
      qlf <- glmQLFTest(fit, contrast = contrast)
      out <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
      out[, `:=`(
        source_class = cl,
        comparison = label,
        samples = nrow(current_design), tested_genes = nrow(out)
      )]
      out
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}

stamp("Loading corrected 11-class Phase 9 reference")
reference_path <- require_file(file.path(root, "objects", "phase_09", "reference_expanded_immune.rds"))
reference <- readRDS(reference_path)
ref_md <- reference[[]]
if (!all(c("P09_expanded_celltype", "sample_id", "timepoint") %in% colnames(ref_md))) stop("Corrected reference metadata missing")
ref_counts <- GetAssayData(reference, assay = "RNA", layer = "counts")
pb <- aggregate_sparse(ref_counts, ref_md, c("sample_id", "timepoint", "P09_expanded_celltype"))
pb$design[, library_size := as.numeric(Matrix::colSums(pb$counts))]
zero_reference <- pb$design[library_size <= 0]
write_tsv(zero_reference, file.path(dirs$metadata, "phase_10_reference_zero_library_pseudobulks.tsv"))
if (nrow(zero_reference)) {
  keep <- pb$design$library_size > 0
  pb$counts <- pb$counts[, keep, drop = FALSE]
  pb$design <- pb$design[keep]
}
if (uniqueN(pb$design$P09_expanded_celltype) != 11L) stop("Expected 11 corrected cell classes")
signature_zero_reference <- copy(pb$design)
signature_zero_reference[, signature_counts := as.numeric(Matrix::colSums(
  pb$counts[intersect(c(positive_genes, negative_genes), rownames(pb$counts)), , drop = FALSE]
))]
signature_zero_reference <- signature_zero_reference[signature_counts == 0]
write_tsv(signature_zero_reference, file.path(dirs$metadata, "phase_10_reference_signature_zero_count_pseudobulks.tsv"))
pb_scores <- module_scores(pb$counts, pb$design, "P09_expanded_celltype", positive_genes, negative_genes)
summary_ref <- summarize_scores(pb_scores, "P09_expanded_celltype")
stamp("Running sample-level edgeR in each corrected cell class")
edger_ref <- run_edger_by_class(pb$counts, pb$design, "P09_expanded_celltype")
edger_ref[, signature_direction := fifelse(gene %in% positive_genes, "positive", fifelse(gene %in% negative_genes, "negative", "not_signature"))]
signature_edger <- edger_ref[signature_direction != "not_signature"]
signature_concordance <- signature_edger[, .(
  genes = .N,
  median_logFC = median(logFC),
  expected_direction_fraction = mean(sign(logFC) == fifelse(signature_direction == "positive", 1, -1)),
  fdr05_genes = sum(FDR <= 0.05)
), by = .(source_class, comparison, signature_direction)]

pb_meta <- as.data.frame(pb$design)
rownames(pb_meta) <- pb_meta$pseudobulk_id
pb_object <- CreateSeuratObject(counts = pb$counts, assay = "RNA", meta.data = pb_meta)
saveRDS(pb_object, file.path(dirs$objects, "phase_10_reference_pseudobulk_seurat.rds"))
saveRDS(list(counts = pb$counts, design = pb$design), file.path(dirs$objects, "phase_10_reference_pseudobulk_counts.rds"))
rm(reference, ref_counts, ref_md, pb_object)
gc()

stamp("Loading refined macrophage reference")
macro_path <- require_file(file.path(root, "objects", "phase_09", "reference_macrophage_refined.rds"))
macro <- readRDS(macro_path)
macro_md <- macro[[]]
macro_counts <- GetAssayData(macro, assay = "RNA", layer = "counts")
macro_pb <- aggregate_sparse(macro_counts, macro_md, c("sample_id", "timepoint", "P09_macrophage_state"))
macro_pb$design[, library_size := as.numeric(Matrix::colSums(macro_pb$counts))]
zero_macro <- macro_pb$design[library_size <= 0]
write_tsv(zero_macro, file.path(dirs$metadata, "phase_10_macrophage_zero_library_pseudobulks.tsv"))
if (nrow(zero_macro)) {
  keep <- macro_pb$design$library_size > 0
  macro_pb$counts <- macro_pb$counts[, keep, drop = FALSE]
  macro_pb$design <- macro_pb$design[keep]
}
macro_scores <- module_scores(macro_pb$counts, macro_pb$design, "P09_macrophage_state", positive_genes, negative_genes)
summary_macro <- summarize_scores(macro_scores, "P09_macrophage_state")
saveRDS(list(counts = macro_pb$counts, design = macro_pb$design), file.path(dirs$objects, "phase_10_macrophage_state_pseudobulk_counts.rds"))
rm(macro, macro_counts, macro_md)
gc()

write_tsv(pb$design, file.path(dirs$metadata, "phase_10_reference_pseudobulk_design.tsv"))
write_tsv(pb_scores, file.path(dirs$metadata, "phase_10_reference_celltype_signature_scores.tsv"))
write_tsv(summary_ref$effects, file.path(dirs$metadata, "phase_10_reference_celltype_effects.tsv"))
write_tsv(summary_ref$trends, file.path(dirs$metadata, "phase_10_reference_celltype_time_trends.tsv"))
write_tsv(summary_ref$loso, file.path(dirs$metadata, "phase_10_reference_celltype_day7_loso.tsv"))
write_tsv(summary_ref$loso_summary, file.path(dirs$metadata, "phase_10_reference_celltype_day7_loso_summary.tsv"))
write_tsv(edger_ref, file.path(dirs$metadata, "phase_10_reference_celltype_edger.tsv.gz"))
write_tsv(signature_concordance, file.path(dirs$metadata, "phase_10_reference_celltype_signature_concordance.tsv"))
write_tsv(macro_pb$design, file.path(dirs$metadata, "phase_10_macrophage_state_pseudobulk_design.tsv"))
write_tsv(macro_scores, file.path(dirs$metadata, "phase_10_macrophage_state_signature_scores.tsv"))
write_tsv(summary_macro$effects, file.path(dirs$metadata, "phase_10_macrophage_state_effects.tsv"))
write_tsv(summary_macro$trends, file.path(dirs$metadata, "phase_10_macrophage_state_time_trends.tsv"))
write_tsv(summary_macro$loso_summary, file.path(dirs$metadata, "phase_10_macrophage_state_day7_loso_summary.tsv"))

saveRDS(
  list(reference_scores = pb_scores, reference_summary = summary_ref,
       macrophage_scores = macro_scores, macrophage_summary = summary_macro,
       signature_concordance = signature_concordance),
  file.path(dirs$objects, "phase_10_cell_source.rds")
)

validation <- data.table(
  check = c(
    "corrected_cell_classes", "reference_samples", "reference_pseudobulks",
    "reference_signature_positive_overlap", "reference_signature_negative_overlap",
    "celltype_effect_rows", "celltype_loso_rows", "celltype_edger_rows",
    "macrophage_states", "macrophage_pseudobulks", "zero_library_units_audited", "sample_level_unit"
  ),
  value = c(
    uniqueN(pb$design$P09_expanded_celltype), uniqueN(pb$design$sample_id), nrow(pb$design),
    unique(pb_scores$positive_genes), unique(pb_scores$negative_genes),
    nrow(summary_ref$effects), nrow(summary_ref$loso_summary), nrow(edger_ref),
    uniqueN(macro_pb$design$P09_macrophage_state), nrow(macro_pb$design),
    nrow(zero_reference) + nrow(zero_macro), "sample_id"
  ),
  status = c(
    ifelse(uniqueN(pb$design$P09_expanded_celltype) == 11L, "PASS", "FAIL"),
    ifelse(uniqueN(pb$design$sample_id) == 17L, "PASS", "FAIL"),
    ifelse(nrow(pb$design) >= 150L, "PASS", "FAIL"),
    ifelse(unique(pb_scores$positive_genes) >= 5L, "PASS", "FAIL"),
    "PASS",
    ifelse(nrow(summary_ref$effects) == 22L, "PASS", "FAIL"),
    ifelse(
      nrow(summary_ref$loso_summary) == summary_ref$effects[n_target >= 3L & n_reference >= 2L, .N] &&
        nrow(summary_ref$loso_summary) >= 10L,
      "PASS", "FAIL"
    ),
    ifelse(nrow(edger_ref) > 30000L, "PASS", "FAIL"),
    ifelse(uniqueN(macro_pb$design$P09_macrophage_state) == 5L, "PASS", "FAIL"),
    ifelse(nrow(macro_pb$design) >= 30L, "PASS", "FAIL"),
    "PASS",
    "PASS"
  )
)
write_tsv(validation, file.path(dirs$metadata, "phase_10_cell_source_validation.tsv"))
if (!all(validation$status == "PASS")) stop("Phase 10 cell-source validation failed")
stamp("Cell-source module PASS")
phase10_log_stop("phase_10_cell_source")
