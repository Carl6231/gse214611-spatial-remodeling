root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
source(file.path(root, "scripts", "phase_10_common.R"))
root <- phase10_root()
dirs <- phase10_dirs(root)
phase10_log_start(root, "phase_10_external_projection")

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(matrixStats)
})

sig <- readRDS(require_file(file.path(dirs$objects, "phase_10_spatial_signature.rds")))
recurrence <- as.data.table(sig$recurrence)
positive <- sig$positive_genes
negative <- sig$negative_genes
raw_dir <- file.path(root, "results", "phase_09", "external_validation", "GSE280373", "raw")
expr_path <- require_file(file.path(raw_dir, "GSE280373_expdata_all_timepoints_LAD_Cryo_Sham.rds"))
meta_path <- require_file(file.path(raw_dir, "GSE280373_mdata_cells_timepoints_annotations_final.rds"))

stamp("Loading independent GSE280373 processed counts and annotations")
expr <- readRDS(expr_path)
meta <- readRDS(meta_path)
if (!inherits(expr, "sparseMatrix") || !identical(colnames(expr), rownames(meta))) stop("External matrix and metadata are not aligned")
if (!all(c("cells_timepoints", "cluster_ann_all") %in% colnames(meta))) stop("External annotation fields missing")
if (any(expr@x < 0) || any(abs(expr@x - round(expr@x)) > 1e-8)) stop("External matrix is not nonnegative integer counts")

positive_use <- intersect(positive, rownames(expr))
negative_use <- intersect(negative, rownames(expr))
genes <- unique(c(positive_use, negative_use))
if (length(genes) < 10L || length(positive_use) < 5L) stop("Too few Phase 10 signature genes overlap GSE280373")

stamp("Computing descriptive cell-level projection scores")
lognorm <- as.matrix(logcpm_sparse(expr[genes, , drop = FALSE]))
z <- safe_row_z(lognorm)
pos_score <- colMeans(z[positive_use, , drop = FALSE])
neg_score <- if (length(negative_use)) colMeans(z[negative_use, , drop = FALSE]) else rep(0, ncol(expr))
score <- pos_score - neg_score
external_group <- as.character(meta$cells_timepoints)
annotation <- as.character(meta$cluster_ann_all)

group_annotation <- data.table(
  external_group = external_group, annotation = annotation,
  positive_score = pos_score, negative_score = neg_score,
  remodeling_signature_score = score
)[, .(
  cells = .N,
  mean_score = mean(remodeling_signature_score),
  median_score = median(remodeling_signature_score),
  score_sd = sd(remodeling_signature_score),
  q25 = quantile(remodeling_signature_score, 0.25),
  q75 = quantile(remodeling_signature_score, 0.75),
  mean_positive = mean(positive_score), mean_negative = mean(negative_score)
), by = .(external_group, annotation)]
setorder(group_annotation, external_group, -mean_score)
write_tsv(group_annotation, file.path(dirs$metadata, "phase_10_external_annotation_projection.tsv.gz"))

lineage <- fifelse(grepl("^(Macro_|mo_Macro_)", annotation), "Macrophage",
                   fifelse(grepl("^(FB|Fbmyo)", annotation), "Fibroblast", "Other"))
lineage_summary <- data.table(
  external_group = external_group, lineage = lineage,
  remodeling_signature_score = score
)[, .(
  cells = .N, mean_score = mean(remodeling_signature_score),
  median_score = median(remodeling_signature_score),
  q25 = quantile(remodeling_signature_score, 0.25),
  q75 = quantile(remodeling_signature_score, 0.75)
), by = .(external_group, lineage)]
write_tsv(lineage_summary, file.path(dirs$metadata, "phase_10_external_lineage_projection.tsv"))

stamp("Computing gene-level descriptive external directions")
cell_library <- Matrix::colSums(expr)
split_cells <- split(colnames(expr), external_group)
gene_group <- rbindlist(lapply(names(split_cells), function(group) {
  cells <- split_cells[[group]]
  sums <- Matrix::rowSums(expr[genes, cells, drop = FALSE])
  data.table(
    gene = genes, external_group = group, cells = length(cells), counts = sums,
    cpm = 1e6 * sums / sum(cell_library[cells]),
    detected_fraction = Matrix::rowMeans(expr[genes, cells, drop = FALSE] > 0)
  )
}))
wide <- dcast(gene_group, gene ~ external_group, value.var = "cpm", fill = 0)
required_groups <- c("LAD_Day3", "LAD_Day7", "Cryo_Day3", "Cryo_Day7")
if (!all(required_groups %in% colnames(wide))) stop("Required GSE280373 groups missing")
wide[, `:=`(
  external_LAD_D7_vs_D3 = log2((LAD_Day7 + 0.1) / (LAD_Day3 + 0.1)),
  external_Cryo_D7_vs_D3 = log2((Cryo_Day7 + 0.1) / (Cryo_Day3 + 0.1)),
  external_LAD_vs_Cryo_D7 = log2((LAD_Day7 + 0.1) / (Cryo_Day7 + 0.1))
)]
wide <- merge(
  wide,
  recurrence[, .(gene, internal_median_rho = median_rho, internal_direction = direction,
                 internal_fisher_fdr = fisher_block_fdr)],
  by = "gene", all.x = TRUE
)
wide[, signature_direction := fifelse(gene %in% positive_use, "positive", "negative")]
write_tsv(gene_group, file.path(dirs$metadata, "phase_10_external_signature_expression.tsv.gz"))
write_tsv(wide, file.path(dirs$metadata, "phase_10_external_signature_directions.tsv"))

concordance <- rbindlist(lapply(
  c("external_LAD_D7_vs_D3", "external_Cryo_D7_vs_D3", "external_LAD_vs_Cryo_D7"),
  function(comparison) {
    rho <- suppressWarnings(cor(wide$internal_median_rho, wide[[comparison]], method = "spearman", use = "complete.obs"))
    data.table(
      comparison = comparison, genes = sum(complete.cases(wide$internal_median_rho, wide[[comparison]])),
      spearman_rho = rho,
      direction_concordance = mean(sign(wide$internal_median_rho) == sign(wide[[comparison]]), na.rm = TRUE),
      inference_scope = "gene-level descriptive; processed metadata lack animal identifiers"
    )
  }
))
write_tsv(concordance, file.path(dirs$metadata, "phase_10_external_gene_concordance.tsv"))

saveRDS(
  list(group_annotation = group_annotation, lineage_summary = lineage_summary,
       gene_directions = wide, concordance = concordance,
       positive_genes = positive_use, negative_genes = negative_use,
       inference_scope = "descriptive_only_no_animal_level_replicates"),
  file.path(dirs$objects, "phase_10_external_projection.rds")
)

validation <- data.table(
  check = c(
    "external_cells", "aligned_cells", "integer_counts", "signature_overlap",
    "positive_overlap", "external_groups", "annotation_groups", "lineages",
    "gene_direction_rows", "descriptive_scope"
  ),
  value = c(
    ncol(expr), identical(colnames(expr), rownames(meta)), TRUE, length(genes),
    length(positive_use), uniqueN(external_group), uniqueN(annotation), uniqueN(lineage),
    nrow(wide), "no_animal_level_inference"
  ),
  status = c(
    ifelse(ncol(expr) == nrow(meta), "PASS", "FAIL"), "PASS", "PASS",
    ifelse(length(genes) >= 10L, "PASS", "FAIL"),
    ifelse(length(positive_use) >= 5L, "PASS", "FAIL"),
    ifelse(all(required_groups %in% external_group), "PASS", "FAIL"),
    ifelse(uniqueN(annotation) >= 10L, "PASS", "FAIL"),
    ifelse(all(c("Macrophage", "Fibroblast", "Other") %in% lineage), "PASS", "FAIL"),
    ifelse(nrow(wide) == length(genes), "PASS", "FAIL"), "PASS"
  )
)
write_tsv(validation, file.path(dirs$metadata, "phase_10_external_projection_validation.tsv"))
if (!all(validation$status == "PASS")) stop("Phase 10 external projection validation failed")
stamp("External projection module PASS")
phase10_log_stop("phase_10_external_projection")
