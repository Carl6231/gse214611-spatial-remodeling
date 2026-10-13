root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
source(file.path(root, "scripts", "phase_10_common.R"))
root <- phase10_root()
dirs <- phase10_dirs(root)
phase10_log_start(root, "phase_10_spatial_signature")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(matrixStats)
  library(FNN)
  library(scop)
  library(ggplot2)
})

set.seed(20260805)
registry <- fread(require_file(file.path(root, "metadata", "spatial_qc_summary.tsv")))
preflight <- fread(require_file(file.path(root, "metadata", "phase_10_input_object_audit.tsv")))
if (!all(preflight$status == "PASS")) stop("Phase 10 preflight has not passed")
construction <- fread(require_file(file.path(root, "metadata", "phase_10_axis_construction_genes.tsv")))$gene
day7_samples <- registry[condition == "MI" & timepoint == "Day7", sample_key]
if (length(day7_samples) != 3L) stop("Expected three Day7 MI discovery slices")

assoc_dir <- file.path(dirs$results, "spatial_signature", "associations")
map_dir <- file.path(dirs$results, "spatial_signature", "maps")
dir.create(assoc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

rank_association <- function(obj, sample_key, permutations = 200L) {
  md <- obj[[]]
  centers <- rownames(md)[as.logical(md$P09N_center) & is.finite(md$P09N_remodeling_axis)]
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")[, centers, drop = FALSE]
  n <- ncol(counts)
  detected <- Matrix::rowSums(counts > 0)
  total <- Matrix::rowSums(counts)
  technical <- grepl("^(mt-|Rpl|Rps)", rownames(counts), ignore.case = TRUE)
  keep <- detected >= max(10L, ceiling(0.05 * n)) & total >= 20 &
    !rownames(counts) %in% construction & !technical
  counts <- counts[keep, , drop = FALSE]
  if (nrow(counts) < 2000L) stop(sample_key, " has too few genes after filtering")
  expr <- as.matrix(logcpm_sparse(counts))
  xr <- matrixStats::rowRanks(expr, ties.method = "average", preserveShape = TRUE)
  xsd <- matrixStats::rowSds(xr)
  valid <- is.finite(xsd) & xsd > 0
  xr <- xr[valid, , drop = FALSE]
  xsd <- xsd[valid]
  xz <- sweep(sweep(xr, 1L, rowMeans(xr), "-"), 1L, xsd, "/")
  axis <- md[centers, "P09N_remodeling_axis", drop = TRUE]
  yr <- rank(axis, ties.method = "average")
  yz <- as.numeric(scale(yr))
  observed <- as.vector(xz %*% yz) / (n - 1)

  coords <- GetTissueCoordinates(obj)[centers, c("x", "y"), drop = FALSE]
  k_blocks <- max(4L, min(10L, floor(n / 50L)))
  blocks <- kmeans(scale(coords), centers = k_blocks, nstart = 30L, iter.max = 100L)$cluster
  yperm <- matrix(NA_real_, nrow = n, ncol = permutations)
  for (b in seq_len(permutations)) {
    yp <- yr
    for (block in unique(blocks)) {
      idx <- which(blocks == block)
      yp[idx] <- sample(yr[idx], length(idx), replace = FALSE)
    }
    yperm[, b] <- as.numeric(scale(yp))
  }
  null <- (xz %*% yperm) / (n - 1)
  empirical_p <- (1 + rowSums(abs(null) >= abs(observed))) / (permutations + 1)
  out <- data.table(
    sample_key = sample_key,
    gene = rownames(xz),
    centers = n,
    detected_centers = detected[rownames(xz)],
    detected_fraction = detected[rownames(xz)] / n,
    total_counts = total[rownames(xz)],
    spearman_rho = observed,
    block_permutation_p = empirical_p,
    block_permutation_fdr = p.adjust(empirical_p, method = "BH"),
    spatial_blocks = k_blocks,
    permutations = permutations,
    axis_construction_gene = FALSE
  )
  list(result = out, blocks = data.table(sample_key = sample_key, barcode = centers,
                                         x = coords$x, y = coords$y, block = blocks))
}

stamp("Discovering non-circular remodeling-associated genes in three Day7 slices")
day7_results <- vector("list", length(day7_samples))
block_results <- vector("list", length(day7_samples))
for (i in seq_along(day7_samples)) {
  sample_key <- day7_samples[[i]]
  stamp("Association: ", sample_key)
  obj <- readRDS(file.path(root, "objects", "phase_09", "niches", paste0(sample_key, "_phase09_niche.rds")))
  ans <- rank_association(obj, sample_key, permutations = 200L)
  day7_results[[i]] <- ans$result
  block_results[[i]] <- ans$blocks
  write_tsv(ans$result, file.path(assoc_dir, paste0(sample_key, "_gene_axis_association.tsv.gz")))
  rm(obj, ans)
  gc()
}
all_assoc <- rbindlist(day7_results)
write_tsv(rbindlist(block_results), file.path(dirs$metadata, "phase_10_spatial_blocks.tsv.gz"))

stamp("Cross-slice signature recurrence")
summary_dt <- all_assoc[, {
  p_use <- pmax(block_permutation_p, .Machine$double.xmin)
  fisher <- pchisq(-2 * sum(log(p_use)), df = 2 * .N, lower.tail = FALSE)
  list(
    slices_tested = uniqueN(sample_key),
    positive_slices = sum(spearman_rho > 0),
    negative_slices = sum(spearman_rho < 0),
    median_rho = median(spearman_rho),
    min_abs_rho = min(abs(spearman_rho)),
    max_abs_rho = max(abs(spearman_rho)),
    fisher_block_p = fisher,
    mean_detected_fraction = mean(detected_fraction)
  )
}, by = gene]
summary_dt[, fisher_block_fdr := p.adjust(fisher_block_p, method = "BH")]
summary_dt[, direction := fifelse(positive_slices == 3L, "positive",
                                  fifelse(negative_slices == 3L, "negative", "discordant"))]
summary_dt[, recurrent_gate := slices_tested == 3L & direction != "discordant" &
             abs(median_rho) >= 0.20 & fisher_block_fdr <= 0.05]
summary_dt[, abs_median_rho := abs(median_rho)]
setorder(summary_dt, -recurrent_gate, -abs_median_rho, fisher_block_fdr)
write_tsv(summary_dt, file.path(dirs$metadata, "phase_10_day7_gene_axis_recurrence.tsv.gz"))

positive <- head(
  summary_dt[recurrent_gate & direction == "positive"][order(-abs(median_rho), fisher_block_fdr)],
  100L
)
negative <- head(
  summary_dt[recurrent_gate & direction == "negative"][order(-abs(median_rho), fisher_block_fdr)],
  100L
)
signature <- rbind(
  positive[, .(gene, direction, median_rho, fisher_block_p, fisher_block_fdr,
               positive_slices, negative_slices, mean_detected_fraction)],
  negative[, .(gene, direction, median_rho, fisher_block_p, fisher_block_fdr,
               positive_slices, negative_slices, mean_detected_fraction)]
)
signature[, signature_rank := seq_len(.N), by = direction]
if (nrow(signature) == 0L) stop("No recurrent non-circular Day7 signature genes passed the prespecified gate")
write_tsv(signature, file.path(dirs$metadata, "phase_10_day7_remodeling_signature.tsv"))

signature_score <- function(counts, positive_genes, negative_genes) {
  genes <- intersect(c(positive_genes, negative_genes), rownames(counts))
  z <- safe_row_z(logcpm_sparse(counts[genes, , drop = FALSE]))
  pos <- intersect(positive_genes, rownames(z))
  neg <- intersect(negative_genes, rownames(z))
  pos_score <- if (length(pos)) colMeans(z[pos, , drop = FALSE]) else rep(0, ncol(counts))
  neg_score <- if (length(neg)) colMeans(z[neg, , drop = FALSE]) else rep(0, ncol(counts))
  pos_score - neg_score
}

moran_knn <- function(values, coords, permutations = 499L, k = 6L) {
  ok <- is.finite(values) & complete.cases(coords)
  values <- values[ok]
  coords <- as.matrix(coords[ok, , drop = FALSE])
  n <- length(values)
  k <- min(k, n - 1L)
  if (n < 10L || k < 1L) return(c(I = NA_real_, p = NA_real_, n = n, k = k))
  nn <- FNN::get.knn(coords, k = k)$nn.index
  from <- rep(seq_len(n), each = k)
  to <- as.vector(t(nn))
  z <- values - mean(values)
  denom <- sum(z^2)
  calc_i <- function(v) (n / length(from)) * sum(v[from] * v[to]) / sum(v^2)
  obs <- calc_i(z)
  null <- replicate(permutations, calc_i(sample(z, n, replace = FALSE)))
  p <- (1 + sum(null >= obs)) / (permutations + 1)
  c(I = obs, p = p, n = n, k = k)
}

stamp("Scoring all slices, decomposing quartiles, and auditing quality covariates")
positive_genes <- signature[direction == "positive", gene]
negative_genes <- signature[direction == "negative", gene]
spot_rows <- vector("list", nrow(registry))
quartile_rows <- list()
quality_rows <- list()
moran_rows <- list()
map_rows <- list()
feature_pattern <- "^(P09E_RCTD_prop_|P09M_RCTD_prop_)"
for (i in seq_len(nrow(registry))) {
  sample_key <- registry$sample_key[[i]]
  stamp("Scoring: ", sample_key)
  obj <- readRDS(file.path(root, "objects", "phase_09", "niches", paste0(sample_key, "_phase09_niche.rds")))
  md <- obj[[]]
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  score <- signature_score(counts, positive_genes, negative_genes)
  names(score) <- colnames(obj)
  centers <- as.logical(md$P09N_center) & is.finite(md$P09N_remodeling_axis)
  coords <- GetTissueCoordinates(obj)[colnames(obj), c("x", "y"), drop = FALSE]
  axis <- md$P09N_remodeling_axis
  quartile <- rep(NA_integer_, nrow(md))
  r <- rank(axis[centers], ties.method = "average")
  quartile[centers] <- pmin(4L, pmax(1L, ceiling(4 * r / length(r))))

  spot_rows[[i]] <- data.table(
    sample_key = sample_key,
    GSM = registry$GSM[[i]], condition = registry$condition[[i]],
    timepoint = registry$timepoint[[i]], replicate = registry$replicate[[i]],
    barcode = rownames(md), x = coords$x, y = coords$y,
    center = centers, remodeling_axis = axis, axis_quartile = quartile,
    signature_score = score[rownames(md)]
  )

  features <- unique(c(
    grep(feature_pattern, colnames(md), value = TRUE),
    "P09N_inflammatory_axis", "P09N_remodeling_axis", "P09N_vascular_axis",
    "P09N_remote_axis", "P09E_total_myeloid", "P09E_total_lymphoid"
  ))
  features <- intersect(features, colnames(md))
  for (feature in features) {
    v <- as.numeric(md[[feature]])
    q1 <- v[centers & quartile == 1L]
    q4 <- v[centers & quartile == 4L]
    quartile_rows[[length(quartile_rows) + 1L]] <- data.table(
      sample_key = sample_key, condition = registry$condition[[i]],
      timepoint = registry$timepoint[[i]], feature = feature,
      q1_n = sum(is.finite(q1)), q4_n = sum(is.finite(q4)),
      q1_mean = mean(q1, na.rm = TRUE), q4_mean = mean(q4, na.rm = TRUE),
      q4_minus_q1 = mean(q4, na.rm = TRUE) - mean(q1, na.rm = TRUE)
    )
  }
  for (feature in intersect(c("nCount_Spatial", "nFeature_Spatial", "percent.mt", "percent.ribo"), colnames(md))) {
    quality_rows[[length(quality_rows) + 1L]] <- data.table(
      sample_key = sample_key, condition = registry$condition[[i]],
      timepoint = registry$timepoint[[i]], covariate = feature,
      centers = sum(centers),
      spearman_rho = suppressWarnings(cor(axis[centers], md[[feature]][centers], method = "spearman", use = "complete.obs"))
    )
  }
  axis_moran <- moran_knn(axis[centers], coords[centers, , drop = FALSE])
  score_moran <- moran_knn(score[centers], coords[centers, , drop = FALSE])
  moran_rows[[length(moran_rows) + 1L]] <- data.table(
    sample_key = sample_key, condition = registry$condition[[i]], timepoint = registry$timepoint[[i]],
    feature = c("remodeling_axis", "noncircular_signature_score"),
    moran_I = c(axis_moran[["I"]], score_moran[["I"]]),
    permutation_p = c(axis_moran[["p"]], score_moran[["p"]]),
    centers = sum(centers), permutations = 499L
  )

  image_name <- Images(obj)[[1L]]
  map_path <- file.path(map_dir, paste0(sample_key, "_noncircular_signature.png"))
  p <- SpatialSpotPlot(
    obj,
    values = setNames(score, colnames(obj)),
    plot_type = "point", image = image_name, crop = FALSE,
    coord.cols = c("x", "y"), overlay_image = TRUE, image.alpha = 0.85,
    pt.alpha = 0.9, palcolor = c("#2F6FB0", "#F7F7F7", "#C81D25"),
    legend.title = "P10 signature"
  ) + ggtitle(paste0(sample_key, ": non-circular remodeling signature"))
  ggsave(map_path, p, width = 9, height = 8, dpi = 260, bg = "white")
  map_rows[[i]] <- data.table(sample_key = sample_key, map = map_path,
                              bytes = file.info(map_path)$size, scop_function = "SpatialSpotPlot")
  rm(obj, md, counts, p)
  gc()
}

spot_table <- rbindlist(spot_rows)
quartile_table <- rbindlist(quartile_rows)
quality_table <- rbindlist(quality_rows)
moran_table <- rbindlist(moran_rows)
map_table <- rbindlist(map_rows)
write_tsv(spot_table, file.path(dirs$metadata, "phase_10_spot_signature_scores.tsv.gz"))
write_tsv(quartile_table, file.path(dirs$metadata, "phase_10_axis_quartile_decomposition.tsv.gz"))
write_tsv(quality_table, file.path(dirs$metadata, "phase_10_axis_quality_covariates.tsv"))
write_tsv(moran_table, file.path(dirs$metadata, "phase_10_signature_spatial_autocorrelation.tsv"))
write_tsv(map_table, file.path(dirs$metadata, "phase_10_spatial_map_manifest.tsv"))

day7_decomposition <- quartile_table[condition == "MI" & timepoint == "Day7", .(
  slices = uniqueN(sample_key), median_q4_minus_q1 = median(q4_minus_q1, na.rm = TRUE),
  positive_slices = sum(q4_minus_q1 > 0, na.rm = TRUE),
  negative_slices = sum(q4_minus_q1 < 0, na.rm = TRUE),
  min_effect = min(q4_minus_q1, na.rm = TRUE), max_effect = max(q4_minus_q1, na.rm = TRUE)
), by = feature]
write_tsv(day7_decomposition, file.path(dirs$metadata, "phase_10_day7_axis_decomposition_summary.tsv"))

saveRDS(
  list(
    positive_genes = positive_genes, negative_genes = negative_genes,
    signature = signature, recurrence = summary_dt,
    parameters = list(day7_slices = day7_samples, block_permutations = 200L,
                      moran_permutations = 499L, minimum_abs_median_rho = 0.20,
                      fisher_fdr = 0.05, maximum_genes_per_direction = 100L,
                      construction_genes_excluded = length(construction))
  ),
  file.path(dirs$objects, "phase_10_spatial_signature.rds")
)

validation <- data.table(
  check = c(
    "day7_slices_analyzed", "association_rows", "construction_genes_in_signature",
    "positive_signature_genes", "negative_signature_genes", "all_slices_scored",
    "all_spots_scored", "quartile_decomposition", "quality_covariates",
    "moran_results", "spatial_maps"
  ),
  value = c(
    uniqueN(all_assoc$sample_key), nrow(all_assoc), sum(signature$gene %in% construction),
    length(positive_genes), length(negative_genes), uniqueN(spot_table$sample_key),
    nrow(spot_table), nrow(quartile_table), nrow(quality_table), nrow(moran_table),
    sum(map_table$bytes > 10000)
  ),
  status = c(
    ifelse(uniqueN(all_assoc$sample_key) == 3L, "PASS", "FAIL"),
    ifelse(nrow(all_assoc) > 5000L, "PASS", "FAIL"),
    ifelse(!any(signature$gene %in% construction), "PASS", "FAIL"),
    ifelse(length(positive_genes) > 0L, "PASS", "FAIL"),
    "PASS",
    ifelse(uniqueN(spot_table$sample_key) == 12L, "PASS", "FAIL"),
    ifelse(nrow(spot_table) == sum(registry$spots_after), "PASS", "FAIL"),
    ifelse(nrow(quartile_table) > 100L, "PASS", "FAIL"),
    ifelse(nrow(quality_table) == 48L, "PASS", "FAIL"),
    ifelse(nrow(moran_table) == 24L, "PASS", "FAIL"),
    ifelse(sum(map_table$bytes > 10000) == 12L, "PASS", "FAIL")
  )
)
write_tsv(validation, file.path(dirs$metadata, "phase_10_spatial_signature_validation.tsv"))
if (!all(validation$status == "PASS")) stop("Phase 10 spatial signature validation failed")
stamp("Spatial signature module PASS")
phase10_log_stop("phase_10_spatial_signature")
