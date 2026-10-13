args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: P03_run_meringue.R ROOT OUT")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- args[[2]]
source(file.path(root, "21_scop_expansion/P01/script/common.R"))
started_at <- now_utc()

.libPaths(unique(c(
  "/home/zhenhanli/R/x86_64-pc-linux-gnu-library/4.5",
  file.path(root, "Rlib"),
  .libPaths()
)))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(digest)
  library(scop)
})
if (pkg_version_or_na("scop") != "0.9.0") stop("P03 must use scop 0.9.0")
if (pkg_version_or_na("MERINGUE") != "1.0") stop("P03 MERINGUE 1.0 backend is unavailable")
if (!"RunMERINGUE" %in% getNamespaceExports("scop")) stop("RunMERINGUE is not exported")

p02 <- file.path(root, "21_scop_expansion/P02")
if (!file.exists(file.path(p02, "_SUCCESS"))) stop("P02 PASS is required")
p01 <- file.path(root, "21_scop_expansion/P01")
discovery <- read_tsv(file.path(p01, "tables/discovery_object_contract.tsv"))
primary <- discovery[discovery$biological_sample_id %in% c(
  "MI_D3_rep1", "MI_D3_rep2", "MI_D3_rep3", "MI_D7_rep1", "MI_D7_rep2", "MI_D7_rep3"
), , drop = FALSE]
primary <- primary[match(c("MI_D3_rep1", "MI_D7_rep1", "MI_D3_rep2", "MI_D3_rep3", "MI_D7_rep2", "MI_D7_rep3"), primary$biological_sample_id), , drop = FALSE]
if (nrow(primary) != 6L || anyNA(primary$biological_sample_id)) stop("Primary discovery contract is not six frozen tissues")

bank_features <- unique(scan(file.path(root, "08_recurrent_svg/BANKSY_discovery_features.txt"), what = character(), quiet = TRUE))
if (length(bank_features) != 931L) stop("Frozen BANKSY exclusion is not 931 features")
reference_markers <- read_tsv(file.path(root, "03_snrna_reference/reference_marker_audit_data.tsv"))
rctd_markers <- unique(reference_markers$features.plot)
myeloid_markers <- read_tsv(file.path(root, "03_snrna_reference/myeloid_refinement/myeloid_top15_cluster_markers.tsv"))
rctd_markers <- unique(c(rctd_markers, myeloid_markers$gene))
anchor <- read_tsv(file.path(root, "06_positive_control/MI_anchor_gene_coverage.tsv"))
endpoint_features <- unique(trimws(unlist(strsplit(paste(anchor$genes_present, collapse = ","), ",", fixed = TRUE))))

detection_rate <- 0.10
detection <- list()
all_features <- character()
for (i in seq_len(nrow(primary))) {
  id <- primary$biological_sample_id[i]
  message("P03 feature freeze: ", id)
  x <- readRDS(primary$object_file[i])
  mat <- LayerData(x[["Spatial"]], layer = "data")
  rate <- Matrix::rowSums(mat != 0) / ncol(mat)
  detection[[id]] <- setNames(as.numeric(rate), rownames(mat))
  all_features <- union(all_features, rownames(mat))
  rm(x, mat, rate); gc()
}
detection_matrix <- matrix(0, nrow = length(all_features), ncol = nrow(primary), dimnames = list(all_features, primary$biological_sample_id))
for (id in names(detection)) detection_matrix[names(detection[[id]]), id] <- detection[[id]]
detected_samples <- rowSums(detection_matrix >= detection_rate)
is_mito <- grepl("^mt-", all_features, ignore.case = TRUE)
is_ribo <- grepl("^Rp[sl][0-9]", all_features)
is_hb <- grepl("^Hb[ab][-_0-9]", all_features)
is_banksy <- all_features %in% bank_features
is_rctd_marker <- all_features %in% rctd_markers
is_endpoint <- all_features %in% endpoint_features
eligible <- detected_samples >= 3L & !is_mito & !is_ribo & !is_hb & !is_banksy & !is_rctd_marker & !is_endpoint
feature_table <- data.frame(
  feature = all_features,
  detected_primary_samples_at_10pct = detected_samples,
  max_detection_rate = apply(detection_matrix, 1L, max),
  mitochondrial = is_mito, ribosomal = is_ribo, hemoglobin = is_hb,
  frozen_BANKSY_feature = is_banksy, frozen_RCTD_marker = is_rctd_marker,
  frozen_corridor_endpoint_feature = is_endpoint,
  eligible = eligible,
  exclusion_reason = ifelse(eligible, "INCLUDED", ifelse(is_mito, "MITOCHONDRIAL", ifelse(is_ribo, "RIBOSOMAL", ifelse(is_hb, "HEMOGLOBIN", ifelse(is_banksy, "FROZEN_BANKSY_931", ifelse(is_rctd_marker, "FROZEN_RCTD_MARKER", ifelse(is_endpoint, "FROZEN_CORRIDOR_ENDPOINT", "DETECTED_IN_LT_3_OF_6"))))))),
  stringsAsFactors = FALSE
)
features_use <- feature_table$feature[feature_table$eligible]
if (length(features_use) < 2000L) stop("Frozen MERINGUE universe has fewer than 2000 eligible features")
write_tsv(feature_table, file.path(out, "tables/feature_universe.tsv"))
rm(detection, detection_matrix); gc()

meringue_formals <- paste(capture.output(dput(formals(scop::RunMERINGUE))), collapse = " ")
meringue_body_hash <- sha256_value(deparse(body(scop::RunMERINGUE)))
exports_hash <- sha256_value(sort(getNamespaceExports("scop")))
parameters <- data.frame(
  parameter = c(
    "scop_version", "MERINGUE_version", "mode", "assay", "layer", "coordinate_space", "feature_detection_rate",
    "feature_min_primary_samples", "excluded_BANKSY_features", "excluded_RCTD_markers", "excluded_endpoint_features",
    "nfeatures", "min_spots", "filterDist", "binary", "alternative", "nperm", "ncores", "seed", "backend",
    "module_match_similarity", "module_match_threshold", "frozen_membership_rule", "recurrence_stage", "recurrence_strong",
    "endpoint_remote_band", "endpoint_injury_band", "external_stage_map", "full_cross_correlation_storage"
  ),
  value = c(
    "0.9.0", "1.0", "autocorrelation;cross_correlation;modules", "Spatial", "data", "raw", "0.10",
    "3 of 6", as.character(length(bank_features)), as.character(length(rctd_markers)), as.character(length(endpoint_features)),
    "2000", "5", "NA", "TRUE", "greater", "0", "1", "11", "cpp",
    "gene-set Jaccard across different tissues within stage", "0.25", "gene present in >=2 component tissues", ">=2/3", "3/3",
    "position <=0.20", "position >=0.80", "Day3->external Day1; Day7->external Day7; Day14 descriptive", "full in Seurat tool; table retains top 1000 absolute cross-correlations per tissue"
  ), stringsAsFactors = FALSE
)
parameter_hash <- sha256_value(parameters[c("parameter", "value")])
parameters$parameter_hash <- parameter_hash
write_tsv(parameters, file.path(out, "parameters.tsv"))

call_provenance <- empty_df(c(
  "biological_sample_id", "call_started_at", "exported", "formals", "function_body_hash", "exports_hash",
  "scop_version", "backend_version", "input_object", "input_sha256", "parameter_hash", "input_class", "expected_output_tool"
))
write_tsv(call_provenance, file.path(out, "tables/call_provenance.tsv"))

native_plot <- function(x, id, top_features) {
  if (!length(top_features)) stop(id, ": no finite MERINGUE features for native plot")
  p <- scop::SpatialSpotPlot(
    srt = x, features = utils::head(top_features, 4L), assay = "Spatial", layer = "data",
    image = Images(x)[[1L]], overlay_image = TRUE, combine = TRUE, ncol = 2L, verbose = FALSE
  )
  grDevices::pdf(file.path(out, "figures", paste0("MERINGUE_", id, "_native_spatial_features.pdf")), width = 10, height = 8, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  TRUE
}

run_one <- function(i, pilot) {
  id <- primary$biological_sample_id[i]
  message("P03 RunMERINGUE: ", id, if (pilot) " [pilot]" else " [full]")
  x <- readRDS(primary$object_file[i])
  if (!inherits(x, "Seurat") || !"Spatial" %in% Assays(x) || length(Images(x)) != 1L) stop(id, ": input integrity gate failed")
  pre <- data.frame(
    biological_sample_id = id, call_started_at = now_utc(), exported = TRUE,
    formals = meringue_formals, function_body_hash = meringue_body_hash, exports_hash = exports_hash,
    scop_version = "0.9.0", backend_version = "1.0",
    input_object = primary$object_file[i], input_sha256 = primary$object_sha256[i], parameter_hash = parameter_hash,
    input_class = paste(class(x), collapse = ";"), expected_output_tool = "MERINGUE", stringsAsFactors = FALSE
  )
  call_provenance <<- rbind(call_provenance, pre)
  write_tsv(call_provenance, file.path(out, "tables/call_provenance.tsv"))
  capture.output(utils::sessionInfo(), file = file.path(out, "logs", paste0("precall_", id, "_sessionInfo.txt")))
  x <- scop::RunMERINGUE(
    srt = x, assay = "Spatial", layer = "data", image = Images(x)[[1L]],
    features = features_use, mode = c("autocorrelation", "cross_correlation", "modules"),
    nfeatures = 2000L, min_spots = 5L, filterDist = NA_real_, binary = TRUE,
    alternative = "greater", nperm = 0L, ncores = 1L, set_variable_features = FALSE,
    store_results = TRUE, verbose = TRUE, seed = 11L, coordinate_space = "raw", backend = "cpp"
  )
  tool <- x@tools$MERINGUE
  if (is.null(tool) || !all(c("autocorrelation", "cross_correlation", "modules", "coords", "weight", "parameters") %in% names(tool))) stop(id, ": pilot output audit failed")
  if (!identical(rownames(tool$coords), colnames(x))) stop(id, ": MERINGUE coordinate order mismatch")
  auto <- tool$autocorrelation
  top <- auto$feature[is.finite(auto$score)]
  top <- top[order(auto$score[is.finite(auto$score)], decreasing = TRUE)]
  plot_ok <- if (pilot) native_plot(x, id, top) else NA
  saveRDS(x, file.path(out, "objects", paste0(id, "_MERINGUE.rds")))
  audit <- data.frame(
    biological_sample_id = id, pilot = pilot, input_class = "Seurat", sparse_input = inherits(LayerData(x[["Spatial"]], layer = "data"), "sparseMatrix"),
    spatial_graph_scope = id, cross_slice_graph = FALSE, tool_slot = "MERINGUE", tool_complete = TRUE,
    n_autocorrelation = nrow(tool$autocorrelation), n_cross_correlation = nrow(tool$cross_correlation), n_module_features = nrow(tool$modules),
    native_plot_status = if (pilot) if (isTRUE(plot_ok)) "PASS" else "FAIL" else "NOT_REQUIRED_AFTER_PILOT",
    object_size_mb = as.numeric(object.size(x)) / 1024^2, stringsAsFactors = FALSE
  )
  rm(x); gc()
  audit
}

pilot_index <- match(c("MI_D3_rep1", "MI_D7_rep1"), primary$biological_sample_id)
pilot_audit <- do.call(rbind, lapply(pilot_index, function(i) run_one(i, TRUE)))
if (!all(pilot_audit$tool_complete) || !all(pilot_audit$native_plot_status == "PASS")) stop("P03 pilot gate failed")
remaining <- setdiff(seq_len(nrow(primary)), pilot_index)
full_audit <- do.call(rbind, lapply(remaining, function(i) run_one(i, FALSE)))
pilot_audit <- rbind(pilot_audit, full_audit)
write_tsv(pilot_audit, file.path(out, "tables/pilot_output_audit.tsv"))

auto_rows <- list(); module_rows <- list(); cross_rows <- list()
for (i in seq_len(nrow(primary))) {
  id <- primary$biological_sample_id[i]
  x <- readRDS(file.path(out, "objects", paste0(id, "_MERINGUE.rds")))
  tool <- x@tools$MERINGUE
  a <- tool$autocorrelation
  a$biological_sample_id <- id; a$timepoint <- primary$timepoint[i]; a$result_type <- "autocorrelation"
  auto_rows[[id]] <- a
  m <- tool$modules
  m$biological_sample_id <- id; m$timepoint <- primary$timepoint[i]
  module_rows[[id]] <- m
  cc <- tool$cross_correlation
  if (nrow(cc)) {
    cc <- cc[order(abs(cc$correlation), decreasing = TRUE), , drop = FALSE]
    cc <- utils::head(cc, 1000L)
    cc$biological_sample_id <- id; cc$timepoint <- primary$timepoint[i]; cc$result_type <- "cross_correlation_top1000"
    cross_rows[[id]] <- cc
  }
  rm(x); gc()
}
autocorrelation <- do.call(rbind, auto_rows)
modules_all <- do.call(rbind, module_rows)
cross_top <- if (length(cross_rows)) do.call(rbind, cross_rows) else empty_df(c("feature1", "feature2", "rank", "correlation", "p_value", "biological_sample_id", "timepoint", "result_type"))
autocorrelation$table_retention <- "complete"
cross_top$table_retention <- "top1000_abs_per_tissue; complete result in object"
per_sample <- rbind(
  data.frame(
    biological_sample_id = autocorrelation$biological_sample_id, timepoint = autocorrelation$timepoint,
    result_type = autocorrelation$result_type, feature1 = autocorrelation$feature, feature2 = NA_character_, module = NA_character_,
    statistic = autocorrelation$statistic, score = autocorrelation$score, correlation = NA_real_, p_value = autocorrelation$p_value,
    q_value = autocorrelation$q_value, rank = autocorrelation$rank, table_retention = autocorrelation$table_retention,
    spot_level_p_used_for_group_inference = FALSE, stringsAsFactors = FALSE
  ),
  data.frame(
    biological_sample_id = cross_top$biological_sample_id, timepoint = cross_top$timepoint,
    result_type = cross_top$result_type, feature1 = cross_top$feature1, feature2 = cross_top$feature2, module = NA_character_,
    statistic = NA_real_, score = NA_real_, correlation = cross_top$correlation, p_value = cross_top$p_value,
    q_value = NA_real_, rank = cross_top$rank, table_retention = cross_top$table_retention,
    spot_level_p_used_for_group_inference = FALSE, stringsAsFactors = FALSE
  ),
  data.frame(
    biological_sample_id = modules_all$biological_sample_id, timepoint = modules_all$timepoint,
    result_type = "module_assignment", feature1 = modules_all$feature, feature2 = NA_character_, module = as.character(modules_all$module),
    statistic = NA_real_, score = NA_real_, correlation = NA_real_, p_value = NA_real_, q_value = NA_real_, rank = modules_all$rank,
    table_retention = "complete", spot_level_p_used_for_group_inference = FALSE, stringsAsFactors = FALSE
  )
)
write_tsv(per_sample, file.path(out, "tables/per_sample_meringue_results.tsv"))

modules_all <- modules_all[!is.na(modules_all$module) & nzchar(as.character(modules_all$module)), , drop = FALSE]
modules_all$local_module_id <- paste(modules_all$biological_sample_id, modules_all$module, sep = "::")
module_sets <- split(modules_all$feature, modules_all$local_module_id)
module_meta <- unique(modules_all[c("local_module_id", "biological_sample_id", "timepoint", "module")])

parent <- setNames(module_meta$local_module_id, module_meta$local_module_id)
find_root <- function(x) { while (parent[[x]] != x) x <- parent[[x]]; x }
union_root <- function(a, b) { ra <- find_root(a); rb <- find_root(b); if (ra != rb) parent[[rb]] <<- ra }
match_rows <- list(); mr <- 0L
for (stage in unique(module_meta$timepoint)) {
  ids <- module_meta$local_module_id[module_meta$timepoint == stage]
  if (length(ids) < 2L) next
  cmb <- utils::combn(ids, 2L, simplify = FALSE)
  for (pair in cmb) {
    sa <- module_meta$biological_sample_id[match(pair[1], module_meta$local_module_id)]
    sb <- module_meta$biological_sample_id[match(pair[2], module_meta$local_module_id)]
    if (sa == sb) next
    ga <- unique(module_sets[[pair[1]]]); gb <- unique(module_sets[[pair[2]]])
    j <- length(intersect(ga, gb)) / length(union(ga, gb))
    matched <- is.finite(j) && j >= 0.25
    if (matched) union_root(pair[1], pair[2])
    mr <- mr + 1L
    match_rows[[mr]] <- data.frame(stage = stage, module_a = pair[1], module_b = pair[2], sample_a = sa, sample_b = sb,
                                   jaccard = j, threshold = 0.25, matched = matched, algorithm = "frozen_graph_components", stringsAsFactors = FALSE)
  }
}
matching <- if (length(match_rows)) do.call(rbind, match_rows) else empty_df(c("stage", "module_a", "module_b", "sample_a", "sample_b", "jaccard", "threshold", "matched", "algorithm"))
module_meta$component_root <- vapply(module_meta$local_module_id, find_root, character(1))
component_members <- split(module_meta$local_module_id, module_meta$component_root)
frozen_rows <- list(); fr <- 0L
for (component in component_members) {
  mm <- module_meta[match(component, module_meta$local_module_id), , drop = FALSE]
  samples <- unique(mm$biological_sample_id)
  if (length(samples) < 2L) next
  gene_sample <- unique(do.call(rbind, lapply(component, function(local_id) data.frame(
    gene = unique(module_sets[[local_id]]),
    sample = module_meta$biological_sample_id[match(local_id, module_meta$local_module_id)],
    stringsAsFactors = FALSE
  ))))
  gene_counts <- table(gene_sample$gene)
  genes <- names(gene_counts)[gene_counts >= 2L]
  if (length(genes) < 2L) next
  fr <- fr + 1L
  frozen_rows[[fr]] <- data.frame(
    stage = mm$timepoint[1], component_root = mm$component_root[1], local_modules = paste(component, collapse = ";"),
    discovery_samples = paste(sort(samples), collapse = ";"), n_discovery_samples = length(samples),
    recurrence_level = if (length(samples) == 3L) "strong_recurrent_3_of_3" else "stage_recurrent_2_of_3",
    genes = paste(sort(genes), collapse = ";"), n_genes = length(genes),
    membership_rule = "gene present in >=2 matched tissue modules", external_used_for_selection = FALSE,
    stringsAsFactors = FALSE
  )
}
frozen <- if (length(frozen_rows)) do.call(rbind, frozen_rows) else empty_df(c("stage", "component_root", "local_modules", "discovery_samples", "n_discovery_samples", "recurrence_level", "genes", "n_genes", "membership_rule", "external_used_for_selection"))
if (nrow(frozen)) {
  frozen <- frozen[order(frozen$stage, -frozen$n_discovery_samples, -frozen$n_genes, frozen$component_root), , drop = FALSE]
  frozen$frozen_module_id <- sprintf("%s_M%02d", gsub("Day", "D", frozen$stage), ave(seq_len(nrow(frozen)), frozen$stage, FUN = seq_along))
  frozen$membership_hash <- vapply(frozen$genes, sha256_value, character(1))
} else {
  frozen$frozen_module_id <- character(); frozen$membership_hash <- character()
}
if (nrow(matching)) {
  matching$component_a <- vapply(matching$module_a, find_root, character(1))
  matching$component_b <- vapply(matching$module_b, find_root, character(1))
}
write_tsv(matching, file.path(out, "tables/module_matching.tsv"))
write_tsv(frozen, file.path(out, "tables/frozen_recurrent_modules.tsv"))

position_for <- function(x, id, cohort, contract_row = NULL) {
  if (cohort == "GSE176092") {
    p <- read_tsv(contract_row$corridor_table)
    if (!identical(p$barcode, colnames(x))) stop(id, ": external corridor order changed")
    return(p$RZ_injury_position)
  }
  meta <- x[[]]
  col <- intersect(c("RZ_injury_position", "FrozenTransfer_position"), colnames(meta))
  if (!length(col)) return(rep(NA_real_, ncol(x)))
  meta[[col[1]]]
}

module_effects_for <- function(contract, cohort, role) {
  if (!nrow(frozen)) return(empty_df(c("biological_sample_id", "cohort", "condition", "timepoint", "analysis_role", "projection_role", "frozen_module_id", "module_stage", "n_genes", "n_genes_present", "mean_module_score", "corridor_spearman_rho", "remote_mean", "injury_mean", "injury_minus_remote", "median_member_moran_statistic", "inference_unit")))
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(contract))) {
    id <- contract$biological_sample_id[i]
    message("P03 frozen module scoring: ", cohort, " / ", id)
    object_path <- if (role == "primary_discovery") file.path(out, "objects", paste0(id, "_MERINGUE.rds")) else contract$object_file[i]
    x <- readRDS(object_path)
    mat <- LayerData(x[["Spatial"]], layer = "data")
    pos <- position_for(x, id, cohort, contract[i, , drop = FALSE])
    auto <- if ("MERINGUE" %in% names(x@tools)) x@tools$MERINGUE$autocorrelation else NULL
    for (j in seq_len(nrow(frozen))) {
      genes <- strsplit(frozen$genes[j], ";", fixed = TRUE)[[1]]
      present <- intersect(genes, rownames(mat))
      score <- if (length(present)) Matrix::colMeans(mat[present, , drop = FALSE]) else rep(NA_real_, ncol(x))
      ok <- is.finite(score) & is.finite(pos)
      rho <- if (sum(ok) >= 10L && stats::sd(score[ok]) > 0 && stats::sd(pos[ok]) > 0) suppressWarnings(stats::cor(score[ok], pos[ok], method = "spearman")) else NA_real_
      remote <- mean(score[is.finite(pos) & pos <= 0.20], na.rm = TRUE)
      injury <- mean(score[is.finite(pos) & pos >= 0.80], na.rm = TRUE)
      if (!is.finite(remote)) remote <- NA_real_
      if (!is.finite(injury)) injury <- NA_real_
      moran <- if (!is.null(auto)) stats::median(auto$statistic[match(present, auto$feature)], na.rm = TRUE) else NA_real_
      if (!is.finite(moran)) moran <- NA_real_
      k <- k + 1L
      rows[[k]] <- data.frame(
        biological_sample_id = id, cohort = cohort, condition = contract$condition[i], timepoint = contract$timepoint[i],
        analysis_role = contract$analysis_role[i], projection_role = role, frozen_module_id = frozen$frozen_module_id[j], module_stage = frozen$stage[j],
        n_genes = length(genes), n_genes_present = length(present), mean_module_score = mean(score, na.rm = TRUE),
        corridor_spearman_rho = rho, remote_mean = remote, injury_mean = injury, injury_minus_remote = injury - remote,
        median_member_moran_statistic = moran, inference_unit = "physical_tissue", stringsAsFactors = FALSE
      )
    }
    rm(x, mat); gc()
  }
  do.call(rbind, rows)
}

primary_contract <- primary[match(c("MI_D3_rep1", "MI_D3_rep2", "MI_D3_rep3", "MI_D7_rep1", "MI_D7_rep2", "MI_D7_rep3"), primary$biological_sample_id), , drop = FALSE]
context_contract <- discovery[!discovery$biological_sample_id %in% primary_contract$biological_sample_id, , drop = FALSE]
primary_effects <- module_effects_for(primary_contract, "GSE214611", "primary_discovery")
context_effects <- module_effects_for(context_contract, "GSE214611", "contextual_descriptive")
module_effects <- rbind(primary_effects, context_effects)
write_tsv(module_effects, file.path(out, "tables/module_sample_effects.tsv"))

# External biological objects are first opened here, after discovery modules and memberships are frozen on disk.
external_contract <- read_tsv(file.path(p01, "tables/external_object_contract.tsv"))
external_effects <- module_effects_for(external_contract, "GSE176092", "external_frozen_transfer")
external_effects$external_used_for_selection <- FALSE
write_tsv(external_effects, file.path(out, "tables/external_module_transfer.tsv"))

summary_rows <- list()
if (nrow(frozen)) for (j in seq_len(nrow(frozen))) {
  fm <- frozen$frozen_module_id[j]
  stage <- frozen$stage[j]
  d <- primary_effects[primary_effects$frozen_module_id == fm & primary_effects$timepoint == stage, , drop = FALSE]
  expected_sign <- sign(stats::median(d$corridor_spearman_rho, na.rm = TRUE))
  concordant <- sum(sign(d$corridor_spearman_rho) == expected_sign, na.rm = TRUE)
  external_stage <- if (stage == "Day3") "Day1" else "Day7"
  e <- external_effects[external_effects$frozen_module_id == fm & external_effects$timepoint == external_stage, , drop = FALSE]
  ext_concordant <- sum(sign(e$corridor_spearman_rho) == expected_sign, na.rm = TRUE)
  summary_rows[[j]] <- data.frame(
    frozen_module_id = fm, stage = stage, n_discovery_tissues = nrow(d), discovery_direction = expected_sign,
    discovery_direction_concordance = paste0(concordant, "/", nrow(d)), median_discovery_rho = stats::median(d$corridor_spearman_rho, na.rm = TRUE),
    external_stage = external_stage, n_external_tissues = nrow(e), external_direction_concordance = paste0(ext_concordant, "/", nrow(e)),
    median_external_rho = stats::median(e$corridor_spearman_rho, na.rm = TRUE),
    main_figure_candidate = frozen$n_discovery_samples[j] == 3L && concordant == 3L && nrow(e) == 3L && ext_concordant >= 2L,
    biological_verdict = if (frozen$n_discovery_samples[j] >= 2L) "RECURRENT_SPATIAL_MODULE_ASSOCIATION" else "DESCRIPTIVE_ONLY",
    inference_unit = "physical_tissue", spot_level_p_used_for_group_inference = FALSE, stringsAsFactors = FALSE
  )
}
stage_summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else empty_df(c("frozen_module_id", "stage", "n_discovery_tissues", "discovery_direction", "discovery_direction_concordance", "median_discovery_rho", "external_stage", "n_external_tissues", "external_direction_concordance", "median_external_rho", "main_figure_candidate", "biological_verdict", "inference_unit", "spot_level_p_used_for_group_inference"))
write_tsv(stage_summary, file.path(out, "tables/module_stage_summary.tsv"))

write_session_info(file.path(out, "sessionInfo.txt"))
input_hash <- sha256_value(sort(primary$object_sha256))
status <- if (nrow(frozen)) "PASS" else "NEGATIVE_COMPLETE"
verdict <- if (nrow(frozen)) paste0("FROZEN_RECURRENT_MODULES_N=", nrow(frozen)) else "NO_STAGE_RECURRENT_MERINGUE_MODULE"
write_phase_status(
  out = out, phase = "P03_MERINGUE", status = status, started_at = started_at,
  input_hash = input_hash, parameter_hash = parameter_hash,
  n_expected = 6L, n_completed = 6L, producer_functions = "scop::RunMERINGUE",
  biological_verdict = verdict, blocker = "", output_root = out
)
file.create(file.path(out, "_SUCCESS"))
message("P03 ", status, ": frozen recurrent modules = ", nrow(frozen))
