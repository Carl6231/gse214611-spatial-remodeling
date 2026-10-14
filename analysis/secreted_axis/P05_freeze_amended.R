#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: P05_freeze_amended.R ROOT OUT")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- normalizePath(args[[2]], mustWork = TRUE)
expected_out <- file.path(root, "21_scop_expansion/P05_SECRETED_AXIS_GATE_AMENDED")
if (!identical(out, expected_out)) stop("Output is outside the authorized amended P05 directory")
source(file.path(root, "21_scop_expansion/P01/script/common.R"))

write_failed <- function(e, started_at) {
  dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
  writeLines(conditionMessage(e), file.path(out, "logs/error_message.txt"))
  try(write_session_info(file.path(out, "sessionInfo.txt")), silent = TRUE)
  s <- phase_status(
    phase = "P05_SECRETED_AXIS_GATE_AMENDED", status = "FAILED",
    started_at = started_at, input_hash = NA_character_, parameter_hash = NA_character_,
    n_expected = 6L, n_completed = 0L, producer_functions = "NONE_REUSE_ONLY",
    biological_verdict = "FAILED_INPUT_MAPPING_OR_OUTPUT_INTEGRITY",
    blocker = conditionMessage(e), output_root = out
  )
  try(write_tsv(s, file.path(out, "STATUS.tsv")), silent = TRUE)
}

main <- function() {
  started_at <- now_utc()
  threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  if (!all(Sys.getenv(threads) == "1")) stop("Thread environment is not fixed to 1")
  suppressPackageStartupMessages({
    library(SeuratObject)
    library(CellChat)
    library(SecAct)
    library(digest)
  })

  p02 <- file.path(dirname(out), "P02")
  if (!file.exists(file.path(p02, "_SUCCESS"))) stop("P02 _SUCCESS is missing")
  p02_status <- read_tsv(file.path(p02, "STATUS.tsv"))
  if (nrow(p02_status) != 1L || !identical(p02_status$status[[1]], "PASS")) {
    stop("P02 STATUS is not uniquely PASS")
  }

  output_files <- file.path(out, c(
    "tables/gate_amendment.tsv", "tables/expression_coverage.tsv",
    "tables/neighborhood_celltype_support.tsv", "tables/gradient_relation_support.tsv",
    "tables/secreted_axis_candidate_matrix.tsv", "tables/frozen_secreted_axis.tsv",
    "tables/rejected_secreted_axes.tsv", "parameters.tsv", "sessionInfo.txt",
    "STATUS.tsv", "_SUCCESS"
  ))
  if (any(file.exists(output_files))) stop("Refusing to overwrite an amended P05 result")

  amendment <- data.frame(
    field = c("removed_hard_dependencies", "required_upstream", "retained_thresholds",
              "meringue_support", "statial_support", "effective_downstream"),
    value = c("P03_MERINGUE;P04_STATIAL", "P02", "ALL_P02_SECRETED_AXIS_THRESHOLDS",
              "NOT_REQUIRED_NOT_READ", "NOT_REQUIRED_NOT_READ", "P06_SECACT_AMENDED")
  )
  write_tsv(amendment, file.path(out, "tables/gate_amendment.tsv"))

  cc_path <- file.path(p02, "tables/existing_spatialcellchat_evidence.tsv")
  grad_path <- file.path(p02, "tables/existing_gradient_evidence.tsv")
  gate_path <- file.path(p02, "tables/secreted_axis_gate_definition.tsv")
  cc <- read_tsv(cc_path)
  grad <- read_tsv(grad_path)
  require_columns(cc, c(
    "biological_sample_id", "sender", "receiver", "ligand", "receptor", "score",
    "interaction_name", "pathway_name", "annotation", "passes_filter", "target_axis",
    "condition", "timepoint", "analysis_role", "stored_object"
  ), "P02 SpatialCellChat evidence")
  require_columns(grad, c(
    "biological_sample_id", "condition", "timepoint", "analysis_role", "cell_type",
    "corridor_spearman_rho", "evidence_source"
  ), "P02 gradient evidence")

  is_true <- function(x) isTRUE(x) || identical(toupper(as.character(x)), "TRUE")
  relevant <- cc[
    cc$condition == "MI" & cc$timepoint %in% c("Day3", "Day7") &
      cc$analysis_role == "primary_replicated" &
      cc$annotation == "Secreted Signaling" &
      vapply(cc$passes_filter, is_true, logical(1)) &
      !is.na(cc$target_axis) & nzchar(trimws(cc$target_axis)),
    , drop = FALSE
  ]
  if (!nrow(relevant)) stop("P02 contains no structurally eligible SpatialCellChat records")
  expected_samples <- c("MI_D3_rep1", "MI_D3_rep2", "MI_D3_rep3",
                        "MI_D7_rep1", "MI_D7_rep2", "MI_D7_rep3")
  if (!setequal(unique(relevant$biological_sample_id), expected_samples)) {
    stop("Filtered SpatialCellChat evidence is not the six primary MI tissues")
  }
  relevant$stage <- relevant$timepoint
  relevant$pathway <- relevant$pathway_name
  relevant$axis_id <- paste(relevant$sender, relevant$ligand, relevant$receptor,
                            relevant$receiver, relevant$pathway, sep = "|")
  duplicate_key <- paste(relevant$biological_sample_id, relevant$axis_id, sep = "|")
  if (anyDuplicated(duplicate_key)) stop("SpatialCellChat sample-axis keys are not unique")

  db <- CellChatDB.mouse
  if (anyDuplicated(db$interaction$interaction_name)) stop("CellChatDB interaction names are not unique")
  axis_meta <- do.call(rbind, lapply(split(relevant, relevant$axis_id), function(d) {
    fields <- c("sender", "ligand", "receptor", "receiver", "pathway", "interaction_name")
    vals <- lapply(fields, function(z) unique(as.character(d[[z]])))
    if (any(lengths(vals) != 1L)) stop("Axis metadata are not unique: ", d$axis_id[[1]])
    names(vals) <- fields
    db_i <- match(vals$interaction_name, db$interaction$interaction_name)
    if (is.na(db_i)) stop("CellChatDB interaction is missing: ", vals$interaction_name)
    db_row <- db$interaction[db_i, , drop = FALSE]
    if (!identical(vals$ligand, as.character(db_row$ligand[[1]])) ||
        !identical(vals$receptor, as.character(db_row$receptor[[1]]))) {
      stop("CellChatDB ligand/receptor mismatch for ", vals$interaction_name)
    }
    receptor_key <- as.character(db_row$receptor[[1]])
    if (receptor_key %in% rownames(db$complex)) {
      subunits <- trimws(as.character(unlist(db$complex[receptor_key, , drop = TRUE])))
      subunits <- unique(subunits[!is.na(subunits) & nzchar(subunits)])
    } else {
      subunits <- receptor_key
    }
    if (!length(subunits)) stop("No receptor subunits resolved for ", vals$interaction_name)
    data.frame(
      axis_id = d$axis_id[[1]], sender = vals$sender, ligand_mouse = vals$ligand,
      receptor = vals$receptor, receptor_subunits_mouse = paste(subunits, collapse = ";"),
      receiver = vals$receiver, pathway = vals$pathway,
      interaction_name = vals$interaction_name
    )
  }))
  rownames(axis_meta) <- NULL

  mapping <- read.csv(system.file("extdata", "Mouse2Human_filter.csv", package = "SecAct"),
                      row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(mapping$mouse)) stop("SecAct mouse symbols are not unique")
  signature_names <- strsplit(
    readLines(gzfile(system.file("extdata", "SecAct.tsv.gz", package = "SecAct")), n = 1L),
    "\t", fixed = TRUE
  )[[1]]
  signature_names <- unique(signature_names[nzchar(signature_names)])
  map_i <- match(axis_meta$ligand_mouse, mapping$mouse)
  axis_meta$ligand_secact_human <- ifelse(is.na(map_i), NA_character_, mapping$human[map_i])
  axis_meta$secact_mapping_unique <- !is.na(axis_meta$ligand_secact_human) &
    nzchar(axis_meta$ligand_secact_human)
  axis_meta$ligand_in_secact_signature <- axis_meta$secact_mapping_unique &
    axis_meta$ligand_secact_human %in% signature_names

  sample_stage <- setNames(c(rep("Day3", 3), rep("Day7", 3)), expected_samples)
  stage_samples <- split(expected_samples, sample_stage)
  coverage_rows <- list()
  neighborhood_rows <- list()
  ci <- 0L
  ni <- 0L

  for (sid in expected_samples) {
    stage <- sample_stage[[sid]]
    stored <- unique(relevant$stored_object[relevant$biological_sample_id == sid])
    if (length(stored) != 1L || !file.exists(stored)) stop("Missing unique stored object for ", sid)
    obj <- readRDS(stored)
    if (!inherits(obj, "Seurat") || !"Spatial" %in% SeuratObject::Assays(obj)) {
      stop("Invalid spatial object for ", sid)
    }
    if (!"counts" %in% SeuratObject::Layers(obj[["Spatial"]])) stop("Spatial/counts missing: ", sid)
    if (!"P8_major_dominant" %in% colnames(obj[[]])) stop("P8_major_dominant missing: ", sid)
    counts <- SeuratObject::LayerData(obj[["Spatial"]], layer = "counts")
    labels <- setNames(as.character(obj$P8_major_dominant), colnames(obj))
    if (!identical(colnames(counts), names(labels))) stop("Counts/metadata order mismatch: ", sid)

    rel_s <- relevant[relevant$biological_sample_id == sid, , drop = FALSE]
    metas <- axis_meta[match(rel_s$axis_id, axis_meta$axis_id), , drop = FALSE]
    requests <- do.call(rbind, lapply(seq_len(nrow(metas)), function(i) {
      rec <- strsplit(metas$receptor_subunits_mouse[[i]], ";", fixed = TRUE)[[1]]
      data.frame(
        axis_id = metas$axis_id[[i]],
        component_role = c("ligand", rep("receptor_subunit", length(rec))),
        component_gene = c(metas$ligand_mouse[[i]], rec),
        dominant_cell_type = c(metas$sender[[i]], rep(metas$receiver[[i]], length(rec)))
      )
    }))
    requests <- unique(requests)
    for (i in seq_len(nrow(requests))) {
      idx <- which(labels == requests$dominant_cell_type[[i]])
      denom <- length(idx)
      nnz <- if (denom && requests$component_gene[[i]] %in% rownames(counts)) {
        sum(as.numeric(counts[requests$component_gene[[i]], idx, drop = TRUE]) > 0)
      } else 0L
      rate <- if (denom) nnz / denom else NA_real_
      ci <- ci + 1L
      coverage_rows[[ci]] <- data.frame(
        biological_sample_id = sid, stage = stage, axis_id = requests$axis_id[[i]],
        component_role = requests$component_role[[i]], component_gene = requests$component_gene[[i]],
        dominant_cell_type = requests$dominant_cell_type[[i]], n_dominant_spots = denom,
        n_nonzero_spots = nnz, coverage_rate = rate,
        coverage_at_least_0_10 = is.finite(rate) && rate >= 0.10,
        at_least_two_nonzero_spots = nnz >= 2L,
        component_pass = is.finite(rate) && rate >= 0.10 && nnz >= 2L
      )
    }

    nb_path <- file.path(root, "11_domain_adjacency/sample_rds",
                         paste0(sid, "_BANKSY_adjacency.rds"))
    if (!file.exists(nb_path)) stop("Neighborhood object missing: ", sid)
    nb <- readRDS(nb_path)
    edge <- nb@tools$SpatialNeighborhood$methods$observed$edge_table
    require_columns(edge, c("cell", "neighbor"), paste("Neighborhood edge table", sid))
    prefix <- paste0(sid, "_")
    strip_prefix <- function(x) ifelse(startsWith(x, prefix), substring(x, nchar(prefix) + 1L), x)
    from_id <- strip_prefix(as.character(edge$cell))
    to_id <- strip_prefix(as.character(edge$neighbor))
    if (any(!from_id %in% names(labels)) || any(!to_id %in% names(labels))) {
      stop("Neighborhood barcodes do not map one-to-one to P8 labels: ", sid)
    }
    from_label <- unname(labels[from_id])
    to_label <- unname(labels[to_id])
    stage_axes <- axis_meta[axis_meta$axis_id %in%
                              unique(relevant$axis_id[relevant$stage == stage]), , drop = FALSE]
    for (i in seq_len(nrow(stage_axes))) {
      ni <- ni + 1L
      ec <- sum(from_label == stage_axes$sender[[i]] & to_label == stage_axes$receiver[[i]])
      neighborhood_rows[[ni]] <- data.frame(
        biological_sample_id = sid, stage = stage, axis_id = stage_axes$axis_id[[i]],
        sender = stage_axes$sender[[i]], receiver = stage_axes$receiver[[i]],
        edge_count = ec, sender_spot_count = sum(labels == stage_axes$sender[[i]]),
        receiver_spot_count = sum(labels == stage_axes$receiver[[i]]),
        directed_adjacency_present = ec > 0L, k = 6L,
        group_label = "P8_major_dominant", graph_source = "existing_SpatialNeighborhood_observed_edge_table"
      )
    }
    rm(obj, counts, nb, edge)
    gc(verbose = FALSE)
  }
  expression_coverage <- do.call(rbind, coverage_rows)
  neighborhood <- do.call(rbind, neighborhood_rows)
  write_tsv(expression_coverage, file.path(out, "tables/expression_coverage.tsv"))

  grad <- grad[
    grad$condition == "MI" & grad$timepoint %in% c("Day3", "Day7") &
      grad$analysis_role == "primary_replicated" & grad$evidence_source == "primary_discovery",
    , drop = FALSE
  ]
  rho_sign <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > 0, "+", ifelse(x < 0, "-", "0")))
  gradient_rows <- list()
  gi <- 0L
  candidate_keys <- unique(relevant[, c("axis_id", "stage")])
  for (j in seq_len(nrow(candidate_keys))) {
    aid <- candidate_keys$axis_id[[j]]
    stage <- candidate_keys$stage[[j]]
    am <- axis_meta[match(aid, axis_meta$axis_id), , drop = FALSE]
    map_gradient <- function(x) if (identical(x, "Macrophage")) "Mononuclear_myeloid" else x
    gs <- map_gradient(am$sender[[1]])
    gr <- map_gradient(am$receiver[[1]])
    local <- list()
    for (sid in stage_samples[[stage]]) {
      a <- grad[grad$biological_sample_id == sid & grad$cell_type == gs, , drop = FALSE]
      b <- grad[grad$biological_sample_id == sid & grad$cell_type == gr, , drop = FALSE]
      if (nrow(a) > 1L || nrow(b) > 1L) stop("Gradient rows are not unique for ", sid, " / ", aid)
      ar <- if (nrow(a)) as.numeric(a$corridor_spearman_rho[[1]]) else NA_real_
      br <- if (nrow(b)) as.numeric(b$corridor_spearman_rho[[1]]) else NA_real_
      sp <- if (is.finite(ar) && is.finite(br)) paste0(rho_sign(ar), ",", rho_sign(br)) else NA_character_
      local[[sid]] <- data.frame(
        biological_sample_id = sid, stage = stage, axis_id = aid,
        sender = am$sender[[1]], receiver = am$receiver[[1]],
        gradient_sender_label = gs, gradient_receiver_label = gr,
        sender_rho = ar, receiver_rho = br, ordered_sign_pair = sp
      )
    }
    local <- do.call(rbind, local)
    evaluable <- !is.na(local$ordered_sign_pair)
    tab <- table(local$ordered_sign_pair[evaluable])
    if (sum(evaluable) < 3L) {
      support <- "NOT_EVALUABLE"
      recurrent_pair <- NA_character_
      recurrent_n <- if (length(tab)) max(tab) else 0L
    } else {
      recurrent_n <- if (length(tab)) max(tab) else 0L
      recurrent_pair <- if (length(tab)) names(tab)[which.max(tab)] else NA_character_
      support <- if (recurrent_n >= 2L) "TRUE" else "FALSE"
    }
    local$n_evaluable_tissues <- sum(evaluable)
    local$recurrent_sign_pair <- recurrent_pair
    local$n_recurrent_sign_pair <- recurrent_n
    local$gradient_support <- support
    gi <- gi + 1L
    gradient_rows[[gi]] <- local
  }
  gradient_relation <- do.call(rbind, gradient_rows)

  nb_summary <- do.call(rbind, lapply(split(neighborhood, paste(neighborhood$axis_id, neighborhood$stage, sep = "||")), function(d) {
    data.frame(axis_id = d$axis_id[[1]], stage = d$stage[[1]],
               n_tissues_with_directed_adjacency = sum(d$directed_adjacency_present),
               neighborhood_support = sum(d$directed_adjacency_present) >= 2L)
  }))
  neighborhood <- merge(neighborhood, nb_summary, by = c("axis_id", "stage"), all.x = TRUE, sort = FALSE)
  write_tsv(neighborhood, file.path(out, "tables/neighborhood_celltype_support.tsv"))
  write_tsv(gradient_relation, file.path(out, "tables/gradient_relation_support.tsv"))

  candidate_rows <- list()
  for (j in seq_len(nrow(candidate_keys))) {
    aid <- candidate_keys$axis_id[[j]]
    stage <- candidate_keys$stage[[j]]
    d <- relevant[relevant$axis_id == aid & relevant$stage == stage, , drop = FALSE]
    am <- axis_meta[match(aid, axis_meta$axis_id), , drop = FALSE]
    tissues <- sort(unique(d$biological_sample_id))
    ec <- expression_coverage[expression_coverage$axis_id == aid &
                                expression_coverage$biological_sample_id %in% tissues, , drop = FALSE]
    nb <- nb_summary[nb_summary$axis_id == aid & nb_summary$stage == stage, , drop = FALSE]
    gd <- gradient_relation[gradient_relation$axis_id == aid & gradient_relation$stage == stage, , drop = FALSE]
    gradient_support <- unique(gd$gradient_support)
    if (length(gradient_support) != 1L || nrow(nb) != 1L) stop("Support summaries are not unique: ", aid, " / ", stage)
    n_support <- length(tissues)
    cellchat_pass <- n_support >= 2L
    coverage_pass <- nrow(ec) > 0L && all(ec$component_pass)
    min_coverage <- if (nrow(ec) && all(is.finite(ec$coverage_rate))) min(ec$coverage_rate) else NA_real_
    additional <- isTRUE(nb$neighborhood_support[[1]]) || identical(gradient_support, "TRUE")
    map_pass <- isTRUE(am$secact_mapping_unique[[1]]) && isTRUE(am$ligand_in_secact_signature[[1]])
    reasons <- c(
      if (!cellchat_pass) "same_ordered_axis_recurrence_below_2_of_3",
      if (!coverage_pass) "expression_coverage_or_two_spot_gate_failed",
      if (!additional) "no_evaluable_passing_additional_spatial_layer",
      if (!isTRUE(am$secact_mapping_unique[[1]])) "ligand_has_no_unique_secact_ortholog",
      if (isTRUE(am$secact_mapping_unique[[1]]) && !isTRUE(am$ligand_in_secact_signature[[1]])) "ligand_not_in_secact_signature"
    )
    candidate_rows[[j]] <- data.frame(
      axis_id = aid, sender = am$sender[[1]], ligand_mouse = am$ligand_mouse[[1]],
      ligand_secact_human = am$ligand_secact_human[[1]], receptor = am$receptor[[1]],
      receptor_subunits_mouse = am$receptor_subunits_mouse[[1]], receiver = am$receiver[[1]],
      pathway = am$pathway[[1]], interaction_name = am$interaction_name[[1]], stage = stage,
      supporting_tissues = paste(tissues, collapse = ";"), n_same_ordered_axis = n_support,
      direction_concordance = n_support / 3, spatialcellchat_support = cellchat_pass,
      spatialcellchat_median_score = median(as.numeric(d$score)),
      expression_coverage_pass = coverage_pass, minimum_component_coverage = min_coverage,
      neighborhood_support = nb$neighborhood_support[[1]],
      n_neighborhood_supporting_tissues = nb$n_tissues_with_directed_adjacency[[1]],
      gradient_support = gradient_support,
      n_gradient_evaluable_tissues = unique(gd$n_evaluable_tissues)[[1]],
      n_gradient_recurrent_sign_pair = unique(gd$n_recurrent_sign_pair)[[1]],
      additional_spatial_layer_pass = additional,
      secact_mapping_unique = am$secact_mapping_unique[[1]],
      ligand_in_secact_signature = am$ligand_in_secact_signature[[1]],
      stage_pass = cellchat_pass && coverage_pass && additional && map_pass,
      rejection_reason = if (length(reasons)) paste(reasons, collapse = ";") else ""
    )
  }
  candidate_matrix <- do.call(rbind, candidate_rows)
  rownames(candidate_matrix) <- NULL
  write_tsv(candidate_matrix, file.path(out, "tables/secreted_axis_candidate_matrix.tsv"))

  passed <- candidate_matrix[candidate_matrix$stage_pass, , drop = FALSE]
  frozen_all <- if (nrow(passed)) do.call(rbind, lapply(split(passed, passed$axis_id), function(d) {
    aid <- d$axis_id[[1]]
    source <- relevant[relevant$axis_id == aid & relevant$stage %in% d$stage, , drop = FALSE]
    data.frame(
      axis_id = aid, sender = d$sender[[1]], ligand_mouse = d$ligand_mouse[[1]],
      ligand_secact_human = d$ligand_secact_human[[1]], receptor = d$receptor[[1]],
      receptor_subunits_mouse = d$receptor_subunits_mouse[[1]], receiver = d$receiver[[1]],
      pathway = d$pathway[[1]], interaction_name = d$interaction_name[[1]],
      supported_stages = paste(sort(unique(d$stage)), collapse = ";"),
      n_supported_stages = length(unique(d$stage)),
      discovery_supporting_tissues = paste(sort(unique(source$biological_sample_id)), collapse = ";"),
      n_discovery_supporting_tissues = sum(d$n_same_ordered_axis),
      n_non_cellchat_spatial_layers = as.integer(any(d$neighborhood_support)) +
        as.integer(any(d$gradient_support == "TRUE")),
      minimum_component_coverage = min(d$minimum_component_coverage),
      spatialcellchat_median_score = median(as.numeric(source$score)),
      external_used_for_selection = FALSE
    )
  })) else data.frame()

  if (nrow(frozen_all)) {
    ord <- order(-frozen_all$n_supported_stages,
                 -frozen_all$n_discovery_supporting_tissues,
                 -frozen_all$n_non_cellchat_spatial_layers,
                 -frozen_all$minimum_component_coverage,
                 -frozen_all$spatialcellchat_median_score,
                 frozen_all$axis_id)
    frozen_all <- frozen_all[ord, , drop = FALSE]
    frozen_all$selection_rank <- seq_len(nrow(frozen_all))
    frozen <- frozen_all[frozen_all$selection_rank <= 10L, , drop = FALSE]
  } else {
    frozen <- data.frame(
      axis_id = character(), sender = character(), ligand_mouse = character(),
      ligand_secact_human = character(), receptor = character(),
      receptor_subunits_mouse = character(), receiver = character(), pathway = character(),
      interaction_name = character(), supported_stages = character(),
      n_supported_stages = integer(), discovery_supporting_tissues = character(),
      n_discovery_supporting_tissues = integer(), n_non_cellchat_spatial_layers = integer(),
      minimum_component_coverage = numeric(), spatialcellchat_median_score = numeric(),
      external_used_for_selection = logical(), selection_rank = integer()
    )
  }
  write_tsv(frozen, file.path(out, "tables/frozen_secreted_axis.tsv"))

  rejected_stage <- candidate_matrix[!candidate_matrix$stage_pass,
                                     c("axis_id", "stage", "rejection_reason"), drop = FALSE]
  rejected_stage$selection_rank <- NA_integer_
  if (nrow(frozen_all)) {
    rejected_rank <- frozen_all[frozen_all$selection_rank > 10L, c("axis_id", "supported_stages", "selection_rank"), drop = FALSE]
    if (nrow(rejected_rank)) {
      names(rejected_rank)[[2]] <- "stage"
      rejected_rank$rejection_reason <- "eligible_but_below_fixed_top_10_rank"
      rejected_rank <- rejected_rank[, names(rejected_stage), drop = FALSE]
      rejected_stage <- rbind(rejected_stage, rejected_rank)
    }
  }
  write_tsv(rejected_stage, file.path(out, "tables/rejected_secreted_axes.tsv"))

  p02_gate <- read_tsv(gate_path)
  parameters <- rbind(
    data.frame(parameter = paste0("P02_", p02_gate$criterion), value = p02_gate$frozen_rule),
    data.frame(parameter = c(
      "operation", "required_upstream", "removed_hard_dependencies", "direction_interpretation",
      "expression_layer", "coverage_denominator", "coverage_threshold", "minimum_nonzero_spots",
      "neighborhood_source", "neighborhood_k", "gradient_crosswalk", "additional_spatial_layer",
      "SecAct_version", "SecAct_signature", "maximum_axes", "external_biology_read",
      "meringue_support", "statial_support", threads
    ), value = c(
      "P05_SECRETED_AXIS_GATE_AMENDED", "P02", "P03_MERINGUE;P04_STATIAL",
      "same_exact_ordered_axis_recurrence", "Spatial/counts",
      "ligand_sender_dominant_spots;receptor_subunits_receiver_dominant_spots",
      "0.10", "2", "existing_SpatialNeighborhood_observed_edge_table", "6",
      "Macrophage=Mononuclear_myeloid;all_others_exact",
      "Neighborhood_TRUE_OR_gradient_TRUE", as.character(packageVersion("SecAct")),
      "SecAct", "10", "FALSE", "NOT_REQUIRED_NOT_READ", "NOT_REQUIRED_NOT_READ",
      rep("1", length(threads))
    ))
  )
  parameter_hash <- sha256_value(parameters)
  parameters$parameter_hash <- parameter_hash
  write_tsv(parameters, file.path(out, "parameters.tsv"))
  write_session_info(file.path(out, "sessionInfo.txt"))

  reread_candidate <- read_tsv(file.path(out, "tables/secreted_axis_candidate_matrix.tsv"))
  reread_frozen <- read_tsv(file.path(out, "tables/frozen_secreted_axis.tsv"))
  if (nrow(reread_candidate) != nrow(candidate_matrix) || nrow(reread_frozen) > 10L ||
      anyDuplicated(reread_frozen$axis_id)) stop("P05 output integrity gate failed")
  if (nrow(reread_frozen) && !identical(reread_frozen$selection_rank, seq_len(nrow(reread_frozen)))) {
    stop("Frozen axis rank integrity failed")
  }
  status <- if (nrow(frozen)) "PASS" else "NEGATIVE_COMPLETE"
  verdict <- if (nrow(frozen)) "PASS_AXIS_FROZEN" else "NEGATIVE_COMPLETE_NO_AXIS"
  input_hash <- sha256_value(c(sha256_file(cc_path), sha256_file(grad_path), sha256_file(gate_path)))
  write_phase_status(
    out = out, phase = "P05_SECRETED_AXIS_GATE_AMENDED", status = status,
    started_at = started_at, input_hash = input_hash, parameter_hash = parameter_hash,
    n_expected = 6L, n_completed = 6L, producer_functions = "NONE_REUSE_ONLY",
    biological_verdict = verdict, blocker = "", output_root = out
  )
  if (!file.create(file.path(out, "_SUCCESS"))) stop("Could not create P05 _SUCCESS")
  message("P05_EVENT status=", status, " frozen_axes=", nrow(frozen),
          " candidates=", length(unique(candidate_matrix$axis_id)),
          " external_used_for_selection=FALSE")
}

started_at <- now_utc()
tryCatch(main(), error = function(e) {
  write_failed(e, started_at)
  message("P05_EVENT status=FAILED error=", conditionMessage(e))
  quit(save = "no", status = 1L)
})
