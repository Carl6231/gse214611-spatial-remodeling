#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Usage: P06_run_amended.R ROOT OUT SPACET_LIB")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- normalizePath(args[[2]], mustWork = TRUE)
spacet_lib <- normalizePath(args[[3]], mustWork = TRUE)
expected_out <- file.path(root, "22_JTM_completion_20261014/P06_SECACT_compat2")
if (!identical(out, expected_out)) stop("Output is outside the authorized amended P06 directory")
.libPaths(c(spacet_lib, .libPaths()))
source(file.path(root, "21_scop_expansion/P01/script/common.R"))

phase_started <- now_utc()
threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")

empty_status <- function() data.frame(
  `function` = character(), sample_id = character(), stage = character(), modality = character(),
  target = character(), signal_mode = character(), status = character(), called = logical(),
  reason = character(), evidence_path = character(), biological_unit = character(),
  check.names = FALSE
)

write_gate_terminal <- function(status, verdict, reason, function_status) {
  write_tsv(function_status, file.path(out, "functions_final_status.tsv"))
  write_session_info(file.path(out, "sessionInfo.txt"))
  write_phase_status(
    out = out, phase = "P06_SECACT_AMENDED", status = status,
    started_at = phase_started, input_hash = NA_character_, parameter_hash = NA_character_,
    n_expected = 12L, n_completed = 0L, producer_functions = "NONE",
    biological_verdict = verdict, blocker = reason, output_root = out
  )
}

write_failure <- function(e, statuses = empty_status()) {
  dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
  writeLines(conditionMessage(e), file.path(out, "logs/error_message.txt"))
  if (nrow(statuses)) try(write_tsv(statuses, file.path(out, "functions_final_status.tsv")), silent = TRUE)
  try(write_session_info(file.path(out, "sessionInfo.txt")), silent = TRUE)
  try(write_phase_status(
    out = out, phase = "P06_SECACT_AMENDED", status = "FAILED",
    started_at = phase_started, input_hash = NA_character_, parameter_hash = NA_character_,
    n_expected = 12L,
    n_completed = length(unique(statuses$sample_id[statuses$status %in% c("PASS", "NEGATIVE_COMPLETE")])),
    producer_functions = paste(unique(statuses[["function"]][statuses$called]), collapse = ";"),
    biological_verdict = "FAILED_AFTER_FIXED_PILOT_OR_FULL_RUN_CALL",
    blocker = conditionMessage(e), output_root = out
  ), silent = TRUE)
}

main <- function() {
  if (!all(Sys.getenv(threads) == "1")) stop("Thread environment is not fixed to 1")
  suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
    library(scop)
    library(SecAct)
    library(ggplot2)
    library(NMF)
    library(png)
    library(digest)
  })
  if (!identical(as.character(packageVersion("scop")), "0.9.1")) stop("scop version drift")
  if (!identical(as.character(packageVersion("SecAct")), "1.1.0")) stop("SecAct version drift")
  producers <- c("RunSecAct", "RunSecActCCC", "RunSecActSignalingPattern",
                 "RunSecActPatternGenes", "RunSecActVelocity")
  if (!all(producers %in% getNamespaceExports("scop"))) stop("A required SCOP SecAct producer is not exported")

  source(file.path(out, "script/secact_compat.R"))
  NMF::nmf.options(pbackend="seq")
  foreach::registerDoSEQ()
  result_files <- file.path(out, c(
    "tables/ortholog_mapping.tsv", "tables/secact_activity.tsv.gz", "tables/secact_ccc.tsv",
    "tables/secact_patterns.tsv", "tables/secact_pattern_genes.tsv",
    "tables/secact_velocity_sample_effects.tsv", "tables/secact_discovery_summary.tsv",
    "tables/pilot_output_audit.tsv", "functions_final_status.tsv", "parameters.tsv",
    "sessionInfo.txt", "STATUS.tsv", "_SUCCESS"
  ))
  existing_results <- result_files[file.exists(result_files)]
  if (length(existing_results)) {
    prior_status_path <- file.path(out, "STATUS.tsv")
    prior_functions_path <- file.path(out, "functions_final_status.tsv")
    resumable_gate_only <- file.exists(prior_status_path) && file.exists(prior_functions_path) &&
      !file.exists(file.path(out, "_SUCCESS")) &&
      all(basename(existing_results) %in% c("parameters.tsv", "sessionInfo.txt", "STATUS.tsv",
                                           "functions_final_status.tsv"))
    if (resumable_gate_only) {
      prior_status <- read_tsv(prior_status_path)
      prior_functions <- read_tsv(prior_functions_path)
      resumable_gate_only <- nrow(prior_status) == 1L &&
        identical(prior_status$status[[1]], "BLOCKED_DEPENDENCY") &&
        identical(as.integer(prior_status$n_samples_completed[[1]]), 0L) &&
        nrow(prior_functions) == 5L && !any(prior_functions$called) &&
        all(prior_functions$status == "BLOCKED_DEPENDENCY") &&
        grepl("SpaCET", prior_status$blocker[[1]], fixed = TRUE)
    }
    if (!resumable_gate_only) stop("Refusing to overwrite or merge an amended P06 result")
    message("P06_EVENT resume=ONE_NARROW_SPACET_DEPENDENCY_RETRY prior_producer_calls=0")
  }

  parameters <- data.frame(
    parameter = c(
      "signature", "input_layer", "scale.factor", "nrand", "lambda", "is.filter.sig",
      "is.group.sig", "is.group.cor", "backend", "rng_method", "cores", "spatial_radius_um",
      "pattern_k_candidates", "pattern_native_seed", "velocity_modes", "biological_unit",
      "source_species", "modeling_namespace", "external_biology_read", "Statial_required",
      threads
    ),
    value = c(
      "SecAct", "raw counts", "1e5", "1000", "5e5", "FALSE", "TRUE", "0.9",
      "cpu-pure", "mt19937", "1", "200", "2:5", "123456", "sending;receiving",
      "physical_tissue/snRNA_sample_id", "mouse", "human_ortholog", "FALSE", "FALSE",
      rep("1", length(threads))
    )
  )
  parameter_hash <- sha256_value(parameters)
  parameters$parameter_hash <- parameter_hash
  write_tsv(parameters, file.path(out, "parameters.tsv"))

  p05 <- file.path(root, "21_scop_expansion/P05_SECRETED_AXIS_GATE_AMENDED")
  p05_status_path <- file.path(p05, "STATUS.tsv")
  spacet_runtime <- dirname(dirname(spacet_lib))
  spacet_status_path <- file.path(spacet_runtime, "STATUS.tsv")
  if (!file.exists(p05_status_path)) {
    fs <- data.frame(
      `function` = producers, sample_id = "ALL", stage = "Day3;Day7", modality = "snRNA;Visium",
      target = "", signal_mode = "", status = "BLOCKED_DEPENDENCY", called = FALSE,
      reason = "P05 amended status is missing", evidence_path = p05_status_path,
      biological_unit = "physical_tissue/biological_sample", check.names = FALSE
    )
    write_gate_terminal("BLOCKED_DEPENDENCY", "NOT_EVALUATED_P05_UNAVAILABLE", fs$reason[[1]], fs)
    message("P06_EVENT status=BLOCKED_DEPENDENCY reason=P05_STATUS_MISSING")
    return(invisible(NULL))
  }
  p05s <- read_tsv(p05_status_path)
  if (nrow(p05s) != 1L) stop("P05 amended STATUS is not unique")
  if (identical(p05s$status[[1]], "NEGATIVE_COMPLETE")) {
    fs <- data.frame(
      `function` = producers, sample_id = "ALL", stage = "Day3;Day7", modality = "snRNA;Visium",
      target = "", signal_mode = "", status = "NOT_APPLICABLE", called = FALSE,
      reason = "P05 completed with no eligible frozen axis", evidence_path = p05_status_path,
      biological_unit = "physical_tissue/biological_sample", check.names = FALSE
    )
    write_gate_terminal("NOT_APPLICABLE", "NOT_APPLICABLE_NO_FROZEN_SECRETED_AXIS", fs$reason[[1]], fs)
    message("P06_EVENT status=NOT_APPLICABLE reason=NO_FROZEN_AXIS")
    return(invisible(NULL))
  }
  if (!identical(p05s$status[[1]], "PASS") || !file.exists(file.path(p05, "_SUCCESS"))) {
    fs <- data.frame(
      `function` = producers, sample_id = "ALL", stage = "Day3;Day7", modality = "snRNA;Visium",
      target = "", signal_mode = "", status = "BLOCKED_DEPENDENCY", called = FALSE,
      reason = paste0("P05 is not completed PASS: ", p05s$status[[1]]), evidence_path = p05_status_path,
      biological_unit = "physical_tissue/biological_sample", check.names = FALSE
    )
    write_gate_terminal("BLOCKED_DEPENDENCY", "NOT_EVALUATED_P05_NOT_PASS", fs$reason[[1]], fs)
    message("P06_EVENT status=BLOCKED_DEPENDENCY reason=P05_NOT_PASS")
    return(invisible(NULL))
  }
  spacet_success <- file.path(spacet_runtime, "_SUCCESS")
  spacet_error <- file.path(spacet_runtime, "error_message.txt")
  if (!file.exists(spacet_success) && file.exists(spacet_error)) {
    first_error <- paste(readLines(spacet_error, warn = FALSE), collapse = " ")
    if (grepl("JEFworks-Lab/MUDAN", first_error, fixed = TRUE) &&
        grepl("HTTP error 404", first_error, fixed = TRUE)) {
      message("SPACET_EVENT retry=1 narrow_fix=correct_MUDAN_repository_owner commit_unchanged")
      retry_rc <- system2(
        file.path(R.home("bin"), "Rscript"),
        c(file.path(out, "script/install_spacet.R"), spacet_lib, spacet_runtime)
      )
      if (!identical(retry_rc, 0L)) message("SPACET_EVENT retry=1 exit=", retry_rc)
    }
  }
  if (!file.exists(spacet_status_path) || !file.exists(spacet_success)) {
    reason <- "SpaCET isolated installation or smoke test did not pass"
    fs <- data.frame(
      `function` = producers, sample_id = "ALL", stage = "Day3;Day7", modality = "snRNA;Visium",
      target = "", signal_mode = "", status = "BLOCKED_DEPENDENCY", called = FALSE,
      reason = reason, evidence_path = spacet_status_path,
      biological_unit = "physical_tissue/biological_sample", check.names = FALSE
    )
    write_gate_terminal("BLOCKED_DEPENDENCY", "NOT_EVALUATED_SPACET_UNAVAILABLE", reason, fs)
    message("P06_EVENT status=BLOCKED_DEPENDENCY reason=SPACET_NOT_PASS")
    return(invisible(NULL))
  }
  spacet_status <- read_tsv(spacet_status_path)
  if (nrow(spacet_status) != 1L || !identical(spacet_status$status[[1]], "PASS")) {
    stop("SpaCET runtime status is not uniquely PASS")
  }
  suppressPackageStartupMessages(library(SpaCET, lib.loc = spacet_lib))
  if (!identical(as.character(packageVersion("SpaCET")), "1.4.0")) stop("SpaCET version drift")

  axes_path <- file.path(p05, "tables/frozen_secreted_axis.tsv")
  axes <- read_tsv(axes_path)
  require_columns(axes, c("axis_id", "ligand_mouse", "ligand_secact_human", "supported_stages"),
                  "Frozen secreted axes")
  if (!nrow(axes) || nrow(axes) > 10L || anyDuplicated(axes$axis_id)) stop("Frozen axis count/uniqueness gate failed")
  ligands <- sort(unique(axes$ligand_secact_human))
  if (!length(ligands) || anyNA(ligands) || any(!nzchar(ligands))) stop("Frozen SecAct ligand namespace is invalid")

  mapping <- read.csv(system.file("extdata", "Mouse2Human_filter.csv", package = "SecAct"),
                      row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(mapping$mouse)) stop("SecAct mouse mapping input is not unique")
  signature_names <- strsplit(
    readLines(gzfile(system.file("extdata", "SecAct.tsv.gz", package = "SecAct")), n = 1L),
    "\t", fixed = TRUE
  )[[1]]
  signature_names <- unique(signature_names[nzchar(signature_names)])
  mapping$in_secact_signature <- mapping$human %in% signature_names
  mapping$mapping_source <- "SecAct_1.1.0_Mouse2Human_filter.csv"
  write_tsv(mapping, file.path(out, "tables/ortholog_mapping.tsv"))
  map_counts <- getFromNamespace("mouse2human_mat", "SecAct")
  if (!is.function(map_counts)) stop("SecAct mouse2human_mat is unavailable")

  statuses <- empty_status()
  activity_rows <- data.frame(
    sample_id = character(), stage = character(), modality = character(),
    celltype_or_scope = character(), secretedProtein = character(), beta = numeric(),
    se = numeric(), zscore = numeric(), pvalue = numeric(), aggregation = character(),
    n_input_units = integer(), inference_scope = character()
  )
  ccc_rows <- data.frame(
    record_type = character(), sample_id = character(), stage = character(),
    sender = character(), secretedProtein = character(), receiver = character(),
    overall_strength = numeric(), status = character(), n_nuclei = integer(), reason = character()
  )
  pattern_rows <- data.frame(
    sample_id = character(), stage = character(), pattern_index = integer(),
    secretedProtein = character(), pattern_weight = numeric(), pattern_scope = character()
  )
  pattern_gene_rows <- data.frame(
    sample_id = character(), stage = character(), pattern_index = integer(), rank = integer(),
    secretedProtein = character(), pattern_weight = numeric(), status = character()
  )
  velocity_rows <- data.frame(
    sample_id = character(), stage = character(), axis_id = character(),
    ligand_mouse = character(), ligand_secact_human = character(), signal_mode = character(),
    n_vectors = integer(), median_vector_length = numeric(), IQR_vector_length = numeric(),
    median_dx = numeric(), median_dy = numeric(), status = character()
  )
  pilot_audit <- data.frame(
    component = character(), sample_id = character(), status = character(), pass = logical(),
    detail = character(), evidence_path = character()
  )

  append_status <- function(fun, sample_id, stage, modality, status, called,
                            target = "", signal_mode = "", reason = "", evidence_path = "") {
    row <- data.frame(
      `function` = fun, sample_id = sample_id, stage = stage, modality = modality,
      target = target, signal_mode = signal_mode, status = status, called = called,
      reason = reason, evidence_path = evidence_path,
      biological_unit = if (modality == "Visium") "physical_tissue" else "snRNA_sample_id",
      check.names = FALSE
    )
    key <- paste(row[["function"]], row$sample_id, row$target, row$signal_mode, sep = "|")
    old_key <- paste(statuses[["function"]], statuses$sample_id, statuses$target, statuses$signal_mode, sep = "|")
    if (key %in% old_key) stop("Duplicate final function status key: ", key)
    statuses <<- rbind(statuses, row)
    write_tsv(statuses, file.path(out, "functions_final_status.tsv"))
  }

  append_audit <- function(component, sample_id, status, pass, detail, evidence_path) {
    pilot_audit <<- rbind(pilot_audit, data.frame(
      component = component, sample_id = sample_id, status = status, pass = pass,
      detail = detail, evidence_path = evidence_path
    ))
    write_tsv(pilot_audit, file.path(out, "tables/pilot_output_audit.tsv"))
  }

  write_gz_tsv <- function(x, path) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }

  flush_tables <- function() {
    write_gz_tsv(activity_rows, file.path(out, "tables/secact_activity.tsv.gz"))
    write_tsv(ccc_rows, file.path(out, "tables/secact_ccc.tsv"))
    write_tsv(pattern_rows, file.path(out, "tables/secact_patterns.tsv"))
    write_tsv(pattern_gene_rows, file.path(out, "tables/secact_pattern_genes.tsv"))
    write_tsv(velocity_rows, file.path(out, "tables/secact_velocity_sample_effects.tsv"))
  }

  extract_activity <- function(obj, modality) {
    res <- if (modality == "snRNA") obj@tools$SecAct$activity else
      obj@results$SecAct_output$SecretedProteinActivity
    required <- c("beta", "se", "zscore", "pvalue")
    if (!is.list(res) || !all(required %in% names(res))) stop("SecAct activity structure is incomplete")
    dims <- lapply(res[required], dim)
    if (!all(vapply(dims, identical, logical(1), dims[[1]]))) stop("SecAct activity dimensions differ")
    for (nm in required) {
      if (!length(res[[nm]]) || any(!is.finite(res[[nm]]))) stop("Non-finite or empty activity matrix: ", nm)
    }
    res
  }

  activity_table <- function(res, sample_id, stage, modality, n_by_column) {
    rn <- rownames(res$zscore)
    cn <- colnames(res$zscore)
    if (is.null(rn) || is.null(cn)) stop("SecAct activity dimnames are missing")
    if (modality == "snRNA") {
      out_df <- data.frame(
        sample_id = sample_id, stage = stage, modality = modality,
        celltype_or_scope = rep(cn, each = length(rn)),
        secretedProtein = rep(rn, times = length(cn)),
        beta = as.vector(res$beta), se = as.vector(res$se),
        zscore = as.vector(res$zscore), pvalue = as.vector(res$pvalue),
        aggregation = "native_celltype_aggregate",
        n_input_units = as.integer(rep(n_by_column[cn], each = length(rn))),
        inference_scope = "within_sample_algorithmic_only"
      )
    } else {
      med <- function(x) apply(x, 1L, median)
      out_df <- data.frame(
        sample_id = sample_id, stage = stage, modality = modality,
        celltype_or_scope = "all_spots_tissue_median", secretedProtein = rn,
        beta = med(res$beta), se = med(res$se), zscore = med(res$zscore),
        pvalue = med(res$pvalue), aggregation = "median_across_spots_within_tissue",
        n_input_units = ncol(res$zscore), inference_scope = "within_tissue_descriptive_only"
      )
    }
    out_df
  }

  sn_path <- file.path(root, "03_snrna_reference/rctd_reference_8types_robust.rds")
  if (!file.exists(sn_path)) stop("Authoritative snRNA RDS is missing")
  sn <- readRDS(sn_path)
  if (!inherits(sn, "Seurat") || !"RNA" %in% Assays(sn) || !"counts" %in% Layers(sn[["RNA"]])) {
    stop("snRNA object lacks RNA/counts")
  }
  require_columns(sn[[]], c("sample_id", "timepoint", "celltype"), "snRNA metadata")
  sn_samples <- c("SRR21777596", "SRR21777597", "SRR21777598",
                  "SRR21777581", "SRR21777582", "SRR21777599")
  if (!all(sn_samples %in% unique(sn$sample_id))) stop("One or more fixed snRNA samples are missing")
  if (sum(sn$sample_id == "SRR21777597") != 746L) stop("Pilot snRNA sample is not exactly 746 nuclei")
  mouse_sn_counts <- LayerData(sn[["RNA"]], layer = "counts")
  if (any(mouse_sn_counts@x < 0) || any(abs(mouse_sn_counts@x - round(mouse_sn_counts@x)) > 1e-8)) {
    stop("snRNA raw counts are not non-negative integers")
  }
  human_sn_counts <- map_counts(mouse_sn_counts)
  if (anyDuplicated(rownames(human_sn_counts)) || !identical(colnames(human_sn_counts), colnames(sn))) {
    stop("snRNA ortholog conversion changed cell order or produced duplicate genes")
  }

  build_sn_sample <- function(sample_id) {
    cells <- colnames(sn)[sn$sample_id == sample_id]
    meta <- sn[[]][cells, , drop = FALSE]
    obj <- CreateSeuratObject(counts = human_sn_counts[, cells, drop = FALSE], assay = "RNA", meta.data = meta)
    if (!identical(colnames(obj), cells)) stop("Derived snRNA object changed cell order: ", sample_id)
    obj@misc$P06_input_contract <- list(
      source_species = "mouse", modeling_namespace = "human_ortholog",
      input_layer = "RNA/counts", biological_unit = "snRNA_sample_id"
    )
    obj
  }

  run_sn <- function(sample_id, pilot = FALSE) {
    stage <- unique(as.character(sn$timepoint[sn$sample_id == sample_id]))
    if (length(stage) != 1L) stop("snRNA sample has non-unique timepoint: ", sample_id)
    base <- build_sn_sample(sample_id)
    group_n <- table(base$celltype)
    activity_path <- file.path(out, "objects", paste0(sample_id, "_snRNA_RunSecAct.rds"))
    act_obj <- if (file.exists(activity_path)) {
      message("Reusing validated snRNA activity: ",sample_id)
      readRDS(activity_path)
    } else tryCatch(
      RunSecAct_compat(
        srt = base, mode = "scRNAseq", assay = "RNA", cellType_meta = "celltype",
        is.singleCellLevel = FALSE, sigMatrix = "SecAct", is.filter.sig = FALSE,
        is.group.sig = TRUE, is.group.cor = 0.9, lambda = 5e5, nrand = 1000,
        cores = 1L, backend = "cpu-pure", rng_method = "mt19937",
        store_assay = FALSE, store_results = TRUE, return_seurat = TRUE
      ),
      error = function(e) {
        append_status("RunSecAct", sample_id, stage, "snRNA", "FAILED", TRUE,
                      reason = conditionMessage(e))
        stop(e)
      }
    )
    act <- extract_activity(act_obj, "snRNA")
    saveRDS(act_obj, activity_path)
    reloaded <- readRDS(activity_path)
    extract_activity(reloaded, "snRNA")
    activity_rows <<- rbind(activity_rows, activity_table(act, sample_id, stage, "snRNA", group_n))
    append_status("RunSecAct", sample_id, stage, "snRNA", "PASS", TRUE,
                  evidence_path = activity_path)
    if (pilot) append_audit("RunSecAct_scRNAseq", sample_id, "PASS", TRUE,
                            paste(dim(act$zscore), collapse = "x"), activity_path)

    eligibility <- data.frame(
      record_type = "celltype_eligibility", sample_id = sample_id, stage = stage,
      sender = names(group_n), secretedProtein = "", receiver = "", overall_strength = NA_real_,
      status = ifelse(group_n >= 30L, "PASS", "NOT_APPLICABLE"),
      n_nuclei = as.integer(group_n),
      reason = ifelse(group_n >= 30L, "eligible_for_native_CCC_source_step", "fewer_than_30_nuclei")
    )
    ccc_rows <<- rbind(ccc_rows, eligibility)
    ccc_path <- file.path(out, "objects", paste0(sample_id, "_snRNA_RunSecActCCC.rds"))
    if (sum(group_n >= 30L) < 2L) {
      append_status("RunSecActCCC", sample_id, stage, "snRNA", "NOT_APPLICABLE", FALSE,
                    reason = "Fewer than two cell types have at least 30 nuclei")
      if (pilot) append_audit("RunSecActCCC_scRNAseq", sample_id, "NOT_APPLICABLE", FALSE,
                              "fewer_than_two_eligible_celltypes", "")
      flush_tables()
      return(invisible(NULL))
    }
    ccc_obj <- if (file.exists(ccc_path)) {
      message("Reusing validated snRNA CCC: ",sample_id)
      readRDS(ccc_path)
    } else tryCatch(
      RunSecActCCC_compat(
        srt = base, mode = "scRNAseq", cellType_meta = "celltype", condition_meta = NULL,
        scale.factor = 1e5, act_diff_cutoff = 2, exp_logFC_cutoff = 0.2,
        exp_mean_all_cutoff = 2, exp_fraction_case_cutoff = 0.1, padj_cutoff = 0.01,
        sigMatrix = "SecAct", is.group.sig = TRUE, is.group.cor = 0.9,
        lambda = 5e5, nrand = 1000, coreNo = 1, store_results = TRUE
      ),
      error = function(e) {
        append_status("RunSecActCCC", sample_id, stage, "snRNA", "FAILED", TRUE,
                      reason = conditionMessage(e))
        stop(e)
      }
    )
    ccc <- ccc_obj@misc$SecAct_output$SecretedProteinCCC
    if (is.null(ccc) || !is.data.frame(ccc)) stop("RunSecActCCC did not return its expected data.frame")
    if (!"SecAct_CCC" %in% names(ccc_obj@tools)) stop("RunSecActCCC SCOP tool record is missing")
    ccc_status <- if (nrow(ccc)) "PASS" else "NEGATIVE_COMPLETE"
    if (nrow(ccc)) {
      needed <- c("sender", "secretedProtein", "receiver", "overall_strength")
      require_columns(ccc, needed, paste("SecAct CCC", sample_id))
      ccc_rows <<- rbind(ccc_rows, data.frame(
        record_type = "inferred_edge", sample_id = sample_id, stage = stage,
        sender = ccc$sender, secretedProtein = ccc$secretedProtein, receiver = ccc$receiver,
        overall_strength = ccc$overall_strength, status = "PASS", n_nuclei = NA_integer_, reason = ""
      ))
    } else {
      ccc_rows <<- rbind(ccc_rows, data.frame(
        record_type = "sample_summary", sample_id = sample_id, stage = stage,
        sender = "", secretedProtein = "", receiver = "", overall_strength = NA_real_,
        status = "NEGATIVE_COMPLETE", n_nuclei = ncol(base),
        reason = "no_CCC_edge_after_fixed_native_thresholds"
      ))
    }
    saveRDS(ccc_obj, ccc_path)
    if (!inherits(readRDS(ccc_path), "Seurat")) stop("RunSecActCCC object failed reload: ", sample_id)
    append_status("RunSecActCCC", sample_id, stage, "snRNA", ccc_status, TRUE,
                  evidence_path = ccc_path)
    if (pilot) append_audit("RunSecActCCC_scRNAseq", sample_id, ccc_status, TRUE,
                            paste0("n_edges=", nrow(ccc)), ccc_path)
    flush_tables()
    message("P06_EVENT module=snRNA sample=", sample_id, " status=PASS ccc_status=", ccc_status)
  }

  visium_file <- function(sample_id) {
    stem <- sub("^MI_D3_", "MI_Day3_", sample_id)
    stem <- sub("^MI_D7_", "MI_Day7_", stem)
    file.path(root, "objects/communication", paste0(stem, "_phase08.rds"))
  }

  build_spacet <- function(sample_id) {
    path <- visium_file(sample_id)
    if (!file.exists(path)) stop("Authoritative Visium object is missing: ", sample_id)
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat") || !"Spatial" %in% Assays(obj) ||
        !"counts" %in% Layers(obj[["Spatial"]])) stop("Visium object lacks Spatial/counts: ", sample_id)
    counts_mouse <- LayerData(obj[["Spatial"]], layer = "counts")
    if (any(counts_mouse@x < 0) || any(abs(counts_mouse@x - round(counts_mouse@x)) > 1e-8)) {
      stop("Visium raw counts are not non-negative integers: ", sample_id)
    }
    counts_human <- map_counts(counts_mouse)
    if (anyDuplicated(rownames(counts_human)) || !identical(colnames(counts_human), colnames(obj))) {
      stop("Visium ortholog conversion changed spot order or produced duplicate genes: ", sample_id)
    }
    image_name <- names(obj@images)
    if (length(image_name) != 1L) stop("Visium image is not unique: ", sample_id)
    raw_dir <- file.path(root, "spatial_raw", image_name, "unpacked")
    positions_path <- list.files(raw_dir, pattern = "^tissue_positions_list\\.csv$",
                                 recursive = TRUE, full.names = TRUE)
    image_path <- list.files(raw_dir, pattern = "^tissue_lowres_image\\.png$",
                             recursive = TRUE, full.names = TRUE)
    if (length(positions_path) != 1L || length(image_path) != 1L) {
      stop("Raw Visium positions/image are not unique: ", sample_id)
    }
    pos <- read.csv(positions_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(pos) <- c("barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres")
    if (anyDuplicated(pos$barcode)) stop("Duplicate raw Visium barcode: ", sample_id)
    pos <- pos[match(colnames(counts_human), pos$barcode), , drop = FALSE]
    if (anyNA(pos$barcode)) stop("Not all modeled spots have raw coordinates: ", sample_id)
    lowres <- obj@images[[1]]@scale.factors$lowres
    coords <- data.frame(
      pixel_row = pos$pxl_row_in_fullres * lowres,
      pixel_col = pos$pxl_col_in_fullres * lowres,
      array_row = pos$array_row, array_col = pos$array_col,
      coordinate_x_um = pos$array_col * 0.5 * 100,
      coordinate_y_um = pos$array_row * 0.5 * sqrt(3) * 100,
      row.names = colnames(counts_human)
    )
    coords$coordinate_y_um <- max(coords$coordinate_y_um) - coords$coordinate_y_um
    raster <- png::readPNG(image_path)
    if (any(coords$pixel_row < 0 | coords$pixel_row > dim(raster)[1]) ||
        any(coords$pixel_col < 0 | coords$pixel_col > dim(raster)[2])) {
      stop("Low-resolution pixel coordinates fall outside the raster: ", sample_id)
    }
    meta <- obj[[]][colnames(counts_human), , drop = FALSE]
    meta$source_species <- "mouse"
    meta$modeling_namespace <- "human_ortholog"
    meta$biological_unit <- "physical_tissue"
    if (!identical(colnames(counts_human), rownames(coords)) ||
        !identical(colnames(counts_human), rownames(meta))) stop("SpaCET input order mismatch: ", sample_id)
    sp <- SpaCET::create.SpaCET.object(
      counts = counts_human, spotCoordinates = coords, metaData = meta,
      imagePath = image_path, platform = "Visium", organism = "human"
    )
    if (!inherits(sp, "SpaCET") || !identical(sp@input$organism, "human")) {
      stop("SpaCET constructor contract failed: ", sample_id)
    }
    list(spacet = sp, n_spots = ncol(counts_human))
  }

  extract_segments <- function(plot) {
    for (layer in plot$layers) {
      d <- layer$data
      if (is.data.frame(d) && all(c("x_start", "y_start", "x_change", "y_change", "x_end", "y_end") %in% colnames(d))) {
        return(d)
      }
    }
    data.frame(x_start = numeric(), y_start = numeric(), x_change = numeric(),
               y_change = numeric(), x_end = numeric(), y_end = numeric())
  }

  render_velocity <- function(native_plot, segments, path, title) {
    if (nrow(segments)) {
      p <- native_plot
    } else {
      p <- ggplot() + annotate("text", x = 0, y = 0,
                               label = "No vectors under fixed parameters") +
        labs(title = title) + theme_void()
    }
    ggsave(path, plot = p, width = 7, height = 6, units = "in", device = "pdf", limitsize = FALSE)
    isTRUE(file.info(path)$size > 500)
  }

  run_st <- function(sample_id, pilot = FALSE) {
    stage <- if (startsWith(sample_id, "MI_D3_")) "Day3" else "Day7"
    built <- build_spacet(sample_id)
    if (pilot && built$n_spots != 2251L) stop("Pilot Visium tissue is not exactly 2251 spots")
    sp <- built$spacet
    object_path <- file.path(out, "objects", paste0(sample_id, "_SpaCET_SecAct.rds"))
    activity_only_path <- file.path(out,"objects",paste0(sample_id,"_SpaCET_activity.rds"))
    sp <- if (file.exists(activity_only_path)) readRDS(activity_only_path) else tryCatch(
      RunSecAct_compat(
        inputProfile = sp, mode = "ST", scale.factor = 1e5, sigMatrix = "SecAct",
        is.filter.sig = FALSE, is.group.sig = TRUE, is.group.cor = 0.9,
        lambda = 5e5, nrand = 1000, cores = 1L, backend = "cpu-pure",
        rng_method = "mt19937"
      ),
      error = function(e) {
        append_status("RunSecAct", sample_id, stage, "Visium", "FAILED", TRUE,
                      reason = conditionMessage(e))
        stop(e)
      }
    )
    act <- extract_activity(sp, "Visium")
    if (!file.exists(activity_only_path)) saveRDS(sp,activity_only_path)
    activity_rows <<- rbind(activity_rows,
                            activity_table(act, sample_id, stage, "Visium", ncol(act$zscore)))
    append_status("RunSecAct", sample_id, stage, "Visium", "PASS", TRUE,
                  evidence_path = object_path)
    if (pilot) append_audit("RunSecAct_ST", sample_id, "PASS", TRUE,
                            paste(dim(act$zscore), collapse = "x"), object_path)

    sp <- tryCatch(
      scop::RunSecActSignalingPattern(sp, scale.factor = 1e5, radius = 200, k = 2:5),
      error = function(e) {
        append_status("RunSecActSignalingPattern", sample_id, stage, "Visium", "FAILED", TRUE,
                      reason = conditionMessage(e))
        stop(e)
      }
    )
    pattern <- sp@results$SecAct_output$pattern
    if (!is.list(pattern) || !all(c("ccc.SP", "weight.W", "signal.H") %in% names(pattern))) {
      stop("SecAct pattern structure is incomplete: ", sample_id)
    }
    W <- pattern$weight.W
    H <- pattern$signal.H
    if (!is.matrix(W) || !is.matrix(H) || ncol(W) < 2L || ncol(W) > 5L ||
        any(!is.finite(W)) || any(!is.finite(H))) stop("SecAct pattern matrices failed integrity: ", sample_id)
    k_selected <- ncol(W)
    pattern_rows <<- rbind(pattern_rows, data.frame(
      sample_id = sample_id, stage = stage,
      pattern_index = rep(seq_len(k_selected), each = nrow(W)),
      secretedProtein = rep(rownames(W), times = k_selected),
      pattern_weight = as.vector(W), pattern_scope = "slice_specific_not_cross_slice_aligned"
    ))
    append_status("RunSecActSignalingPattern", sample_id, stage, "Visium", "PASS", TRUE,
                  target = paste0("k=", k_selected), evidence_path = object_path)
    if (pilot) append_audit("RunSecActSignalingPattern", sample_id, "PASS", TRUE,
                            paste0("k_selected=", k_selected, ";W=", paste(dim(W), collapse = "x"),
                                   ";H=", paste(dim(H), collapse = "x")), object_path)

    pg_list <- vector("list", k_selected)
    any_pg <- FALSE
    for (n in seq_len(k_selected)) {
      pg <- tryCatch(
        scop::RunSecActPatternGenes(sp, n = n),
        error = function(e) stop("PatternGenes n=", n, ": ", conditionMessage(e))
      )
      if (!is.matrix(pg) && !is.data.frame(pg)) stop("PatternGenes returned an invalid structure: ", sample_id)
      pg_list[[n]] <- pg
      if (nrow(pg)) {
        any_pg <- TRUE
        pattern_gene_rows <<- rbind(pattern_gene_rows, data.frame(
          sample_id = sample_id, stage = stage, pattern_index = n,
          rank = seq_len(nrow(pg)), secretedProtein = rownames(pg),
          pattern_weight = as.numeric(pg[, n]), status = "PASS"
        ))
      } else {
        pattern_gene_rows <<- rbind(pattern_gene_rows, data.frame(
          sample_id = sample_id, stage = stage, pattern_index = n, rank = NA_integer_,
          secretedProtein = "", pattern_weight = NA_real_, status = "NEGATIVE_COMPLETE"
        ))
      }
    }
    pg_status <- if (any_pg) "PASS" else "NEGATIVE_COMPLETE"
    pg_path <- file.path(out, "objects", paste0(sample_id, "_pattern_genes.rds"))
    saveRDS(pg_list, pg_path)
    if (length(readRDS(pg_path)) != k_selected) stop("PatternGenes object failed reload: ", sample_id)
    append_status("RunSecActPatternGenes", sample_id, stage, "Visium", pg_status, TRUE,
                  target = paste0("all_patterns_1_to_", k_selected), evidence_path = pg_path)
    if (pilot) append_audit("RunSecActPatternGenes", sample_id, pg_status, TRUE,
                            paste0("patterns_called=", k_selected), pg_path)

    saveRDS(sp, object_path)
    sp_reload <- readRDS(object_path)
    if (!inherits(sp_reload, "SpaCET") || is.null(sp_reload@results$SecAct_output$pattern)) {
      stop("SpaCET SecAct object failed reload: ", sample_id)
    }

    for (gene in ligands) {
      axis_gene <- axes[axes$ligand_secact_human == gene, , drop = FALSE]
      for (signal_mode in c("sending", "receiving")) {
        safe_gene <- gsub("[^A-Za-z0-9_.-]", "_", gene)
        plot_path <- file.path(out, "figures", paste0("SECACT_", sample_id, "_", safe_gene,
                                                       "_", signal_mode, "_velocity.pdf"))
        velocity_path <- file.path(out, "objects", paste0(sample_id, "_", safe_gene,
                                                            "_", signal_mode, "_velocity.rds"))
        native_plot <- tryCatch(
          scop::RunSecActVelocity(
            sp_reload, mode = "spotST", scale.factor = 1e5, gene = gene,
            signalMode = signal_mode, radius = 200, contourMap = FALSE, animated = FALSE
          ),
          error = function(e) {
            append_status("RunSecActVelocity", sample_id, stage, "Visium", "FAILED", TRUE,
                          target = gene, signal_mode = signal_mode, reason = conditionMessage(e))
            stop(e)
          }
        )
        segments <- extract_segments(native_plot)
        if (nrow(segments) && any(!is.finite(as.matrix(segments[, c(
          "x_start", "y_start", "x_change", "y_change", "x_end", "y_end"
        )])))) stop("Velocity segment data are non-finite: ", sample_id, " / ", gene, " / ", signal_mode)
        vstatus <- if (nrow(segments)) "PASS" else "NEGATIVE_COMPLETE"
        rendered <- render_velocity(native_plot, segments, plot_path,
                                    paste(gene, signal_mode, sample_id))
        if (!rendered) stop("Velocity PDF did not render: ", plot_path)
        saveRDS(list(
          spots = as.data.frame(native_plot$data), segments = segments, status = vstatus,
          format = "compact_native_display_source_v1",
          sample_id = sample_id, stage = stage, gene = gene, signal_mode = signal_mode,
          radius_um = 200, scale_factor = 1e5,
          interpretation = "ligand_level_computational_vector_not_receptor_specific"
        ), velocity_path)
        vr <- readRDS(velocity_path)
        if (!identical(vr$status, vstatus)) stop("Velocity object failed reload: ", velocity_path)
        if (nrow(segments)) {
          lengths <- sqrt(segments$x_change^2 + segments$y_change^2)
          med_len <- median(lengths)
          iqr_len <- IQR(lengths)
          med_dx <- median(segments$x_change)
          med_dy <- median(segments$y_change)
        } else {
          med_len <- iqr_len <- med_dx <- med_dy <- 0
        }
        velocity_rows <<- rbind(velocity_rows, data.frame(
          sample_id = sample_id, stage = stage, axis_id = axis_gene$axis_id,
          ligand_mouse = axis_gene$ligand_mouse, ligand_secact_human = gene,
          signal_mode = signal_mode, n_vectors = nrow(segments),
          median_vector_length = med_len, IQR_vector_length = iqr_len,
          median_dx = med_dx, median_dy = med_dy, status = vstatus
        ))
        append_status("RunSecActVelocity", sample_id, stage, "Visium", vstatus, TRUE,
                      target = gene, signal_mode = signal_mode, evidence_path = velocity_path)
        if (pilot) append_audit(paste0("RunSecActVelocity_", gene, "_", signal_mode),
                                sample_id, vstatus, TRUE, paste0("n_vectors=", nrow(segments)), plot_path)
      }
    }
    flush_tables()
    message("P06_EVENT module=Visium sample=", sample_id, " status=PASS patterns=", k_selected,
            " ligands=", length(ligands))
  }

  write_phase_status(
    out = out, phase = "P06_SECACT_AMENDED", status = "RUNNING",
    started_at = phase_started, input_hash = NA_character_, parameter_hash = parameter_hash,
    n_expected = 12L, n_completed = 0L, producer_functions = paste(producers, collapse = ";"),
    biological_verdict = "PILOT_RUNNING", blocker = "", output_root = out
  )
  flush_tables()
  write_tsv(pilot_audit, file.path(out, "tables/pilot_output_audit.tsv"))
  write_tsv(statuses, file.path(out, "functions_final_status.tsv"))
  write_session_info(file.path(out, "sessionInfo.txt"))

  run_sn("SRR21777597", pilot = TRUE)
  run_st("MI_D3_rep1", pilot = TRUE)
  required_pilot <- c("RunSecAct_scRNAseq", "RunSecActCCC_scRNAseq", "RunSecAct_ST",
                      "RunSecActSignalingPattern", "RunSecActPatternGenes",
                      as.vector(outer(paste0("RunSecActVelocity_", ligands, "_"),
                                      c("sending", "receiving"), paste0)))
  if (!all(required_pilot %in% pilot_audit$component) ||
      !all(pilot_audit$pass[pilot_audit$component %in% required_pilot])) {
    stop("Pilot output audit did not pass all fixed modules")
  }
  message("P06_EVENT pilot_gate=PASS snRNA=SRR21777597 Visium=MI_D3_rep1")

  for (sample_id in setdiff(sn_samples, "SRR21777597")) run_sn(sample_id, pilot = FALSE)
  visium_samples <- c("MI_D3_rep1", "MI_D3_rep2", "MI_D3_rep3",
                      "MI_D7_rep1", "MI_D7_rep2", "MI_D7_rep3")
  for (sample_id in setdiff(visium_samples, "MI_D3_rep1")) run_st(sample_id, pilot = FALSE)

  activity <- activity_rows
  activity_groups <- split(activity, paste(activity$modality, activity$stage,
                                             activity$celltype_or_scope,
                                             activity$secretedProtein, sep = "||"))
  activity_summary <- do.call(rbind, lapply(activity_groups, function(d) data.frame(
    record_type = "activity_zscore", modality = d$modality[[1]], stage = d$stage[[1]],
    feature = d$secretedProtein[[1]], context = d$celltype_or_scope[[1]],
    n_biological_samples = length(unique(d$sample_id)), median_effect = median(d$zscore),
    IQR_effect = IQR(d$zscore), n_positive = sum(d$zscore > 0), n_negative = sum(d$zscore < 0),
    inference = "descriptive_sample_level_no_spot_or_nucleus_P_value"
  )))
  velocity_groups <- split(velocity_rows, paste(velocity_rows$stage, velocity_rows$axis_id,
                                                 velocity_rows$signal_mode, sep = "||"))
  velocity_summary <- do.call(rbind, lapply(velocity_groups, function(d) data.frame(
    record_type = "velocity_median_vector_length", modality = "Visium", stage = d$stage[[1]],
    feature = d$axis_id[[1]], context = d$signal_mode[[1]],
    n_biological_samples = length(unique(d$sample_id)),
    median_effect = median(d$median_vector_length), IQR_effect = IQR(d$median_vector_length),
    n_positive = sum(d$median_vector_length > 0), n_negative = 0L,
    inference = "descriptive_physical_tissue_level_not_receptor_specific"
  )))
  discovery_summary <- rbind(activity_summary, velocity_summary)
  rownames(discovery_summary) <- NULL
  write_tsv(discovery_summary, file.path(out, "tables/secact_discovery_summary.tsv"))
  flush_tables()

  expected_status <- c(
    RunSecAct = 12L, RunSecActCCC = 6L, RunSecActSignalingPattern = 6L,
    RunSecActPatternGenes = 6L, RunSecActVelocity = 6L * length(ligands) * 2L
  )
  observed_status <- table(factor(statuses[["function"]], levels = names(expected_status)))
  if (!identical(as.integer(observed_status), as.integer(expected_status))) {
    stop("Final function-status cardinality gate failed")
  }
  if (any(statuses$status %in% c("FAILED", "BLOCKED_DEPENDENCY")) ||
      any(!statuses$status %in% c("PASS", "NEGATIVE_COMPLETE", "NOT_APPLICABLE"))) {
    stop("Final function-status terminal-state gate failed")
  }
  if (!all(producers %in% unique(statuses[["function"]][statuses$called]))) {
    stop("Not all five SCOP SecAct producers were truly called")
  }
  required_tables <- file.path(out, c(
    "tables/ortholog_mapping.tsv", "tables/secact_activity.tsv.gz", "tables/secact_ccc.tsv",
    "tables/secact_patterns.tsv", "tables/secact_pattern_genes.tsv",
    "tables/secact_velocity_sample_effects.tsv", "tables/secact_discovery_summary.tsv",
    "tables/pilot_output_audit.tsv", "functions_final_status.tsv"
  ))
  if (any(!file.exists(required_tables) | file.info(required_tables)$size <= 0)) stop("A final P06 table is absent or empty")
  for (p in required_tables) {
    con <- if (grepl("\\.gz$", p)) gzfile(p, open = "rt") else file(p, open = "rt")
    header <- tryCatch(readLines(con, n = 1L, warn = FALSE), finally = close(con))
    if (length(header) != 1L || !nzchar(header)) stop("A final P06 table cannot be reread: ", p)
  }
  pdfs <- list.files(file.path(out, "figures"), pattern = "^SECACT_.*\\.pdf$", full.names = TRUE)
  if (length(pdfs) != 6L * length(ligands) * 2L || any(file.info(pdfs)$size <= 500)) {
    stop("Velocity PDF count/render gate failed")
  }
  ccc_called <- statuses$sample_id[statuses[["function"]] == "RunSecActCCC" & statuses$called]
  object_paths <- c(
    file.path(out, "objects", paste0(sn_samples, "_snRNA_RunSecAct.rds")),
    file.path(out, "objects", paste0(ccc_called, "_snRNA_RunSecActCCC.rds")),
    file.path(out, "objects", paste0(visium_samples, "_SpaCET_SecAct.rds"))
  )
  if (any(!file.exists(object_paths) | file.info(object_paths)$size <= 0)) stop("A required final object is absent")
  for (p in file.path(out, "objects", paste0(visium_samples, "_SpaCET_SecAct.rds"))) {
    if (!inherits(readRDS(p), "SpaCET")) stop("Final SpaCET reload gate failed: ", p)
  }
  write_session_info(file.path(out, "sessionInfo.txt"))
  input_hash <- sha256_value(c(sha256_file(sn_path), sha256_file(axes_path), sha256_file(spacet_status_path)))
  write_phase_status(
    out = out, phase = "P06_SECACT_AMENDED", status = "PASS",
    started_at = phase_started, input_hash = input_hash, parameter_hash = parameter_hash,
    n_expected = 12L, n_completed = 12L, producer_functions = paste(producers, collapse = ";"),
    biological_verdict = "SECACT_FIXED_PILOT_AND_SIX_SAMPLE_RUN_COMPLETE",
    blocker = "", output_root = out
  )
  if (!file.create(file.path(out, "_SUCCESS"))) stop("Could not create P06 _SUCCESS")
  message("P06_EVENT status=PASS samples=12 producers=5 frozen_ligands=", length(ligands))
}

tryCatch(main(), error = function(e) {
  status_path <- file.path(out, "functions_final_status.tsv")
  current <- tryCatch(if (file.exists(status_path)) read_tsv(status_path) else empty_status(),
                      error = function(...) empty_status())
  write_failure(e, current)
  message("P06_EVENT status=FAILED error=", conditionMessage(e))
  quit(save = "no", status = 1L)
})
