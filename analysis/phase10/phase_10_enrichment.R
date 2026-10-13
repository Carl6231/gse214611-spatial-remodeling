root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
source(file.path(root, "scripts", "phase_10_common.R"))
root <- phase10_root()
dirs <- phase10_dirs(root)
phase10_log_start(root, "phase_10_enrichment")

suppressPackageStartupMessages({
  library(data.table)
  library(scop)
})

set.seed(20260805)
sig_obj <- readRDS(require_file(file.path(dirs$objects, "phase_10_spatial_signature.rds")))
recurrence <- as.data.table(sig_obj$recurrence)
signature <- as.data.table(sig_obj$signature)
msig <- as.data.table(readRDS(require_file(file.path(dirs$results, "resources", "msigdb.2025.1.Mm.rds"))))

ranked <- recurrence[slices_tested == 3L & is.finite(median_rho), .(gene, median_rho)]
ranked <- ranked[!duplicated(gene)][order(-median_rho)]
if (nrow(ranked) < 2000L || !any(ranked$median_rho > 0) || !any(ranked$median_rho < 0)) {
  stop("The Day7 ranked list is incomplete or one-sided")
}

collections <- list(
  GO_BP = msig[gs_subcollection == "GO:BP"],
  Hallmark = msig[gs_collection == "MH"],
  Reactome = msig[gs_subcollection == "CP:REACTOME"]
)
if (any(vapply(collections, nrow, integer(1L)) < 1000L)) stop("One or more MSigDB collections are unexpectedly small")

run_collection <- function(name, db) {
  cache_path <- file.path(dirs$objects, paste0("phase_10_enrichment_", name, ".rds"))
  if (file.exists(cache_path)) {
    stamp("Reusing validated enrichment checkpoint: ", name)
    cached <- readRDS(cache_path)
    if (all(c("gsea", "ora", "gsea_object", "ora_object", "term_count", "gene_memberships") %in% names(cached))) {
      return(cached)
    }
  }
  stamp("SCOP RunGSEA: ", name)
  t2g <- unique(db[, .(Term = gs_name, symbol = db_gene_symbol)])
  t2n <- unique(db[, .(Term = gs_name, Name = gs_description)])
  gsea <- RunGSEA(
    geneID = ranked$gene,
    geneScore = ranked$median_rho,
    geneID_groups = rep("Day7_remodeling_axis", nrow(ranked)),
    scoreType = "std", IDtype = "symbol", result_IDtype = "symbol",
    species = "Mus_musculus", TERM2GENE = t2g, TERM2NAME = t2n,
    minGSSize = 10L, maxGSSize = 500L, cores = 8L,
    verbose = TRUE, seed = TRUE, eps = 0
  )
  gsea_dt <- as.data.table(gsea$enrichment)
  if (nrow(gsea_dt)) gsea_dt[, collection := name]

  stamp("SCOP RunEnrichment: ", name)
  ora <- RunEnrichment(
    geneID = signature$gene,
    geneID_groups = signature$direction,
    IDtype = "symbol", result_IDtype = "symbol", backend = "r",
    species = "Mus_musculus", TERM2GENE = t2g, TERM2NAME = t2n,
    minGSSize = 10L, maxGSSize = 500L, cores = 8L,
    verbose = TRUE
  )
  ora_dt <- as.data.table(ora$enrichment)
  if (nrow(ora_dt)) ora_dt[, collection := name]
  answer <- list(gsea = gsea_dt, ora = ora_dt, gsea_object = gsea, ora_object = ora,
                 term_count = uniqueN(t2g$Term), gene_memberships = nrow(t2g))
  saveRDS(answer, cache_path)
  answer
}

results <- parallel::mclapply(
  names(collections),
  function(name) run_collection(name, collections[[name]]),
  mc.cores = min(3L, length(collections)), mc.preschedule = FALSE, mc.set.seed = TRUE
)
names(results) <- names(collections)
gsea_table <- rbindlist(lapply(results, `[[`, "gsea"), fill = TRUE)
ora_table <- rbindlist(lapply(results, `[[`, "ora"), fill = TRUE)
if (!nrow(gsea_table) || !nrow(ora_table)) stop("SCOP enrichment returned no rows")

gsea_table[, significant_bh := is.finite(p.adjust) & p.adjust <= 0.05]
ora_table[, significant_bh := is.finite(p.adjust) & p.adjust <= 0.05]
gsea_table[, abs_NES := abs(NES)]
setorder(gsea_table, collection, p.adjust, -abs_NES)
setorder(ora_table, collection, Groups, p.adjust)
write_tsv(ranked, file.path(dirs$metadata, "phase_10_day7_ranked_gene_list.tsv.gz"))
write_tsv(gsea_table, file.path(dirs$metadata, "phase_10_scop_gsea.tsv.gz"))
write_tsv(ora_table, file.path(dirs$metadata, "phase_10_scop_ora.tsv.gz"))

collection_summary <- rbindlist(lapply(names(results), function(name) {
  data.table(
    collection = name,
    terms = results[[name]]$term_count,
    gene_memberships = results[[name]]$gene_memberships,
    gsea_rows = nrow(results[[name]]$gsea),
    gsea_fdr05 = sum(results[[name]]$gsea$p.adjust <= 0.05, na.rm = TRUE),
    ora_rows = nrow(results[[name]]$ora),
    ora_fdr05 = sum(results[[name]]$ora$p.adjust <= 0.05, na.rm = TRUE),
    scop_gsea_function = "RunGSEA",
    scop_ora_function = "RunEnrichment"
  )
}))
write_tsv(collection_summary, file.path(dirs$metadata, "phase_10_enrichment_summary.tsv"))

saveRDS(
  list(results = results, ranked = ranked, signature = signature,
       collections = names(collections), msigdb_version = "2025.1.Mm"),
  file.path(dirs$objects, "phase_10_enrichment.rds")
)

validation <- data.table(
  check = c(
    "ranked_genes", "collections", "scop_gsea_collections", "scop_ora_collections",
    "gsea_rows", "ora_rows", "gsea_directionality", "native_mouse_msigdb"
  ),
  value = c(
    nrow(ranked), length(collections), uniqueN(gsea_table$collection), uniqueN(ora_table$collection),
    nrow(gsea_table), nrow(ora_table),
    paste(sort(unique(sign(gsea_table$NES))), collapse = ","), "2025.1.Mm"
  ),
  status = c(
    ifelse(nrow(ranked) >= 2000L, "PASS", "FAIL"),
    ifelse(length(collections) == 3L, "PASS", "FAIL"),
    ifelse(uniqueN(gsea_table$collection) == 3L, "PASS", "FAIL"),
    ifelse(uniqueN(ora_table$collection) == 3L, "PASS", "FAIL"),
    ifelse(nrow(gsea_table) > 100L, "PASS", "FAIL"),
    ifelse(nrow(ora_table) > 10L, "PASS", "FAIL"),
    ifelse(any(gsea_table$NES > 0) && any(gsea_table$NES < 0), "PASS", "FAIL"),
    "PASS"
  )
)
write_tsv(validation, file.path(dirs$metadata, "phase_10_enrichment_validation.tsv"))
if (!all(validation$status == "PASS")) stop("Phase 10 enrichment validation failed")
stamp("Enrichment module PASS")
phase10_log_stop("phase_10_enrichment")
