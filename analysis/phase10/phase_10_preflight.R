root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1L]] else "/home/zhenhanli/myocardial/scop/GSE214611_spatial"
source(file.path(root, "scripts", "phase_10_common.R"))
root <- phase10_root()
dirs <- phase10_dirs(root)
phase10_log_start(root, "phase_10_preflight")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(scop)
})

stamp("Runtime and function audit")
required_packages <- c(
  "scop", "Seurat", "Matrix", "data.table", "ggplot2", "patchwork",
  "matrixStats", "FNN", "edgeR", "clusterProfiler", "fgsea",
  "AnnotationDbi", "org.Mm.eg.db", "GO.db", "reactome.db"
)
package_audit <- rbindlist(lapply(required_packages, function(pkg) {
  available <- requireNamespace(pkg, quietly = TRUE)
  data.table(
    package = pkg,
    available = available,
    version = if (available) as.character(packageVersion(pkg)) else NA_character_
  )
}))
if (!all(package_audit$available)) {
  stop("Missing Phase 10 packages: ", paste(package_audit$package[!package_audit$available], collapse = ", "))
}
write_tsv(package_audit, file.path(dirs$metadata, "phase_10_runtime_packages.tsv"))

required_scop <- c(
  "RunGSEA", "RunEnrichment", "SpatialSpotPlot", "GroupHeatmap",
  "FeatureStatPlot", "CellStatPlot"
)
scop_exports <- getNamespaceExports("scop")
scop_audit <- data.table(
  function_name = required_scop,
  exported = required_scop %in% scop_exports,
  scop_version = as.character(packageVersion("scop"))
)
if (!all(scop_audit$exported)) {
  stop("Required SCOP exports absent: ", paste(scop_audit$function_name[!scop_audit$exported], collapse = ", "))
}
write_tsv(scop_audit, file.path(dirs$metadata, "phase_10_scop_function_audit.tsv"))

stamp("Input resource audit")
registry <- fread(require_file(file.path(root, "metadata", "spatial_qc_summary.tsv")))
expected <- unique(registry$sample_key)
if (length(expected) != 12L) stop("Expected 12 physical slice keys, found ", length(expected))
phase9_summary <- fread(require_file(file.path(root, "metadata", "phase_09_niche_continuous_slice_summary.tsv")))
if (!setequal(expected, phase9_summary$sample_key)) stop("Phase 9 continuous-axis sample set mismatch")

msig_path <- require_file(file.path(dirs$results, "resources", "msigdb.2025.1.Mm.rds"))
msig_md5 <- unname(tools::md5sum(msig_path))
if (!identical(msig_md5, "e3b2d3a2aa3bd77b0b22fdf185207afd")) stop("MSigDB resource checksum mismatch")

module_audit <- fread(require_file(file.path(root, "metadata", "phase_09_corrected_module_audit.tsv")))
construction_genes <- sort(unique(unlist(strsplit(module_audit$genes, ";", fixed = TRUE))))
construction_genes <- construction_genes[nzchar(construction_genes)]
construction_table <- data.table(
  gene = construction_genes,
  exclusion_reason = "used_in_phase09_axis_construction"
)
write_tsv(construction_table, file.path(dirs$metadata, "phase_10_axis_construction_genes.tsv"))

stamp("Reading and auditing 12 Phase 9 niche objects")
celltypes <- c(
  "Cardiomyocyte", "Endothelial", "Fibroblast", "Macrophage", "Monocyte",
  "Neutrophil", "Dendritic", "B_cell", "T_NK", "Pericyte", "Schwann_cell"
)
object_rows <- vector("list", length(expected))
input_rows <- vector("list", length(expected) + 4L)
for (i in seq_along(expected)) {
  sample_key <- expected[[i]]
  object_path <- require_file(file.path(root, "objects", "phase_09", "niches", paste0(sample_key, "_phase09_niche.rds")))
  obj <- readRDS(object_path)
  if (!inherits(obj, "Seurat")) stop("Not a Seurat object: ", sample_key)
  md <- obj[[]]
  required_meta <- c(
    "P09N_center", "P09N_remodeling_axis", "P09E_total_myeloid",
    paste0("P09E_RCTD_prop_", celltypes)
  )
  missing_meta <- setdiff(required_meta, colnames(md))
  if (length(missing_meta)) stop(sample_key, " missing metadata: ", paste(missing_meta, collapse = ", "))
  coords <- GetTissueCoordinates(obj)
  center <- as.logical(md$P09N_center) & is.finite(md$P09N_remodeling_axis)
  pmat <- as.matrix(md[, paste0("P09E_RCTD_prop_", celltypes), drop = FALSE])
  row_sums <- rowSums(pmat)
  observed_mean <- mean(md$P09N_remodeling_axis[center])
  summary_mean <- phase9_summary$remodeling_axis[match(sample_key, phase9_summary$sample_key)]
  object_rows[[i]] <- data.table(
    sample_key = sample_key,
    object_path = object_path,
    genes = nrow(obj),
    spots = ncol(obj),
    centers = sum(center),
    axis_finite = sum(is.finite(md$P09N_remodeling_axis)),
    axis_mean_object = observed_mean,
    axis_mean_phase09_table = summary_mean,
    axis_mean_abs_difference = abs(observed_mean - summary_mean),
    coordinates = nrow(coords),
    coordinate_barcode_match = setequal(rownames(coords), colnames(obj)),
    images = length(Images(obj)),
    min_rctd_row_sum = min(row_sums),
    max_rctd_row_sum = max(row_sums),
    status = if (
      nrow(coords) == ncol(obj) && setequal(rownames(coords), colnames(obj)) &&
      sum(center) == phase9_summary$centers[match(sample_key, phase9_summary$sample_key)] &&
      abs(observed_mean - summary_mean) < 1e-10 &&
      all(abs(row_sums - 1) < 1e-6)
    ) "PASS" else "FAIL"
  )
  fi <- file.info(object_path)
  input_rows[[i]] <- data.table(
    input = paste0("phase09_niche_", sample_key), path = object_path,
    bytes = fi$size, mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%S %z")
  )
  rm(obj, md, coords, pmat)
  gc()
}
object_audit <- rbindlist(object_rows)
if (!all(object_audit$status == "PASS")) stop("One or more Phase 9 objects failed Phase 10 preflight")
write_tsv(object_audit, file.path(dirs$metadata, "phase_10_input_object_audit.tsv"))

extra_inputs <- c(
  reference = file.path(root, "objects", "phase_09", "reference_expanded_immune.rds"),
  macrophage_reference = file.path(root, "objects", "phase_09", "reference_macrophage_refined.rds"),
  external_expression = file.path(root, "results", "phase_09", "external_validation", "GSE280373", "raw", "GSE280373_expdata_all_timepoints_LAD_Cryo_Sham.rds"),
  original_srt = "/home/zhenhanli/myocardial/scop/srt.rds"
)
for (j in seq_along(extra_inputs)) {
  path <- require_file(extra_inputs[[j]])
  fi <- file.info(path)
  input_rows[[length(expected) + j]] <- data.table(
    input = names(extra_inputs)[[j]], path = path,
    bytes = fi$size, mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%S %z")
  )
}
input_signatures <- rbindlist(input_rows)
write_tsv(input_signatures, file.path(dirs$metadata, "phase_10_input_signatures.tsv"))

resource_manifest <- data.table(
  resource = "MSigDB mouse native",
  version = "2025.1.Mm",
  url = "https://zenodo.org/records/15800824/files/msigdb.2025.1.Mm.rds?download=1",
  md5 = msig_md5,
  bytes = file.info(msig_path)$size,
  local_path = msig_path
)
write_tsv(resource_manifest, file.path(dirs$metadata, "phase_10_resource_manifest.tsv"))

validation <- data.table(
  check = c(
    "phase09_objects", "object_validation", "axis_values_unchanged",
    "scop_exports", "required_packages", "construction_gene_exclusion",
    "msigdb_checksum"
  ),
  value = c(
    nrow(object_audit), sum(object_audit$status == "PASS"),
    max(object_audit$axis_mean_abs_difference), sum(scop_audit$exported),
    sum(package_audit$available), length(construction_genes), msig_md5
  ),
  status = c(
    ifelse(nrow(object_audit) == 12L, "PASS", "FAIL"),
    ifelse(all(object_audit$status == "PASS"), "PASS", "FAIL"),
    ifelse(max(object_audit$axis_mean_abs_difference) < 1e-10, "PASS", "FAIL"),
    ifelse(all(scop_audit$exported), "PASS", "FAIL"),
    ifelse(all(package_audit$available), "PASS", "FAIL"),
    ifelse(length(construction_genes) >= 30L, "PASS", "FAIL"),
    ifelse(msig_md5 == "e3b2d3a2aa3bd77b0b22fdf185207afd", "PASS", "FAIL")
  )
)
write_tsv(validation, file.path(dirs$metadata, "phase_10_preflight_validation.tsv"))
saveRDS(
  list(object_audit = object_audit, construction_genes = construction_genes,
       package_audit = package_audit, scop_audit = scop_audit),
  file.path(dirs$objects, "phase_10_preflight.rds")
)
if (!all(validation$status == "PASS")) stop("Phase 10 preflight validation failed")
stamp("Phase 10 preflight PASS")
phase10_log_stop("phase_10_preflight")
