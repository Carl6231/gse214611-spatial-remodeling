#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

OUT <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/21_scop_expansion/X01_OMICVERSE_CELLORACLE"
INPUT_RDS <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/03_snrna_reference/myeloid_refinement/myeloid_reclustered_annotated.rds"
PHASE13 <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/results/phase13/recurrent_program_genes.tsv"
EXPECTED_RAW <- 4897L
EXPECTED_BALANCED <- 4659L
EXPECTED_MONOCYTE <- 2305L
EXPECTED_SPP1 <- 2354L
SEED <- 20260903L
CAP <- 300L

stopifnot(dir.exists(OUT), file.exists(INPUT_RDS), file.exists(PHASE13))
stopifnot(startsWith(normalizePath(OUT), "/home/zhenhanli/myocardial/scop/GSE214611_spatial/21_scop_expansion/"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(Matrix)
})

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " INPUT_READ_START")
obj <- readRDS(INPUT_RDS)
stopifnot(inherits(obj, "Seurat"))
stopifnot("decontXcounts" %in% Assays(obj))
stopifnot(all(c("sample_id", "timepoint", "myeloid_state") %in% colnames(obj[[]])))
stopifnot(all(c("Harmony", "HarmonyUMAP2D") %in% names(obj@reductions)))

DefaultAssay(obj) <- "decontXcounts"
stored_hvgs <- VariableFeatures(obj, assay = "decontXcounts")
stopifnot(length(stored_hvgs) == 2000L, !anyDuplicated(stored_hvgs))

target_states <- c("Monocyte_Ccr2", "Macrophage_Spp1")
raw_cells <- colnames(obj)[obj$myeloid_state %in% target_states & !is.na(obj$sample_id)]
stopifnot(length(raw_cells) == EXPECTED_RAW)
raw_meta <- obj[[]][raw_cells, c("sample_id", "timepoint", "myeloid_state"), drop = FALSE]
stopifnot(length(unique(raw_meta$sample_id)) == 16L)

set.seed(SEED)
split_cells <- split(rownames(raw_meta), interaction(raw_meta$sample_id, raw_meta$myeloid_state, drop = TRUE, lex.order = TRUE))
chosen <- unlist(lapply(sort(names(split_cells)), function(group_name) {
  cells <- sort(split_cells[[group_name]])
  if (length(cells) <= CAP) cells else sample(cells, CAP, replace = FALSE)
}), use.names = FALSE)
keep_cells <- raw_cells[raw_cells %in% chosen]
stopifnot(length(keep_cells) == EXPECTED_BALANCED)

obj <- subset(obj, cells = keep_cells)
DefaultAssay(obj) <- "decontXcounts"
obj$cell_id <- colnames(obj)

state_counts <- table(obj$myeloid_state)
stopifnot(unname(state_counts[["Monocyte_Ccr2"]]) == EXPECTED_MONOCYTE)
stopifnot(unname(state_counts[["Macrophage_Spp1"]]) == EXPECTED_SPP1)
stopifnot(length(unique(obj$sample_id)) == 16L)

counts <- LayerData(obj, assay = "decontXcounts", layer = "counts")
norm <- LayerData(obj, assay = "decontXcounts", layer = "data")
stopifnot(inherits(counts, "sparseMatrix"), inherits(norm, "sparseMatrix"))
stopifnot(all(counts@x >= 0), all(counts@x == floor(counts@x)))
stopifnot(all(is.finite(norm@x)), all(norm@x >= 0))
cp10k <- Matrix::colSums(expm1(norm))
stopifnot(max(abs(cp10k - 10000)) < 1e-5)

meta <- obj[[]][, c("cell_id", "sample_id", "timepoint", "myeloid_state"), drop = FALSE]
meta <- cbind(data.frame(cell_id = rownames(meta), stringsAsFactors = FALSE), meta[, -1, drop = FALSE])
write_tsv(meta, file.path(OUT, "tables", "selected_cells.tsv"))
write_tsv(data.frame(gene = stored_hvgs), file.path(OUT, "tables", "stored_variable_features.tsv"))

harmony <- Embeddings(obj, reduction = "Harmony")
umap <- Embeddings(obj, reduction = "HarmonyUMAP2D")
stopifnot(identical(rownames(harmony), colnames(obj)), identical(rownames(umap), colnames(obj)))
write_tsv(cbind(data.frame(cell_id = rownames(harmony)), as.data.frame(harmony)), file.path(OUT, "tables", "harmony_embedding.tsv"))
write_tsv(cbind(data.frame(cell_id = rownames(umap)), as.data.frame(umap)), file.path(OUT, "tables", "harmony_umap.tsv"))

sample_state <- as.data.frame(table(sample_id = obj$sample_id, myeloid_state = obj$myeloid_state), stringsAsFactors = FALSE)
sample_state <- sample_state[sample_state$Freq > 0, , drop = FALSE]
write_tsv(sample_state, file.path(OUT, "tables", "sample_state_counts.tsv"))

contract <- data.frame(
  field = c("source_rds", "phase13_programs", "raw_target_cells", "balanced_cells", "biological_samples", "Monocyte_Ccr2", "Macrophage_Spp1", "genes", "stored_variable_features", "normalization", "sample_state_cap", "seed"),
  value = c(INPUT_RDS, PHASE13, EXPECTED_RAW, ncol(obj), length(unique(obj$sample_id)), state_counts[["Monocyte_Ccr2"]], state_counts[["Macrophage_Spp1"]], nrow(obj), length(stored_hvgs), "log1p(CP10K=10000)", CAP, SEED),
  stringsAsFactors = FALSE
)
write_tsv(contract, file.path(OUT, "tables", "input_contract.tsv"))

h5s <- file.path(OUT, "objects", "monocyte_spp1_axis.h5Seurat")
h5ad <- file.path(OUT, "objects", "monocyte_spp1_axis.h5ad")
stopifnot(!file.exists(h5s), !file.exists(h5ad))
message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " H5SEURAT_WRITE_START")
SaveH5Seurat(obj, filename = h5s, overwrite = FALSE, verbose = TRUE)
message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " H5AD_CONVERT_START")
Convert(h5s, dest = "h5ad", assay = "decontXcounts", overwrite = FALSE, verbose = TRUE)
stopifnot(file.exists(h5ad), file.info(h5ad)$size > 0)

writeLines(c("R SESSION", capture.output(sessionInfo())), file.path(OUT, "sessionInfo.txt"))
message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " INPUT_EXTRACTION_COMPLETE cells=", ncol(obj), " genes=", nrow(obj), " samples=", length(unique(obj$sample_id)))
