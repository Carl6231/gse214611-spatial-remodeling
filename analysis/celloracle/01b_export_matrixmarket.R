#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

OUT <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/21_scop_expansion/X01_OMICVERSE_CELLORACLE"
INPUT_RDS <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/03_snrna_reference/myeloid_refinement/myeloid_reclustered_annotated.rds"
EXPECTED_BALANCED <- 4659L

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

selected_path <- file.path(OUT, "tables", "selected_cells.tsv")
stopifnot(dir.exists(OUT), file.exists(INPUT_RDS), file.exists(selected_path))
selected <- read.delim(selected_path, check.names = FALSE)
stopifnot(nrow(selected) == EXPECTED_BALANCED, !anyDuplicated(selected$cell_id))

message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " MATRIXMARKET_FALLBACK_START")
obj <- readRDS(INPUT_RDS)
stopifnot(all(selected$cell_id %in% colnames(obj)))
obj <- subset(obj, cells = selected$cell_id)
obj <- obj[, selected$cell_id]
DefaultAssay(obj) <- "decontXcounts"

counts <- LayerData(obj, assay = "decontXcounts", layer = "counts")
norm <- LayerData(obj, assay = "decontXcounts", layer = "data")
stopifnot(identical(colnames(counts), selected$cell_id))
stopifnot(identical(colnames(norm), selected$cell_id))
stopifnot(identical(rownames(counts), rownames(norm)))
stopifnot(all(counts@x >= 0), all(counts@x == floor(counts@x)))
stopifnot(max(abs(Matrix::colSums(expm1(norm)) - 10000)) < 1e-5)

paths <- file.path(OUT, "objects", c("axis_counts.mtx", "axis_normalized.mtx", "axis_genes.tsv", "axis_barcodes.tsv"))
stopifnot(!any(file.exists(paths)))
Matrix::writeMM(counts, paths[[1]])
Matrix::writeMM(norm, paths[[2]])
write.table(data.frame(gene = rownames(counts)), paths[[3]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(cell_id = colnames(counts)), paths[[4]], sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c("R SESSION", capture.output(sessionInfo())), file.path(OUT, "sessionInfo.txt"))
message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " MATRIXMARKET_FALLBACK_COMPLETE cells=", ncol(counts), " genes=", nrow(counts), " nnz_counts=", length(counts@x), " nnz_normalized=", length(norm@x))
