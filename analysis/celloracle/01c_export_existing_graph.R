#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

OUT <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/21_scop_expansion/X01_OMICVERSE_CELLORACLE"
INPUT_RDS <- "/home/zhenhanli/myocardial/scop/GSE214611_spatial/03_snrna_reference/myeloid_refinement/myeloid_reclustered_annotated.rds"

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

selected <- read.delim(file.path(OUT, "tables", "selected_cells.tsv"), check.names = FALSE)
stopifnot(nrow(selected) == 4659L, !anyDuplicated(selected$cell_id))
obj <- readRDS(INPUT_RDS)
stopifnot("Harmony_SNN" %in% names(obj@graphs))
graph <- obj@graphs[["Harmony_SNN"]][selected$cell_id, selected$cell_id, drop = FALSE]
stopifnot(nrow(graph) == 4659L, ncol(graph) == 4659L, isTRUE(Matrix::isSymmetric(graph)))
path <- file.path(OUT, "objects", "axis_harmony_snn.mtx")
stopifnot(!file.exists(path))
Matrix::writeMM(graph, path)
cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " EXISTING_HARMONY_SNN_EXPORTED cells=", nrow(graph), " nnz=", length(graph@x), "\n", sep = "")
