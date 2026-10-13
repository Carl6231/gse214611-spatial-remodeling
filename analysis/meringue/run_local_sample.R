args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, args[[2]] %in% c("check", "run"))
run_dir <- normalizePath(args[[1]], mustWork = TRUE)
mode <- args[[2]]
id <- "MI_D3_rep1"
lib <- file.path(dirname(run_dir), "Rlib")
.libPaths(c(lib, .libPaths()))
options(warn = 1)
stamp <- function(...) cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), ..., "\n")
options(error = function() {
  stamp("FAILED", "mode=", mode, "pid=", Sys.getpid())
  traceback(12)
  quit(save = "no", status = 1L)
})

threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
stopifnot(all(Sys.getenv(threads) == "1"))
packages <- c("Seurat", "Matrix", "scop", "MERINGUE", "geometry", "dynamicTreeCut")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(scop) })
stopifnot(as.character(getRversion()) == "4.5.1",
          as.character(packageVersion("scop")) == "0.9.0",
          as.character(packageVersion("MERINGUE")) == "1.0",
          identical(find.package("scop"), file.path(lib, "scop")),
          identical(find.package("MERINGUE"), file.path(lib, "MERINGUE")))
stopifnot(all(c("RunMERINGUE", "SpatialCoordinates", "SpatialSpotPlot") %in%
                getNamespaceExports("scop")),
          "backend" %in% names(formals(scop::RunMERINGUE)))
stamp("RUNTIME_OK", "pid=", Sys.getpid(), "platform=", R.version$platform)
print(vapply(packages, function(p) as.character(packageVersion(p)), character(1)))

read_tsv <- function(name) read.delim(file.path(run_dir, "inputs", name),
  stringsAsFactors = FALSE, check.names = FALSE, quote = "", na.strings = "")
contract <- read_tsv("discovery_object_contract.tsv")
sample_row <- contract[contract$biological_sample_id == id, , drop = FALSE]
stopifnot(nrow(sample_row) == 1L, sample_row$cohort == "GSE214611",
          sample_row$analysis_role == "primary_replicated",
          sample_row$condition == "MI", sample_row$timepoint == "Day3",
          sample_row$n_spots == 2251L, sample_row$image == "GSM6613084")
parameters <- read_tsv("parameters.tsv")
fixed <- c(scop_version = "0.9.0", MERINGUE_version = "1.0",
  mode = "autocorrelation;cross_correlation;modules", assay = "Spatial", layer = "data",
  coordinate_space = "raw", feature_detection_rate = "0.10",
  feature_min_primary_samples = "3 of 6", nfeatures = "2000", min_spots = "5",
  filterDist = "NA", binary = "TRUE", alternative = "greater", nperm = "0",
  ncores = "1", seed = "11", backend = "cpp")
stopifnot(!anyDuplicated(parameters$parameter),
          all(names(fixed) %in% parameters$parameter),
          identical(unname(fixed), parameters$value[match(names(fixed), parameters$parameter)]))
feature_table <- read_tsv("feature_universe.tsv")
stopifnot(is.logical(feature_table$eligible), !anyNA(feature_table$eligible))
features_use <- feature_table$feature[feature_table$eligible]
stopifnot(length(features_use) == 7844L, !anyNA(features_use),
          !anyDuplicated(features_use), all(nzchar(features_use)))
stamp("FROZEN_INPUTS_OK", "eligible_features=", length(features_use))
print(parameters[parameters$parameter %in% names(fixed), c("parameter", "value")], row.names = FALSE)

input <- file.path(run_dir, "inputs", basename(sample_row$object_file))
stopifnot(file.info(input)$size == 167376724)
stamp("READ_INPUT", input, "server_source=", sample_row$object_file)
x <- readRDS(input)
stopifnot(inherits(x, "Seurat"), ncol(x) == 2251L,
          "Spatial" %in% Assays(x), identical(Images(x), "GSM6613084"),
          "biological_sample_id" %in% colnames(x[[]]),
          identical(unique(as.character(x[[]]$biological_sample_id)), id),
          all(c("counts", "data") %in% Layers(x[["Spatial"]])),
          is.null(x@tools$MERINGUE))
mat <- LayerData(x[["Spatial"]], layer = "data")
stopifnot(inherits(mat, "sparseMatrix"), identical(colnames(mat), colnames(x)),
          all(is.finite(mat@x)), all(features_use %in% rownames(mat)))
coords <- scop::SpatialCoordinates(x, image = Images(x)[[1L]], space = "raw")$data
stopifnot(identical(rownames(coords), colnames(x)),
          all(is.finite(as.matrix(coords[, c("x", "y"), drop = FALSE]))))
stamp("INPUT_OK", "dimensions=", paste(dim(mat), collapse = "x"),
      "object_GiB=", round(as.numeric(object.size(x)) / 1024^3, 3),
      "image=", Images(x), "biological_sample=", id)
cat("METADATA_COLUMNS\t", paste(colnames(x[[]]), collapse = ";"), "\n", sep = "")
cat("EXISTING_TOOLS\t", paste(names(x@tools), collapse = ";"), "\n", sep = "")
print(formals(scop::RunMERINGUE))
rm(mat, feature_table, coords, contract); invisible(gc())
if (mode == "check") {
  stamp("CHECK_COMPLETE", "no_analysis_executed=TRUE")
  quit(save = "no", status = 0L)
}

output_rds <- file.path(run_dir, "objects", paste0(id, "_MERINGUE.rds"))
plot_file <- file.path(run_dir, "figures", paste0("MERINGUE_", id, "_native_spatial_features.pdf"))
table_files <- file.path(run_dir, "tables", paste0(id,
  c("_autocorrelation.tsv", "_cross_correlation_top1000.tsv", "_modules.tsv")))
session_file <- file.path(run_dir, "logs", "sessionInfo.txt")
stopifnot(!any(file.exists(c(output_rds, plot_file, table_files, session_file))))
capture.output(sessionInfo(), file = session_file)
stamp("RUN_START", id, "pid=", Sys.getpid(), "single_discovery_sample_only=TRUE")
flush.console()

# Same complete SCOP call as the server pilot; reuse its frozen feature universe.
timing <- system.time(x <- scop::RunMERINGUE(
  srt = x, assay = "Spatial", layer = "data", image = Images(x)[[1L]],
  features = features_use, mode = c("autocorrelation", "cross_correlation", "modules"),
  nfeatures = 2000L, min_spots = 5L, filterDist = NA_real_, binary = TRUE,
  alternative = "greater", nperm = 0L, ncores = 1L, set_variable_features = FALSE,
  store_results = TRUE, verbose = TRUE, seed = 11L, coordinate_space = "raw", backend = "cpp"
))
stamp("PRODUCER_RETURNED", id, "elapsed_seconds=", timing[["elapsed"]],
      "user_seconds=", timing[["user.self"]], "system_seconds=", timing[["sys.self"]])
tool <- x@tools$MERINGUE
stopifnot(all(c("autocorrelation", "cross_correlation", "modules", "coords",
                "weight", "parameters") %in% names(tool)),
          identical(rownames(tool$coords), colnames(x)),
          all(dim(tool$weight) == ncol(x)),
          is.data.frame(tool$autocorrelation), is.data.frame(tool$cross_correlation),
          is.data.frame(tool$modules),
          nrow(tool$cross_correlation) == choose(length(tool$pairwise_features), 2L))
saveRDS(x, output_rds)
rm(x, tool); invisible(gc())
x <- readRDS(output_rds)
tool <- x@tools$MERINGUE
stopifnot(inherits(x, "Seurat"), ncol(x) == 2251L,
          identical(rownames(tool$coords), colnames(x)))
stamp("OBJECT_SAVED_AND_READABLE", output_rds)

cross <- tool$cross_correlation
cross <- utils::head(cross[order(abs(cross$correlation), decreasing = TRUE), , drop = FALSE], 1000L)
tables <- list(tool$autocorrelation, cross, tool$modules)
for (i in seq_along(tables)) {
  tab <- tables[[i]]
  tab$biological_sample_id <- rep(id, nrow(tab))
  tab$analysis_scope <- rep("within_sample_descriptive", nrow(tab))
  tab$spot_level_p_used_for_group_inference <- rep(FALSE, nrow(tab))
  write.table(tab, table_files[[i]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
auto <- tool$autocorrelation
top <- auto$feature[is.finite(auto$score)]
top <- top[order(auto$score[is.finite(auto$score)], decreasing = TRUE)]
stopifnot(length(top) > 0L)
p <- scop::SpatialSpotPlot(
  srt = x, features = utils::head(top, 4L), assay = "Spatial", layer = "data",
  image = Images(x)[[1L]], overlay_image = TRUE, combine = TRUE, ncol = 2L, verbose = FALSE
)
grDevices::pdf(plot_file, width = 10, height = 8, useDingbats = FALSE)
print(p)
grDevices::dev.off()
stopifnot(all(file.info(c(output_rds, table_files, plot_file))$size > 0))
stamp("SINGLE_SAMPLE_COMPLETE", id, "n_autocorrelation=", nrow(tool$autocorrelation),
      "n_cross_correlation=", nrow(tool$cross_correlation), "n_module_features=", nrow(tool$modules),
      "other_samples_started=FALSE", "full_P03_complete=FALSE", "figure_visual_review_pending=TRUE")
