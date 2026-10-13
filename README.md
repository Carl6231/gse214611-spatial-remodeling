# GSE214611 spatial remodeling code archive

Available analysis and figure-generation scripts associated with **Spatial remodeling and tissue composition after myocardial infarction in mice**.

## Archive scope

This repository preserves the research scripts available in the local manuscript workspace on 14 October 2026, together with six machine-readable supplementary summary tables. It is a partial archival release. It has not been validated as an end-to-end reproduction of the manuscript.

- `analysis/preprocessing`: spatial-object construction and quality-control scripts.
- `analysis/meringue`: retained MERINGUE execution scripts.
- `analysis/phase10`: spatial-signature, enrichment, cell-source, external-projection and competing-explanation scripts, plus their shared helper and preflight.
- `figures/historical_v2`: an earlier figure-generation source set. The script `render_figure4_v2.R` uses the original 5,937-center reference and corresponds to sensitivity material, not the corrected 5,151-center main Figure 4. Historical panel titles do not supersede the manuscript's interpretation.
- `tables`: Supplementary Tables S1–S6, preserving the precision reported in the manuscript.

## Data and retained inputs

Public source datasets are available from NCBI GEO:

- [GSE214611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214611)
- [GSE280373](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280373), used for a descriptive external comparison.

The scripts preserve their original paths and workflow contracts. Several consume frozen Seurat objects, metadata tables and saved results in the original project directory. Raw sequencing files, Seurat objects and these full intermediate results are not included in this repository. The CSV summary tables do not replace those inputs.

The local archive did not contain `figure_common_v2.R`, the corrected-neighborhood analysis and final Figure 4 source, final Figure 5 source, or the complete CellOracle and secreted-axis execution sources. Some included scripts depend on those missing resources. Re-running them requires the corresponding inputs and dependencies; do not interpret the presence of a plotting script as evidence that every manuscript result can currently be reproduced from this repository alone.

## Software and validation

The manuscript reports R 4.5.1 and SCOP 0.9.0 for later analyses, with earlier outputs originating from SCOP 0.8.9. Individual scripts identify additional packages and fixed settings. Preserve those analysis versions when reconstructing the original environment; this release does not install dependencies.

The deposited R and Python source files were checked for syntactic parseability. Scientific analyses were not rerun for this release. No raw or newly computed scientific results were generated during repository preparation.

## Scientific interpretation

Physical tissues and biological samples remain the relevant summary units. Thousands of spots or nuclei do not supply independent biological replication. Original-reference maps, corrected-reference maps, ligand–receptor candidates and computational perturbation outputs have distinct provenance and should not be combined as interchangeable validation evidence.
