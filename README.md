# GSE214611 spatial remodeling code archive

Analysis and figure-generation sources associated with **Spatial remodeling and tissue composition after myocardial infarction in mice**.

## Scope and entry points

The archive includes sources recovered from the original analysis server on 14 October 2026, the subsequent Statial/SecAct completion scripts, revised publication plotting sources, and six machine-readable supplementary summary tables. It is a source archive with retained input contracts, not a standalone workflow that reconstructs every result directly from raw GEO downloads.

- `analysis/preprocessing`: spatial-object construction and quality control.
- `analysis/corrected_reference`: corrected reference construction, deconvolution, macrophage mapping and the 5,151-center neighborhood analysis.
- `analysis/meringue`: retained spatial autocorrelation and module execution sources.
- `analysis/phase10`: fixed expression signature, enrichment, cell-source, external projection and covariate challenge.
- `analysis/celloracle`: the original extraction and computational perturbation sources. Existing negative results were retained rather than rerun or replaced.
- `analysis/secreted_axis`: the frozen candidate-axis screen.
- `analysis/completion_20261014`: Statial contextual associations, SecAct inference, project-local compatibility fixes, the actual sample-resume dispatcher, and final result collection.
- `figures/jtm_20261014`: revised main and supplementary figures. `render_corrected_figure4.R` is the corrected primary neighborhood figure; `render_sensitivity_S14.R` renders the earlier 5,937-center analysis as sensitivity material. `finalize_figure5_6.R` separates the expression-axis question from spatial recurrence, communication candidates and CellOracle summaries.
- `figures/historical_v2`: earlier plotting sources, retained for provenance. Their titles and reference definitions do not supersede the revised manuscript.
- `tables`: Supplementary Tables S1–S6. Their numerical precision is retained; thousands separators are removed in the revised table formatting.

## Required data

Public source datasets:

- [GSE214611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214611), the primary single-nucleus and spatial resource.
- [GSE176092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176092), ten frozen wild-type spatial tissues used for descriptive contextual analysis.
- [GSE280373](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280373), the earlier descriptive external expression comparison.

The scripts preserve their original project paths and workflow contracts. They consume frozen Seurat/SpaCET objects, tissue/coordinate contracts, saved signature and communication tables, and other results in the original project directory. Raw sequencing files and full intermediate objects are not bundled. CSV summaries do not replace those inputs. The recovered sources remove several source-code gaps from the earlier archive but do not make this repository an independently validated end-to-end raw-data reproduction.

## Completion workflow

`run_statial_completion.R` accepts the original project root and a new output directory. It retains the 12 fixed directed relationships and 200-µm radius, constructs joint dominant-class-by-region labels for the contextual parent, preserves unassigned spots, and checks that the coordinate conversion and a one-context label equivalence hold. A relative 1e-8 convex-window expansion prevents numerical exclusion of boundary points. It writes all evaluable and unavailable relationships for 22 tissues.

`P06_run_resumable.R` is the full SecAct producer source. `P06_run_sample.R`, `dispatch_remaining_secact.py` and `collect_secact.R` document the actual recovery of missing samples from an already validated pilot and saved successful outputs. The dispatcher is intentionally tied to that saved run, not a generic scheduler. Sample jobs retain single-thread numerical settings; at most four independent processes are scheduled. The collector checks all 12 samples, 102 producer records and 72 distinct tissue–ligand–mode vector outputs before marking the combined output complete.

`secact_compat.R` creates project-local function copies without modifying installed namespaces. It maps SCOP's backend thread argument to SecAct's `ncores` argument and preserves native zero-edge CCC results using typed empty tables and zero-row guards. No scientific cutoff, ligand set or randomization count is relaxed by these fixes. NMF is explicitly attached and uses a sequential backend within each sample, with the native seed and rank candidates retained.

The native SecAct visualization independently rescales x and y arrow components. Stored arrow lengths, dx and dy therefore describe display geometry, not comparable signaling magnitudes or physical movement. The revised main display uses the fraction of modeled spots with vectors; the complete native display fields remain supplementary outputs.

The collector also exports spatial pattern loadings and the complete native vector display coordinates as numerical figure sources. The revised SecAct figure script reads these TSV sources directly, without requiring the full model objects. Compact vector objects retain the plotting data and segments while omitting the plot environment, which otherwise serializes a copy of the full model for each field. Original model objects are retained on the analysis server; they are not included in the public code archive.

## Software and interpretation

The original outputs used SCOP 0.8.9 and 0.9.0. The completion and revised plotting sources use R 4.5.1, SCOP 0.9.1, Statial 1.12.0, SecAct 1.1.0 and SpaCET 1.4.0. SpaCET was installed in an isolated project library at a fixed commit; its installation source is retained. The historical perturbation environment used OmicVerse 2.3.1 and CellOracle 0.20.0. Individual source files and the saved analysis session information define other dependencies.

Physical tissues and single-nucleus biological samples remain the summary units. Spots, overlapping neighborhoods, nuclei and duplicated ligand-to-receptor mappings are not independent biological replicates. Dominant-mixture contextual relationships, original-label communication candidates, corrected-reference neighborhoods and simulated perturbations retain distinct provenance. Unavailable relationships are not zero effects; completed computations do not establish causality or experimental validation.

## Final figure layout

`figures/jtm_final_20261014` contains the final R/SCOP rendering sources for Figures 1–8 and Supplementary Figures S1–S18. Related views share grouped panels and legends, with a consistent blue, turquoise, coral and violet palette. Figure S18 now spans six pages to show all 72 native vector fields at a readable size. These scripts reuse completed scientific results; no analysis thresholds, candidate sets or statistical values were changed for the layout revision. The frozen canonical-marker heatmap preserves its original color mapping and matching scale. Earlier plotting directories remain historical versions.
