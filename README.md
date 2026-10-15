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

## Final main and supplementary figure structure

`figures/jtm_final_20261014` contains the current R/SCOP rendering sources for Figures 1–8 and Supplementary Figures S1–S18. Related views share grouped panels and legends with a consistent blue, turquoise, coral and violet palette. The scripts reuse completed scientific results; scientific thresholds, candidate sets, input objects and statistical values are preserved. The frozen canonical-marker heatmap retains its original color mapping and scale.

Main figures present the study sequence; supplementary figures provide complete tissue coverage, model sensitivity and evidence limits:

| Main figure | Supplementary support |
| --- | --- |
| 1: tissue registry | S1 complete images; S2 QC and compatibility |
| 2: initial deconvolution | S3 all tissues, both methods and confidence |
| 3: corrected reference | S6 lineage/marker/sample support plus complete 11-class maps; S7 fine-state uncertainty |
| 4: corrected niches | S4–S5 spatial structure/features; S8 corrected primary atlas; S14 original-reference sensitivity |
| 5: expression and composition | S9 models A–D, cellular sources, axis decomposition, fixed genes and descriptive external challenge |
| 6: recurrence and candidates | S10 spatial genes; S11 network/module limits; S12 candidate evidence; S13 all perturbations |
| 7: contextual associations | S15 formal, descriptive/auxiliary and external tissues in separate page groups |
| 8: SecAct model outputs | S16 activity; S17 tissue-local patterns; S18 paired sending/receiving for each ligand |

Figure 5e directly displays saved joint-model R-squared, spatial-block cross-validated R-squared and residual variance; the former cellular-source and axis-decomposition panels are preserved in S9. S6 appends three pages containing RCTD, SpatialDWLS and mapping-confidence views for all 12 tissues. S8 now displays the corrected 5,151-center primary analysis; the earlier 5,937-center atlas is retained together with the prior main display in S14. S18 uses one ligand per page, with six sending fields above the corresponding six receiving fields, retaining all 72 fields.

Current sources use the rendering root `22_JTM_completion_20261014/figures_reviewed`; unchanged completed figures can be reused from the preceding render. `render_support_atlases.R` adds S6 maps and renders S8; `render_original_sensitivity.R` renders the three atlas/summary pages appended to the prior S14 display. `finalize_figure5_6.R` generates the new Figure 5 and the model/source page placed before the three S9 spatial pages. Earlier rendering directories remain historical versions. Full scientific input objects are required and are not included in this source archive.

The final consistency pass uses corrected S054 center accounting in Figure 1h, preserves the initial six-class reference in Figures 1f–g and 2, and orders Figure 7d consistently with S15. Numerically constant CellOracle transition summaries are displayed with a neutral color rather than amplifying floating-point noise during within-metric standardization. The corrected niche atlas uses the same uncertainty color as Figure 4; Day 7 replicate colors are consistent across the S9 covariate displays.

## SpatialEcoTyper reference projection (15 October 2026)

`analysis/spatialecotyper_20261015/run_spatialecotyper.R` was executed through SCOP 0.9.1 and SpatialEcoTyper 1.0.3 for all 12 GSE214611 Visium tissues (25,548 spots). `figures/spatialecotyper_20261015/plot_spatialecotyper.R` renders the overview and complete 120-view atlas used in Supplementary Figures S19–S20. The final captioned PDF assembly preserves the numerical figures and adds supplementary numbering and legends.

This is an exploratory mouse-to-human projection onto the packaged tumor reference, not de novo cardiac ecotype discovery. The inference uses raw Spatial counts from the existing phase_09 tissue objects, SecAct 1.1.0 mouse-to-human mapping, log1p CP10K normalization on the full mapped transcriptome, native within-tissue scaling, and the fixed Bulk_SE_Recovery_W matrix. It retains native seed 39, fitting/stopping rules, batches of 500 spots, two workers and sum-to-one coefficients. Each tissue maps 462 of 499 reference genes; 401–448 have nonzero counts. All spots are retained, and finite/non-negative coefficients and row sums are checked. Reference labels are not cardiac cell identities or measured cell fractions.

The scripts accept the original project root and an existing output directory:

```sh
mkdir -p "$PROJECT_OUTPUT"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript analysis/spatialecotyper_20261015/run_spatialecotyper.R "$PROJECT_ROOT" "$PROJECT_OUTPUT"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript figures/spatialecotyper_20261015/plot_spatialecotyper.R "$PROJECT_ROOT" "$PROJECT_OUTPUT"
```

The input contract is `objects/phase_09/deconvolution/<sample>_phase09_expanded.rds` under the project root. Full input expression objects and the package's reference matrix are not bundled in this repository. The figure script also reads the input tissue objects; standalone output tables alone do not replace those dependencies. It reuses saved coefficients and does not refit the model. The largest displayed component is a relative maximum, without a confidence threshold; continuous maps share a 0–1 scale. Spatial component views use native SCOP plotting and the categorical overview uses ggplot2 to avoid duplicated legends.

The software reference is fixed at SpatialEcoTyper commit `64c52f3fb79efcbf1bf8aea7d80d7d0b1bb505c4`. Supplementary Additional file 16 contains numerical outputs, coordinates, coverage, execution log and session information. Tissue means and medians are descriptive; neither reference coverage nor computational completion validates tumor-to-myocardium transfer. Existing rule-defined cardiac neighborhoods remain separate.

## Figure refinement (15 October 2026)

`figures/refinement_20261015` contains the final display-only pass. `refine_vector_typography.R` strengthens small labels in the existing SVG figures while retaining their text contents, explicit text widths and scientific graphical elements. It consumes the assembled project figure folders; it is not an analysis producer. `refine_model_atlases.R` uses saved SecAct loadings and SpatialEcoTyper RDS/tissue summaries. S17 has fixed-size map slots and a shared within-factor loading key. S19 separates formal-comparison tissues from descriptive/auxiliary tissues. S20 is organized by reference component (10 pages, 12 tissues per page), retains all 120 spatial views and adds the already-saved tissue medians and interquartile ranges. IQR describes within-tissue spot heterogeneity, not biological-replicate uncertainty. Reference coefficients and the shared 0–1 spatial color scale are unchanged.

Run the typography script first, then the model-atlas script, then `assemble_supplement.R`; each accepts the saved project root and/or destination as specified in its argument block. The final supplement has S1–S20 in 78 pages including the guide and six tables. Typography rendering requires R xml2/rsvg/qpdf; model-atlas rendering uses R ggplot2/patchwork and saved numerical outputs. No scientific model is refitted by these scripts. The manuscript/package assembly uses the resulting figures without replacing primary scientific conclusions.
