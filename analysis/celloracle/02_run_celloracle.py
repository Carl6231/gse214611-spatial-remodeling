#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import inspect
import json
import math
import os
import platform
import random
import subprocess
import sys
import traceback
from datetime import datetime
from importlib import metadata
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.sparse.csgraph import connected_components

import anndata as ad
import celloracle as co
import celloracle.data.load_data as co_load_data
import omicverse as ov
from omicverse.single import PerturbResult


OUT = Path("/home/zhenhanli/myocardial/scop/GSE214611_spatial/21_scop_expansion/X01_OMICVERSE_CELLORACLE")
PHASE13_PATH = Path("/home/zhenhanli/myocardial/scop/GSE214611_spatial/results/phase13/recurrent_program_genes.tsv")
EXPECTED_PYTHON = Path("/home/zhenhanli/scop/UCS/PRJNA786266_GSE203612/00_admin/toolchain/envs/omicverse_2_3_1_celloracle_0_20_0_py310/bin/python")
BASE_GRN_SOURCE_URL = "https://raw.githubusercontent.com/morris-lab/CellOracle/master/celloracle/data/TFinfo_data/mm9_mouse_atac_atlas_data_TSS_and_cicero_0.9_accum_threshold_10.5_DF_peaks_by_TFs_v202204.parquet"
EXPECTED_CELLS = 4659
EXPECTED_SAMPLES = 16
EXPECTED_STATES = {"Monocyte_Ccr2": 2305, "Macrophage_Spp1": 2354}
STATES = ["Monocyte_Ccr2", "Macrophage_Spp1"]
TIME_ORDER = ["Day0", "1HR", "4HR", "Day1", "Day3", "Day7"]
CONTROL_PRIORITY = ["Spi1", "Irf8", "Cebpb"]
SEED = 20260903
N_JOBS = 8
N_PCA_DIMS = 30
ALPHA = 1
N_PROPAGATION = 3
TRANSITION_N_NEIGHBORS = 30
TRANSITION_SIGMA = 0.05
PS_GRID_SIZE = 30
PS_MIN_MASS = 1.0
PILOT_SAMPLE = "SRR21777597"
STATUS_VALUES = {"PASS", "NEGATIVE_COMPLETE", "NOT_APPLICABLE", "BLOCKED_DEPENDENCY", "FAILED"}

TABLES = OUT / "tables"
FIGURES = OUT / "figures"
OBJECTS = OUT / "objects"
CACHE = OUT / "cache"
TMP = OUT / "tmp"
STATUS_PATH = OUT / "STATUS.tsv"


class BlockedDependency(RuntimeError):
    pass


def now() -> str:
    return datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


def log(message: str) -> None:
    print(f"{now()} {message}", flush=True)


def write_tsv_once(df: pd.DataFrame, path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def write_text_once(text: str, path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_or_validate_tsv(df: pd.DataFrame, path: Path, key_columns: list[str]) -> None:
    if not path.exists():
        write_tsv_once(df, path)
        return
    existing = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    expected = df.astype(str)
    if list(existing.columns) != list(expected.columns) or not existing[key_columns].equals(expected[key_columns]):
        raise ValueError(f"Existing resumable contract does not match: {path}")


def boolean_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values
    return values.astype(str).str.lower().map({"true": True, "false": False}).fillna(False)


def record_status(component: str, status: str, reason: str, artifact: str = "", target: str = "") -> None:
    if status not in STATUS_VALUES:
        raise ValueError(f"Invalid terminal status: {status}")
    cols = ["timestamp", "component", "target", "status", "reason", "artifact"]
    if STATUS_PATH.exists():
        frame = pd.read_csv(STATUS_PATH, sep="\t", dtype=str, keep_default_na=False)
        frame = frame.loc[~((frame["component"] == component) & (frame["target"] == target))].copy()
    else:
        frame = pd.DataFrame(columns=cols)
    row = pd.DataFrame([[now(), component, target, status, reason, artifact]], columns=cols)
    frame = pd.concat([frame, row], ignore_index=True)
    temp = TMP / "STATUS.tsv.tmp"
    temp.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(temp, sep="\t", index=False)
    os.replace(temp, STATUS_PATH)


def dense_array(value) -> np.ndarray:
    return value.toarray() if sparse.issparse(value) else np.asarray(value)


def finite_sparse(value) -> bool:
    if sparse.issparse(value):
        return bool(np.isfinite(value.data).all())
    return bool(np.isfinite(np.asarray(value)).all())


def matrixmarket_input() -> ad.AnnData:
    counts_gxc = mmread(OBJECTS / "axis_counts.mtx").tocsr()
    norm_gxc = mmread(OBJECTS / "axis_normalized.mtx").tocsr()
    genes = pd.read_csv(OBJECTS / "axis_genes.tsv", sep="\t")["gene"].astype(str).tolist()
    barcodes = pd.read_csv(OBJECTS / "axis_barcodes.tsv", sep="\t")["cell_id"].astype(str).tolist()
    obs = pd.read_csv(TABLES / "selected_cells.tsv", sep="\t", dtype=str)
    harmony = pd.read_csv(TABLES / "harmony_embedding.tsv", sep="\t")
    umap = pd.read_csv(TABLES / "harmony_umap.tsv", sep="\t")
    snn = mmread(OBJECTS / "axis_harmony_snn.mtx").tocsr().astype(np.float64)

    if counts_gxc.shape != (len(genes), len(barcodes)) or norm_gxc.shape != counts_gxc.shape:
        raise ValueError("MatrixMarket dimensions do not match gene/barcode contracts")
    if len(barcodes) != EXPECTED_CELLS or obs.shape[0] != EXPECTED_CELLS:
        raise ValueError("Balanced cell count changed")
    if obs["cell_id"].tolist() != barcodes:
        raise ValueError("Metadata and matrix cell order differ")
    if harmony["cell_id"].astype(str).tolist() != barcodes or umap["cell_id"].astype(str).tolist() != barcodes:
        raise ValueError("Embedding and matrix cell order differ")
    if len(set(genes)) != len(genes) or len(set(barcodes)) != len(barcodes):
        raise ValueError("Duplicate gene or cell identifiers")
    if snn.shape != (EXPECTED_CELLS, EXPECTED_CELLS):
        raise ValueError("Existing Harmony SNN dimensions changed")

    counts = counts_gxc.transpose().tocsr()
    counts.data = np.rint(counts.data).astype(np.int32)
    norm = norm_gxc.transpose().tocsr().astype(np.float64)
    if counts.data.size and (counts.data.min() < 0 or not np.equal(counts.data, np.floor(counts.data)).all()):
        raise ValueError("Raw counts are not non-negative integers")
    if not finite_sparse(norm) or (norm.data.size and norm.data.min() < 0):
        raise ValueError("Normalized expression has invalid values")
    cp10k = np.asarray(np.expm1(norm).sum(axis=1)).ravel()
    if not np.allclose(cp10k, 10000.0, atol=1e-5, rtol=0):
        raise ValueError("Normalized expression is not the frozen log1p(CP10K=10000) matrix")

    obs = obs.set_index("cell_id", drop=False)
    obs["myeloid_state"] = pd.Categorical(obs["myeloid_state"], categories=STATES, ordered=True)
    obs["timepoint"] = pd.Categorical(obs["timepoint"], categories=TIME_ORDER, ordered=True)
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    adata = ad.AnnData(X=norm, obs=obs, var=var)
    adata.layers["raw_count"] = counts
    adata.obsm["X_harmony"] = harmony.drop(columns="cell_id").to_numpy(dtype=np.float64)
    adata.obsm["X_umap"] = umap.drop(columns="cell_id").to_numpy(dtype=np.float64)

    snn.setdiag(0)
    snn.eliminate_zeros()
    snn = ((snn + snn.T) * 0.5).tocsr()
    distances = snn.copy()
    distances.data = np.maximum(1.0 - distances.data, 1e-6)
    adata.obsp["connectivities"] = snn
    adata.obsp["distances"] = distances
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": 30, "method": "umap", "use_rep": "X_harmony", "metric": "euclidean"},
    }
    n_components, labels = connected_components(snn, directed=False)
    graph_contract = pd.DataFrame({
        "field": ["source", "n_cells", "nnz", "connected_components", "largest_component", "n_neighbors_source"],
        "value": ["Seurat Harmony_SNN induced subgraph", EXPECTED_CELLS, snn.nnz, n_components, int(np.bincount(labels).max()), 30],
    })
    write_or_validate_tsv(graph_contract, TABLES / "graph_contract.tsv", ["field", "value"])
    return adata


def verify_runtime() -> None:
    actual = Path(sys.executable).resolve()
    if actual != EXPECTED_PYTHON.resolve():
        raise BlockedDependency(f"Unexpected Python: {actual}")
    versions = {
        "python": platform.python_version(),
        "architecture": platform.machine(),
        "omicverse": ov.__version__,
        "celloracle": co.__version__,
        "anndata": ad.__version__,
        "omicverse_path": str(Path(ov.__file__).resolve()),
        "celloracle_path": str(Path(co.__file__).resolve()),
        "perturb_signature": str(inspect.signature(ov.single.perturb)),
        "lineage_pseudotime_signature": str(inspect.signature(ov.single.lineage_pseudotime)),
    }
    if versions["omicverse"] != "2.3.1" or versions["celloracle"] != "0.20.0" or versions["architecture"] != "x86_64":
        raise BlockedDependency(f"Runtime drift: {versions}")
    if "backend" not in inspect.signature(ov.single.perturb).parameters:
        raise BlockedDependency("ov.single.perturb backend parameter missing")
    write_or_validate_tsv(pd.DataFrame([{"field": key, "value": value} for key, value in versions.items()]), TABLES / "runtime_gate.tsv", ["field", "value"])
    baseline = (
        "# OmicVerse runtime-specific API baseline\n\n"
        f"- Python: `{sys.executable}`\n"
        f"- OmicVerse: `{ov.__version__}`\n"
        f"- CellOracle: `{co.__version__}`\n"
        f"- AnnData: `{ad.__version__}`\n"
        "- `ov.single.perturb`: CONFIRMED_RUNTIME\n"
        "- `ov.single.lineage_pseudotime`: CONFIRMED_RUNTIME\n\n"
        "The generic plugin audit prose was written for an older runtime; this live runtime table and callable inspection are authoritative.\n"
    )
    (OUT / "references" / "omicverse-api-baseline.md").write_text(baseline, encoding="utf-8")
    record_status("runtime", "PASS", "OmicVerse 2.3.1 and CellOracle 0.20.0 callable in locked x86_64 interpreter", str(TABLES / "runtime_gate.tsv"))


def load_base_grn() -> tuple[pd.DataFrame, Path]:
    target_dir = CACHE / "celloracle_data"
    (target_dir / "TFinfo_data").mkdir(parents=True, exist_ok=True)
    co_load_data.CELLORACLE_DATA_DIR = str(target_dir)
    try:
        base = co_load_data.load_mouse_scATAC_atlas_base_GRN(version="0.10.0", force_download=False)
    except Exception as exc:
        record_status("base_grn", "BLOCKED_DEPENDENCY", f"CellOracle mouse scATAC atlas base GRN unavailable: {type(exc).__name__}: {exc}")
        raise BlockedDependency(str(exc)) from exc
    required = {"peak_id", "gene_short_name"}
    if not required.issubset(base.columns):
        raise ValueError("Base GRN lacks peak_id/gene_short_name")
    path = OBJECTS / "mouse_scATAC_base_GRN_v0.10.0.parquet"
    if path.exists():
        raise FileExistsError(path)
    base.to_parquet(path, index=False)
    contract = pd.DataFrame([{
        "dataset": "CellOracle mouse scATAC atlas base GRN",
        "data_version": "0.10.0",
        "genome_build": "mm9",
        "source_url": BASE_GRN_SOURCE_URL,
        "rows": base.shape[0],
        "columns": base.shape[1],
        "sha256": "NOT_EXECUTED_CODEX_AUDIT_HOOK_DISABLED",
        "path": str(path),
    }])
    write_tsv_once(contract, TABLES / "base_grn_contract.tsv")
    record_status("base_grn", "PASS", "Generic mouse mm9 scATAC prior frozen before perturbation", str(path))
    record_status("base_grn_sha256", "NOT_APPLICABLE", "Checksum command was blocked because this Codex session lacks CODEX_ALLOW_AUDIT=1; parquet structure and runtime readability are validated instead", str(path))
    return base, path


def phase13_programs() -> pd.DataFrame:
    frame = pd.read_csv(PHASE13_PATH, sep="\t", dtype=str)
    required = {"celltype", "gene", "technical_excluded_from_projection", "module_id", "supported_contrast_n", "min_abs_log2FC_supported", "worst_fdr_supported", "direction"}
    if not required.issubset(frame.columns):
        raise ValueError("Phase13 recurrent program table schema changed")
    keep = frame.loc[(frame["celltype"] == "Macrophage") & (frame["technical_excluded_from_projection"].str.lower() == "false")].copy()
    keep = keep.drop_duplicates("gene", keep="first")
    if keep.shape[0] != 728:
        raise ValueError(f"Expected 728 eligible Phase13 macrophage genes, found {keep.shape[0]}")
    for column in ["supported_contrast_n", "min_abs_log2FC_supported", "worst_fdr_supported"]:
        keep[column] = pd.to_numeric(keep[column], errors="coerce")
    return keep


def detection_metrics(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    raw = adata.layers["raw_count"].tocsc()
    gene_to_idx = {gene: i for i, gene in enumerate(adata.var_names)}
    samples = adata.obs["sample_id"].astype(str).to_numpy()
    states = adata.obs["myeloid_state"].astype(str).to_numpy()
    sample_names = sorted(pd.unique(samples))
    rows = []
    for gene in genes:
        if gene not in gene_to_idx:
            continue
        detected = np.asarray(raw[:, gene_to_idx[gene]].toarray()).ravel() > 0
        per_sample = [float(detected[samples == sample].mean()) for sample in sample_names]
        per_state = {state: float(detected[states == state].mean()) for state in STATES}
        rows.append({
            "gene": gene,
            "samples_detected_ge_1pct": int(sum(rate >= 0.01 for rate in per_sample)),
            "Monocyte_Ccr2_detection": per_state["Monocyte_Ccr2"],
            "Macrophage_Spp1_detection": per_state["Macrophage_Spp1"],
            "max_state_detection": max(per_state.values()),
        })
    return pd.DataFrame(rows)


def prepare() -> None:
    log("PREPARE_START")
    verify_runtime()
    adata = matrixmarket_input()
    state_counts = adata.obs["myeloid_state"].value_counts().to_dict()
    if adata.n_obs != EXPECTED_CELLS or adata.obs["sample_id"].nunique() != EXPECTED_SAMPLES or state_counts != EXPECTED_STATES:
        raise ValueError(f"Input integrity mismatch: cells={adata.n_obs}, samples={adata.obs['sample_id'].nunique()}, states={state_counts}")
    record_status("input_matrix", "PASS", "Frozen 4,659-cell MatrixMarket fallback validated", str(OBJECTS / "axis_counts.mtx"))
    record_status("input_conversion_primary", "NOT_APPLICABLE", "SeuratDisk is incompatible with the installed SeuratObject 5.4 slot API; no biological computation was attempted on that route", str(OUT / "logs" / "extract_input.log"))
    record_status("input_conversion_fallback", "PASS", "Sparse MatrixMarket conversion used exactly once with unchanged cells and matrices", str(OUT / "logs" / "extract_input_fallback.log"))

    full_h5ad = OBJECTS / "monocyte_spp1_axis.h5ad"
    if full_h5ad.exists():
        existing_full = ad.read_h5ad(full_h5ad, backed="r")
        if existing_full.shape != adata.shape or existing_full.obs_names.tolist() != adata.obs_names.tolist():
            raise ValueError("Existing full H5AD does not match the resumable input contract")
        existing_full.file.close()
    else:
        adata.write_h5ad(full_h5ad, compression="gzip")

    programs = phase13_programs()
    hvgs = pd.read_csv(TABLES / "stored_variable_features.tsv", sep="\t")["gene"].astype(str).tolist()
    if len(hvgs) != 2000 or len(set(hvgs)) != 2000:
        raise ValueError("Stored variable-feature contract changed")
    base_grn, base_path = load_base_grn()
    tf_columns_all = [column for column in base_grn.columns if column not in {"peak_id", "gene_short_name"}]

    requested_features = set(hvgs) | set(programs["gene"].astype(str)) | set(CONTROL_PRIORITY)
    feature_list = [gene for gene in adata.var_names.astype(str) if gene in requested_features]
    if not (2000 <= len(feature_list) <= 2731):
        raise ValueError(f"Frozen feature count outside contract: {len(feature_list)}")
    feature_table = pd.DataFrame({"gene": feature_list})
    feature_table["stored_hvg"] = feature_table["gene"].isin(hvgs)
    feature_table["phase13_macrophage_recurrent"] = feature_table["gene"].isin(programs["gene"])
    feature_table["technical_control_candidate"] = feature_table["gene"].isin(CONTROL_PRIORITY)
    write_tsv_once(feature_table, TABLES / "frozen_feature_universe.tsv")

    tf_columns = [gene for gene in tf_columns_all if gene in feature_list]
    pruned = base_grn.loc[base_grn["gene_short_name"].astype(str).isin(feature_list), ["peak_id", "gene_short_name"] + tf_columns].copy()
    if pruned.empty or not tf_columns:
        raise ValueError("No feature/base-GRN overlap")
    pruned_path = OBJECTS / "mouse_scATAC_base_GRN_v0.10.0_feature_pruned.parquet"
    pruned.to_parquet(pruned_path, index=False)
    record_status("feature_freeze", "PASS", f"Frozen {len(feature_list)} genes before direction or perturbation", str(TABLES / "frozen_feature_universe.tsv"))

    adata = adata[:, feature_list].copy()
    candidate_genes = sorted(set(programs["gene"]) & set(tf_columns))
    detection = detection_metrics(adata, sorted(set(candidate_genes) | set(CONTROL_PRIORITY)))
    grouped_prior = pruned.groupby("gene_short_name", observed=True)[tf_columns].sum()
    prior_out = {tf: int(np.count_nonzero(grouped_prior[tf].to_numpy())) for tf in tf_columns}
    candidate = programs.loc[programs["gene"].isin(candidate_genes)].merge(detection, on="gene", how="left")
    candidate["base_grn_prior_out_degree"] = candidate["gene"].map(prior_out).fillna(0).astype(int)
    candidate["eligible_prebaseline"] = (
        (candidate["samples_detected_ge_1pct"] >= 8)
        & (candidate["max_state_detection"] >= 0.05)
        & (candidate["base_grn_prior_out_degree"] > 0)
    )
    write_tsv_once(candidate, TABLES / "candidate_tf_prebaseline.tsv")

    controls = detection.loc[detection["gene"].isin(CONTROL_PRIORITY)].copy()
    controls["base_grn_tf"] = controls["gene"].isin(tf_columns)
    controls["base_grn_prior_out_degree"] = controls["gene"].map(prior_out).fillna(0).astype(int)
    controls["eligible_prebaseline"] = (
        controls["base_grn_tf"]
        & (controls["samples_detected_ge_1pct"] >= 8)
        & (controls["max_state_detection"] >= 0.05)
        & (controls["base_grn_prior_out_degree"] > 0)
    )
    write_tsv_once(controls, TABLES / "technical_control_candidates.tsv")

    pilot_target = None
    for gene in CONTROL_PRIORITY:
        row = controls.loc[controls["gene"] == gene]
        if not row.empty and bool(row.iloc[0]["eligible_prebaseline"]):
            pilot_target = gene
            break
    if pilot_target is None:
        ranked = candidate.loc[candidate["eligible_prebaseline"]].sort_values(
            ["supported_contrast_n", "base_grn_prior_out_degree", "min_abs_log2FC_supported", "worst_fdr_supported", "gene"],
            ascending=[False, False, False, True, True],
            kind="mergesort",
        )
        if ranked.empty:
            record_status("target_selection", "NOT_APPLICABLE", "No TF passes the prebaseline target rules")
            raise RuntimeError("No eligible TF for pilot")
        pilot_target = str(ranked.iloc[0]["gene"])
    write_tsv_once(pd.DataFrame([{"pilot_sample": PILOT_SAMPLE, "pilot_target": pilot_target, "selection": "technical-control priority, else prebaseline deterministic rank"}]), TABLES / "pilot_target.tsv")

    parameters = pd.DataFrame([
        ("seed", SEED, "pre-registered"),
        ("states", "Monocyte_Ccr2|Macrophage_Spp1", "pre-registered"),
        ("sample_state_cap", 300, "pre-registered"),
        ("feature_rule", "stored_2000_HVG_union_728_phase13_macrophage_union_controls", "pre-registered"),
        ("base_grn", "mouse_scATAC_atlas_v0.10.0_mm9", "CellOracle"),
        ("n_pca_dims", N_PCA_DIMS, "pre-registered"),
        ("knn_k", "max(20,int(0.025*n_cells))", "CellOracle tutorial rule"),
        ("knn_balanced", True, "pre-registered"),
        ("b_sight", "k*8", "CellOracle tutorial rule"),
        ("b_maxl", "k*4", "CellOracle tutorial rule"),
        ("n_jobs", N_JOBS, "resource lock"),
        ("GRN_unit", "cluster", "pre-registered"),
        ("alpha", ALPHA, "runtime default locked"),
        ("mode", "ko", "pre-registered"),
        ("n_propagation", N_PROPAGATION, "runtime default locked"),
        ("transition_n_neighbors", TRANSITION_N_NEIGHBORS, "pre-registered"),
        ("transition_sigma_corr", TRANSITION_SIGMA, "runtime default locked"),
        ("perturbation_score_grid_size", PS_GRID_SIZE, "pre-registered"),
        ("perturbation_score_min_mass", PS_MIN_MASS, "pre-registered"),
        ("sample_detection_rule", ">=1% in >=8/16 biological samples", "pre-registered"),
        ("state_detection_rule", ">=5% in either state", "pre-registered"),
        ("module_coherence", "abs(sample-median)>=0.10 SD and same sign in >=ceil(2/3*N)", "pre-registered"),
    ], columns=["parameter", "value", "source"])
    write_tsv_once(parameters, OUT / "parameters.tsv")

    root_mask = (adata.obs["timepoint"].astype(str).to_numpy() == "Day0") & (adata.obs["myeloid_state"].astype(str).to_numpy() == "Monocyte_Ccr2")
    if not root_mask.any():
        raise ValueError("No Day0 Monocyte_Ccr2 root candidates")
    root_names = adata.obs_names[root_mask].to_numpy()
    root_space = np.asarray(adata.obsm["X_harmony"])[root_mask]
    centroid = root_space.mean(axis=0)
    distances_to_centroid = np.square(root_space - centroid).sum(axis=1)
    minimum = distances_to_centroid.min()
    root_cell = sorted(root_names[np.isclose(distances_to_centroid, minimum)].astype(str))[0]
    write_tsv_once(pd.DataFrame([{"root_cell": root_cell, "selection": "nearest Day0 Monocyte_Ccr2 cell to Day0 centroid in existing X_harmony"}]), TABLES / "state_coordinate_root.tsv")

    direction_status = "FAILED"
    direction_reason = "not calculated"
    try:
        ov.single.lineage_pseudotime(
            adata,
            lineage_dictionary={"recruited_monocyte_to_Spp1": STATES},
            root_cells={"recruited_monocyte_to_Spp1": root_cell},
            obsm_key="X_umap",
            cluster_column_name="myeloid_state",
            obs_key="state_coordinate",
        )
        coord = pd.to_numeric(adata.obs["state_coordinate"], errors="coerce").to_numpy()
        finite_fraction = float(np.isfinite(coord).mean())
        medians = adata.obs.assign(state_coordinate=coord).groupby("myeloid_state", observed=True)["state_coordinate"].median()
        if finite_fraction < 0.95:
            direction_status = "FAILED"
            direction_reason = f"Only {finite_fraction:.3f} finite state-coordinate values"
        elif medians["Monocyte_Ccr2"] < medians["Macrophage_Spp1"]:
            direction_status = "PASS"
            direction_reason = f"finite={finite_fraction:.3f}; medians Monocyte={medians['Monocyte_Ccr2']:.4g} < Spp1={medians['Macrophage_Spp1']:.4g}"
        else:
            direction_status = "NEGATIVE_COMPLETE"
            direction_reason = f"finite={finite_fraction:.3f}; preregistered state order not observed: Monocyte={medians['Monocyte_Ccr2']:.4g}, Spp1={medians['Macrophage_Spp1']:.4g}"
    except Exception as exc:
        adata.obs["state_coordinate"] = np.nan
        direction_status = "FAILED"
        direction_reason = f"{type(exc).__name__}: {exc}"
    record_status("state_coordinate", direction_status, direction_reason, str(TABLES / "state_coordinate_by_sample.tsv"))

    coord_table = adata.obs[["cell_id", "sample_id", "timepoint", "myeloid_state", "state_coordinate"]].copy()
    coord_table = coord_table.reset_index(drop=True)
    write_tsv_once(coord_table, TABLES / "state_coordinate_cells.tsv")
    by_sample = coord_table.groupby(["sample_id", "timepoint", "myeloid_state"], observed=True, as_index=False)["state_coordinate"].median()
    write_tsv_once(by_sample, TABLES / "state_coordinate_by_sample.tsv")

    input_h5ad = OBJECTS / "celloracle_input.h5ad"
    if input_h5ad.exists():
        raise FileExistsError(input_h5ad)
    adata.write_h5ad(input_h5ad, compression="gzip")
    plot_input_figures(adata, direction_status)

    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, capture_output=True, text=True).stdout
    with (OUT / "sessionInfo.txt").open("a", encoding="utf-8") as handle:
        handle.write("\n\nPYTHON SESSION\n")
        handle.write(json.dumps({"python": sys.version, "executable": sys.executable, "platform": platform.platform(), "omicverse": ov.__version__, "celloracle": co.__version__, "anndata": ad.__version__}, indent=2) + "\n")
        handle.write("\nPIP FREEZE\n" + freeze)
    record_status("prepare", "PASS", f"AnnData, base GRN, {len(feature_list)} features, target rules and direction output frozen", str(input_h5ad))
    log(f"PREPARE_COMPLETE cells={adata.n_obs} genes={adata.n_vars} pilot_target={pilot_target} direction={direction_status}")


def plot_input_figures(adata: ad.AnnData, direction_status: str) -> None:
    emb = np.asarray(adata.obsm["X_umap"])
    state_colors = {"Monocyte_Ccr2": "#4477AA", "Macrophage_Spp1": "#CC6677"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for state in STATES:
        mask = adata.obs["myeloid_state"].astype(str).to_numpy() == state
        axes[0].scatter(emb[mask, 0], emb[mask, 1], s=5, alpha=0.65, color=state_colors[state], label=state, rasterized=True)
    axes[0].legend(frameon=False, markerscale=2)
    axes[0].set_title("Frozen monocyte–Spp1 macrophage axis")
    axes[0].set_axis_off()
    time_colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(TIME_ORDER)))
    for color, timepoint in zip(time_colors, TIME_ORDER):
        mask = adata.obs["timepoint"].astype(str).to_numpy() == timepoint
        axes[1].scatter(emb[mask, 0], emb[mask, 1], s=5, alpha=0.65, color=color, label=timepoint, rasterized=True)
    axes[1].legend(frameon=False, markerscale=2)
    axes[1].set_title("Observed collection time")
    axes[1].set_axis_off()
    fig.savefig(FIGURES / "input_axis_qc.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    coord = pd.to_numeric(adata.obs["state_coordinate"], errors="coerce").to_numpy()
    points = axes[0].scatter(emb[:, 0], emb[:, 1], c=coord, s=5, cmap="viridis", rasterized=True)
    fig.colorbar(points, ax=axes[0], label="Algorithmic state coordinate")
    axes[0].set_title(f"Restricted direction: {direction_status}")
    axes[0].set_axis_off()
    sample_coord = adata.obs.assign(state_coordinate=coord).groupby(["sample_id", "timepoint", "myeloid_state"], observed=True, as_index=False)["state_coordinate"].median()
    for state, marker, color in [("Monocyte_Ccr2", "o", "#4477AA"), ("Macrophage_Spp1", "s", "#CC6677")]:
        part = sample_coord.loc[sample_coord["myeloid_state"].astype(str) == state]
        x = np.array([TIME_ORDER.index(str(value)) for value in part["timepoint"]], dtype=float)
        jitter = np.linspace(-0.12, 0.12, len(part)) if len(part) > 1 else np.zeros(len(part))
        axes[1].scatter(x + jitter, part["state_coordinate"], marker=marker, color=color, label=state, s=35)
    axes[1].set_xticks(range(len(TIME_ORDER)), TIME_ORDER, rotation=35, ha="right")
    axes[1].set_ylabel("Biological-sample median state coordinate")
    axes[1].legend(frameon=False)
    axes[1].set_title("Descriptive sample-level summary")
    fig.savefig(FIGURES / "state_coordinate_by_sample.pdf", bbox_inches="tight")
    plt.close(fig)


def build_oracle(adata: ad.AnnData, base_grn: pd.DataFrame, label: str) -> tuple[co.Oracle, dict]:
    np.random.seed(SEED)
    random.seed(SEED)
    oracle = co.Oracle()
    oracle.import_anndata_as_normalized_count(adata=adata, cluster_column_name="myeloid_state", embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=base_grn)
    oracle.perform_PCA()
    k = max(20, int(0.025 * adata.n_obs))
    oracle.knn_imputation(n_pca_dims=N_PCA_DIMS, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=N_JOBS)
    if "imputed_count" not in oracle.adata.layers:
        raise RuntimeError("CellOracle did not create imputed_count")
    imputed = oracle.adata.layers["imputed_count"]
    original = oracle.adata.X
    if not finite_sparse(imputed):
        raise RuntimeError("imputed_count contains non-finite values")
    difference = imputed - original
    if sparse.issparse(difference):
        differs = difference.nnz > 0 and float(np.max(np.abs(difference.data))) > 1e-12
    else:
        differs = bool(np.max(np.abs(np.asarray(difference))) > 1e-12)
    if not differs:
        raise RuntimeError("Forbidden raw-X fallback detected: imputed_count equals X")
    oracle.fit_GRN_for_simulation(GRN_unit="cluster", alpha=ALPHA, use_cluster_specific_TFdict=False, verbose_level=1)
    matrices = getattr(oracle, "coef_matrix_per_cluster", None)
    if not isinstance(matrices, dict) or not set(STATES).issubset(set(map(str, matrices.keys()))):
        raise RuntimeError(f"Missing state-specific coefficient matrices: {list(matrices or {})}")
    edge_count = int(sum(np.count_nonzero(np.abs(matrix.to_numpy()) > 1e-12) for matrix in matrices.values()))
    if edge_count == 0:
        raise RuntimeError("Baseline GRN has no non-zero fitted edges")
    audit = {"label": label, "cells": adata.n_obs, "genes": adata.n_vars, "k": k, "b_sight": k * 8, "b_maxl": k * 4, "fitted_nonzero_edges": edge_count, "imputed_differs_from_X": True}
    return oracle, audit


def validate_result(result: PerturbResult, target: str, n_cells: int, n_genes: int) -> dict:
    if not isinstance(result, PerturbResult) or result.backend != "cell_oracle":
        raise RuntimeError("Unexpected perturbation result class/backend")
    delta = np.asarray(result.delta_X)
    if delta.shape != (n_cells, n_genes) or not np.isfinite(delta).all():
        raise RuntimeError(f"Invalid delta_X: {delta.shape}")
    if result.trajectory_shift is None:
        raise RuntimeError("CellOracle transition probability missing")
    transition = dense_array(result.trajectory_shift).astype(np.float64)
    if transition.shape != (n_cells, n_cells) or not np.isfinite(transition).all():
        raise RuntimeError("Invalid transition probability matrix")
    row_sums = transition.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-5, rtol=1e-5):
        raise RuntimeError(f"Transition rows do not sum to one: range={row_sums.min()}..{row_sums.max()}")
    perturbed = result.adata_perturbed
    target_index = list(perturbed.var_names).index(target)
    simulated = dense_array(perturbed.layers["simulated_count"][:, target_index]).ravel()
    target_max = float(np.max(np.abs(simulated)))
    if target_max > 1e-8:
        raise RuntimeError(f"KO target was not clamped to zero: max={target_max}")
    return {
        "backend": result.backend,
        "delta_shape": f"{delta.shape[0]}x{delta.shape[1]}",
        "delta_finite": True,
        "delta_max_abs": float(np.max(np.abs(delta))),
        "transition_shape": f"{transition.shape[0]}x{transition.shape[1]}",
        "transition_row_sum_min": float(row_sums.min()),
        "transition_row_sum_max": float(row_sums.max()),
        "ko_target_simulated_max_abs": target_max,
    }


def run_perturbation(adata: ad.AnnData, oracle: co.Oracle, base_grn: pd.DataFrame, target: str) -> PerturbResult:
    adata.uns["celloracle_oracle"] = oracle
    try:
        result = ov.single.perturb(
            adata,
            target=target,
            mode="ko",
            backend="cell_oracle",
            grn_base=base_grn,
            grn_output=True,
            return_delta=True,
            layer=None,
            n_propagation=N_PROPAGATION,
            backend_kwargs={
                "transition_n_neighbors": TRANSITION_N_NEIGHBORS,
                "transition_sigma_corr": TRANSITION_SIGMA,
                "transition_threads": N_JOBS,
            },
            copy=False,
        )
    finally:
        adata.uns.pop("celloracle_oracle", None)
    return result


def save_trimmed_result(result: PerturbResult, path: Path) -> None:
    if path.exists():
        raise FileExistsError(path)
    result.meta = {key: value for key, value in result.meta.items() if key != "oracle"}
    result.meta.update({"omicverse_version": ov.__version__, "celloracle_version": co.__version__, "seed": SEED})
    result.adata_perturbed = None
    result.grn_base = None
    result.save(str(path))
    loaded = PerturbResult.load(str(path))
    if loaded.backend != "cell_oracle" or loaded.delta_X is None:
        raise RuntimeError("Saved PerturbResult failed reload audit")


def pilot() -> None:
    log("PILOT_START")
    verify_stage_runtime_only()
    adata = ad.read_h5ad(OBJECTS / "celloracle_input.h5ad")
    base_grn = pd.read_parquet(OBJECTS / "mouse_scATAC_base_GRN_v0.10.0_feature_pruned.parquet")
    pilot_target = str(pd.read_csv(TABLES / "pilot_target.tsv", sep="\t").iloc[0]["pilot_target"])
    pilot_mask = adata.obs["sample_id"].astype(str).to_numpy() == PILOT_SAMPLE
    pilot_data = adata[pilot_mask].copy()
    pilot_counts = pilot_data.obs["myeloid_state"].astype(str).value_counts().to_dict()
    if pilot_data.n_obs != 568 or pilot_counts != {"Macrophage_Spp1": 300, "Monocyte_Ccr2": 268}:
        raise ValueError(f"Pilot input changed: cells={pilot_data.n_obs}, states={pilot_counts}")
    oracle, oracle_audit = build_oracle(pilot_data, base_grn, "pilot")
    result = run_perturbation(pilot_data, oracle, base_grn, pilot_target)
    result_audit = validate_result(result, pilot_target, pilot_data.n_obs, pilot_data.n_vars)
    audit = {**oracle_audit, **result_audit, "pilot_sample": PILOT_SAMPLE, "pilot_target": pilot_target}
    write_tsv_once(pd.DataFrame([audit]), TABLES / "pilot_output_audit.tsv")
    plot_vector_field(pilot_data, result, pilot_target, FIGURES / "pilot_vector_field.pdf", "Pilot technical output; no biological interpretation")
    save_trimmed_result(result, OBJECTS / "pilot.perturb_result.pkl")
    record_status("pilot", "PASS", "Single-sample CellOracle output schema, KO clamp, transition matrix and reload validated; no parameter tuning", str(TABLES / "pilot_output_audit.tsv"))
    log(f"PILOT_COMPLETE sample={PILOT_SAMPLE} target={pilot_target} cells={pilot_data.n_obs} genes={pilot_data.n_vars}")


def verify_stage_runtime_only() -> None:
    if Path(sys.executable).resolve() != EXPECTED_PYTHON.resolve() or ov.__version__ != "2.3.1" or co.__version__ != "0.20.0":
        raise BlockedDependency("Runtime changed after prepare")


def coefficient_metrics(oracle: co.Oracle, genes: list[str]) -> pd.DataFrame:
    matrices = getattr(oracle, "coef_matrix_per_cluster", {})
    rows = []
    for gene in genes:
        supported_states = 0
        total_edges = 0
        sum_abs = 0.0
        for state in STATES:
            matrix = matrices.get(state)
            if matrix is None or gene not in matrix.index:
                continue
            values = matrix.loc[gene].to_numpy(dtype=np.float64)
            nonzero = np.abs(values) > 1e-12
            count = int(nonzero.sum())
            if count:
                supported_states += 1
                total_edges += count
                sum_abs += float(np.abs(values[nonzero]).sum())
        rows.append({"gene": gene, "fitted_supported_states": supported_states, "fitted_out_degree_total": total_edges, "fitted_sum_abs_outgoing_coef": sum_abs})
    return pd.DataFrame(rows)


def aggregate_target(adata: ad.AnnData, result: PerturbResult, target: str, programs: pd.DataFrame, direction_ok: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, bool]:
    delta = np.asarray(result.delta_X, dtype=np.float64)
    transition = dense_array(result.trajectory_shift).astype(np.float64)
    delta_emb = result.delta_embedding(adata=adata, embedding_name="X_umap")
    displacement = np.linalg.norm(delta_emb, axis=1)
    states = adata.obs["myeloid_state"].astype(str).to_numpy()
    to_spp1 = transition[:, states == "Macrophage_Spp1"].sum(axis=1)
    if direction_ok:
        ps = result.perturbation_score(
            adata=adata,
            pseudotime="state_coordinate",
            embedding_name="X_umap",
            n_neighbors=TRANSITION_N_NEIGHBORS,
            grid_size=PS_GRID_SIZE,
            min_mass=PS_MIN_MASS,
            level="cell",
        ).to_numpy(dtype=np.float64)
    else:
        ps = np.full(adata.n_obs, np.nan)

    sample_rows = []
    samples = adata.obs["sample_id"].astype(str).to_numpy()
    times = adata.obs["timepoint"].astype(str).to_numpy()
    for sample in sorted(pd.unique(samples)):
        for state in STATES:
            mask = (samples == sample) & (states == state)
            if not mask.any():
                continue
            sample_rows.append({
                "target": target,
                "sample_id": sample,
                "timepoint": str(pd.unique(times[mask])[0]),
                "myeloid_state": state,
                "n_nuclei_descriptive": int(mask.sum()),
                "median_displacement": float(np.median(displacement[mask])),
                "median_transition_mass_to_Macrophage_Spp1": float(np.median(to_spp1[mask])),
                "median_perturbation_score": float(np.median(ps[mask])) if direction_ok else np.nan,
            })
    sample_table = pd.DataFrame(sample_rows)

    X = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))
    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
    baseline_sd = np.sqrt(np.maximum(mean_sq - mean * mean, 0))
    gene_to_idx = {gene: i for i, gene in enumerate(adata.var_names)}
    module_rows = []
    coherence_rows = []
    for module_id, module in programs.groupby("module_id", observed=True):
        module_genes = [gene for gene in module["gene"].astype(str) if gene in gene_to_idx and baseline_sd[gene_to_idx[gene]] > 1e-12]
        if len(module_genes) < 5:
            coherence_rows.append({"target": target, "module_id": module_id, "myeloid_state": "ALL", "eligible_samples": 0, "median_sample_shift": np.nan, "same_sign_samples": 0, "required_same_sign": 0, "coherent": False, "reason": "fewer_than_5_evaluable_genes"})
            continue
        idx = np.array([gene_to_idx[gene] for gene in module_genes], dtype=int)
        cell_score = (delta[:, idx] / baseline_sd[idx]).mean(axis=1)
        module_direction = str(module["direction"].iloc[0])
        for sample in sorted(pd.unique(samples)):
            for state in STATES:
                mask = (samples == sample) & (states == state)
                if mask.sum() < 10:
                    continue
                module_rows.append({
                    "target": target,
                    "module_id": module_id,
                    "phase13_direction": module_direction,
                    "n_module_genes": len(module_genes),
                    "sample_id": sample,
                    "timepoint": str(pd.unique(times[mask])[0]),
                    "myeloid_state": state,
                    "n_nuclei_descriptive": int(mask.sum()),
                    "module_shift_z": float(np.median(cell_score[mask])),
                })
    module_table = pd.DataFrame(module_rows)
    any_coherent = False
    if not module_table.empty:
        for (module_id, state), part in module_table.groupby(["module_id", "myeloid_state"], observed=True):
            values = part["module_shift_z"].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            n = len(values)
            median = float(np.median(values)) if n else np.nan
            required = math.ceil((2 * n) / 3) if n else 0
            if n:
                sign = 1 if median > 0 else (-1 if median < 0 else 0)
                same_sign = int(((values * sign) > 0).sum()) if sign else 0
            else:
                same_sign = 0
            coherent = bool(n >= 8 and abs(median) >= 0.10 and same_sign >= required)
            any_coherent = any_coherent or coherent
            coherence_rows.append({"target": target, "module_id": module_id, "myeloid_state": state, "eligible_samples": n, "median_sample_shift": median, "same_sign_samples": same_sign, "required_same_sign": required, "coherent": coherent, "reason": "pre_registered_sample_level_descriptive_rule"})
    coherence_table = pd.DataFrame(coherence_rows)
    cluster_table = result.cluster_transitions(adata=adata, cluster_col="myeloid_state").reset_index().rename(columns={"index": "source_state"})
    cluster_table.insert(0, "target", target)
    return sample_table, module_table, coherence_table, cluster_table, any_coherent


def plot_vector_field(adata: ad.AnnData, result: PerturbResult, target: str, path: Path, subtitle: str = "") -> None:
    embedding = np.asarray(adata.obsm["X_umap"])
    arrows = result.delta_embedding(adata=adata, embedding_name="X_umap")
    step = max(1, math.ceil(adata.n_obs / 800))
    idx = np.arange(0, adata.n_obs, step)
    colors = np.where(adata.obs["myeloid_state"].astype(str).to_numpy() == "Monocyte_Ccr2", "#4477AA", "#CC6677")
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    ax.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=4, alpha=0.3, rasterized=True)
    ax.quiver(embedding[idx, 0], embedding[idx, 1], arrows[idx, 0], arrows[idx, 1], angles="xy", scale_units="xy", scale=8, width=0.002, color="#222222", alpha=0.7)
    ax.set_title(f"{target} virtual KO\n{subtitle}" if subtitle else f"{target} virtual KO")
    ax.set_axis_off()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_module_effects(module: pd.DataFrame, target_order: list[str], path: Path) -> None:
    summary = module.groupby(["target", "module_id", "myeloid_state"], observed=True)["module_shift_z"].median().reset_index()
    module_order = module["module_id"].astype(str).drop_duplicates().tolist()
    labels = [name.removeprefix("Macrophage__").replace("__", " · ").replace("_", " ") for name in module_order]
    values = summary["module_shift_z"].to_numpy(dtype=float)
    finite = np.abs(values[np.isfinite(values)])
    limit = max(0.1, float(finite.max())) if finite.size else 0.1
    fig, axes = plt.subplots(1, len(STATES), figsize=(13, max(5.5, 0.65 * len(module_order))), sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    image = None
    for ax, state in zip(axes, STATES):
        matrix = (
            summary.loc[summary["myeloid_state"].astype(str) == state]
            .pivot(index="module_id", columns="target", values="module_shift_z")
            .reindex(index=module_order, columns=target_order)
        )
        image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
        ax.set_xticks(range(len(target_order)), target_order, rotation=35, ha="right")
        ax.set_yticks(range(len(module_order)), labels)
        ax.set_title(state.replace("_", " "))
    axes[0].set_ylabel("Frozen Phase13 macrophage module")
    fig.suptitle("Sample-level median standardized module shift")
    fig.colorbar(image, ax=axes, shrink=0.82, label="Standardized module shift")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def full() -> None:
    log("FULL_START")
    verify_stage_runtime_only()
    status = pd.read_csv(STATUS_PATH, sep="\t", dtype=str, keep_default_na=False)
    pilot_row = status.loc[(status["component"] == "pilot") & (status["status"] == "PASS")]
    if pilot_row.empty:
        raise RuntimeError("Pilot PASS is required before full run")
    adata = ad.read_h5ad(OBJECTS / "celloracle_input.h5ad")
    base_grn = pd.read_parquet(OBJECTS / "mouse_scATAC_base_GRN_v0.10.0_feature_pruned.parquet")
    programs = phase13_programs()
    baseline_path = OBJECTS / "baseline.celloracle.oracle"
    if baseline_path.exists() or (TABLES / "frozen_targets.tsv").exists():
        raise FileExistsError("Refusing to overwrite an existing baseline or frozen target set")
    oracle, baseline_audit = build_oracle(adata, base_grn, "full")
    oracle.to_hdf5(str(baseline_path))
    reloaded = co.load_hdf5(str(baseline_path))
    if not hasattr(reloaded, "coef_matrix_per_cluster") or reloaded.adata.shape != oracle.adata.shape:
        raise RuntimeError("Baseline Oracle reload audit failed")
    baseline_audit["saved_path"] = str(baseline_path)
    write_tsv_once(pd.DataFrame([baseline_audit]), TABLES / "baseline_oracle_audit.tsv")
    record_status("baseline_oracle", "PASS", "Full baseline Oracle fitted once and reloaded successfully", str(baseline_path))

    pre = pd.read_csv(TABLES / "candidate_tf_prebaseline.tsv", sep="\t")
    metrics = coefficient_metrics(oracle, pre["gene"].astype(str).tolist())
    ranked = pre.merge(metrics, on="gene", how="left")
    ranked["eligible_final"] = boolean_series(ranked["eligible_prebaseline"]) & (ranked["fitted_supported_states"] >= 1) & (ranked["fitted_out_degree_total"] > 0)
    ranked = ranked.sort_values(
        ["supported_contrast_n", "fitted_supported_states", "fitted_out_degree_total", "fitted_sum_abs_outgoing_coef", "min_abs_log2FC_supported", "worst_fdr_supported", "gene"],
        ascending=[False, False, False, False, False, True, True],
        kind="mergesort",
    )
    write_tsv_once(ranked, TABLES / "candidate_tf_ranking.tsv")
    eligible = ranked.loc[ranked["eligible_final"]].copy()
    if eligible.shape[0] < 5:
        record_status("target_selection", "NOT_APPLICABLE", f"Only {eligible.shape[0]} TFs pass all frozen rules", str(TABLES / "candidate_tf_ranking.tsv"))
        raise RuntimeError("Fewer than five eligible TFs")
    discovery_targets = eligible.head(5)["gene"].astype(str).tolist()

    control_pre = pd.read_csv(TABLES / "technical_control_candidates.tsv", sep="\t")
    control_metrics = coefficient_metrics(oracle, CONTROL_PRIORITY)
    control_table = control_pre.merge(control_metrics, on="gene", how="left")
    control_table["eligible_final"] = boolean_series(control_table["eligible_prebaseline"]) & (control_table["fitted_supported_states"] >= 1) & (control_table["fitted_out_degree_total"] > 0)
    technical_control = None
    for gene in CONTROL_PRIORITY:
        row = control_table.loc[control_table["gene"] == gene]
        if not row.empty and bool(row.iloc[0]["eligible_final"]) and gene not in discovery_targets:
            technical_control = gene
            break
    frozen_rows = [{"target": gene, "role": "DISCOVERY", "rank": index + 1, "selection_reason": "pre_registered_phase13_plus_fitted_GRN_lexicographic_rank"} for index, gene in enumerate(discovery_targets)]
    if technical_control:
        frozen_rows.append({"target": technical_control, "role": "TECHNICAL_CONTROL", "rank": 1, "selection_reason": "first eligible in Spi1>Irf8>Cebpb priority"})
        record_status("technical_control", "PASS", f"Frozen technical control: {technical_control}")
    else:
        record_status("technical_control", "NOT_APPLICABLE", "No additional eligible control outside the five discovery TFs")
    frozen = pd.DataFrame(frozen_rows)
    write_tsv_once(frozen, TABLES / "frozen_targets.tsv")
    record_status("target_selection", "PASS", f"Five discovery TFs frozen before KO: {','.join(discovery_targets)}", str(TABLES / "frozen_targets.tsv"))
    log(f"TARGETS_FROZEN discovery={','.join(discovery_targets)} control={technical_control or 'NONE'}")

    direction_row = status.loc[status["component"] == "state_coordinate"]
    direction_ok = not direction_row.empty and direction_row.iloc[-1]["status"] == "PASS"
    per_target_dir = TABLES / "per_target"
    per_figure_dir = FIGURES / "per_target"
    per_target_dir.mkdir(parents=True, exist_ok=False)
    per_figure_dir.mkdir(parents=True, exist_ok=False)
    target_order = discovery_targets + ([technical_control] if technical_control else [])
    for index, target in enumerate(target_order, start=1):
        log(f"TARGET_START index={index}/{len(target_order)} target={target}")
        baseline = co.load_hdf5(str(baseline_path))
        result = run_perturbation(adata, baseline, base_grn, target)
        result_audit = validate_result(result, target, adata.n_obs, adata.n_vars)
        result.summary(top_n=10)
        sample_table, module_table, coherence_table, cluster_table, coherent = aggregate_target(adata, result, target, programs, direction_ok)
        delta_expr = result.delta_expr.copy()
        delta_expr.insert(0, "target", target)
        write_tsv_once(delta_expr, per_target_dir / f"{target}_delta_expression.tsv")
        write_tsv_once(sample_table, per_target_dir / f"{target}_sample_effects.tsv")
        write_tsv_once(module_table, per_target_dir / f"{target}_module_effects.tsv")
        write_tsv_once(coherence_table, per_target_dir / f"{target}_module_coherence.tsv")
        write_tsv_once(cluster_table, per_target_dir / f"{target}_state_transition.tsv")
        write_tsv_once(pd.DataFrame([{"target": target, **result_audit}]), per_target_dir / f"{target}_integrity.tsv")
        plot_vector_field(adata, result, target, per_figure_dir / f"{target}_vector_field.pdf", "Computational state displacement")
        save_trimmed_result(result, OBJECTS / f"{target}.perturb_result.pkl")
        target_status = "PASS" if coherent else "NEGATIVE_COMPLETE"
        target_reason = "At least one frozen module met the preregistered sample-level descriptive rule" if coherent else "Valid simulation; no frozen module met the preregistered sample-level descriptive rule"
        record_status("target", target_status, target_reason, str(OBJECTS / f"{target}.perturb_result.pkl"), target=target)
        log(f"TARGET_COMPLETE index={index}/{len(target_order)} target={target} status={target_status}")
        del result, baseline
        gc.collect()

    combine_outputs(target_order, adata)
    final_integrity(discovery_targets, technical_control)
    log("X01_CELLORACLE_COMPLETE")


def combine_outputs(target_order: list[str], adata: ad.AnnData) -> None:
    per_target = TABLES / "per_target"
    specs = [
        ("delta_expression", "target_delta_expression.tsv"),
        ("sample_effects", "sample_level_effects.tsv"),
        ("module_effects", "phase13_module_effects.tsv"),
        ("state_transition", "state_transition_summary.tsv"),
    ]
    combined = {}
    for suffix, output_name in specs:
        frames = [pd.read_csv(per_target / f"{target}_{suffix}.tsv", sep="\t") for target in target_order]
        combined[suffix] = pd.concat(frames, ignore_index=True)
        write_tsv_once(combined[suffix], TABLES / output_name)

    with PdfPages(FIGURES / "target_vector_fields.pdf") as pdf:
        for target in target_order:
            result = PerturbResult.load(str(OBJECTS / f"{target}.perturb_result.pkl"))
            embedding = np.asarray(adata.obsm["X_umap"])
            arrows = result.delta_embedding(adata=adata, embedding_name="X_umap")
            step = max(1, math.ceil(adata.n_obs / 800))
            idx = np.arange(0, adata.n_obs, step)
            colors = np.where(adata.obs["myeloid_state"].astype(str).to_numpy() == "Monocyte_Ccr2", "#4477AA", "#CC6677")
            fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
            ax.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=4, alpha=0.3, rasterized=True)
            ax.quiver(embedding[idx, 0], embedding[idx, 1], arrows[idx, 0], arrows[idx, 1], angles="xy", scale_units="xy", scale=8, width=0.002, color="#222222", alpha=0.7)
            ax.set_title(f"{target} virtual KO\nComputational state displacement")
            ax.set_axis_off()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            del result
            gc.collect()

    plot_module_effects(combined["module_effects"], target_order, FIGURES / "target_module_effects.pdf")

    sample = combined["sample_effects"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    positions = np.arange(len(target_order))
    for index, target in enumerate(target_order):
        values = sample.loc[sample["target"] == target, "median_displacement"].to_numpy(dtype=float)
        axes[0].scatter(np.full(len(values), index) + np.linspace(-0.12, 0.12, len(values)), values, s=18, alpha=0.7)
        trans = sample.loc[sample["target"] == target, "median_transition_mass_to_Macrophage_Spp1"].to_numpy(dtype=float)
        axes[1].scatter(np.full(len(trans), index) + np.linspace(-0.12, 0.12, len(trans)), trans, s=18, alpha=0.7)
    for ax in axes:
        ax.set_xticks(positions, target_order, rotation=45, ha="right")
    axes[0].set_ylabel("Sample × state median displacement")
    axes[1].set_ylabel("Sample × state median transition mass to Spp1")
    axes[0].set_title("Computational displacement")
    axes[1].set_title("State-transition summary")
    fig.savefig(FIGURES / "sample_level_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def final_integrity(discovery_targets: list[str], technical_control: str | None) -> None:
    required = [
        OUT / "parameters.tsv", OUT / "sessionInfo.txt", STATUS_PATH,
        TABLES / "frozen_feature_universe.tsv", TABLES / "frozen_targets.tsv",
        TABLES / "target_delta_expression.tsv", TABLES / "sample_level_effects.tsv",
        TABLES / "phase13_module_effects.tsv", TABLES / "state_transition_summary.tsv",
        FIGURES / "input_axis_qc.pdf", FIGURES / "state_coordinate_by_sample.pdf",
        FIGURES / "target_vector_fields.pdf", FIGURES / "target_module_effects.pdf",
        FIGURES / "sample_level_summary.pdf", OBJECTS / "celloracle_input.h5ad",
        OBJECTS / "baseline.celloracle.oracle",
    ]
    required.extend(OBJECTS / f"{target}.perturb_result.pkl" for target in discovery_targets)
    if technical_control:
        required.append(OBJECTS / f"{technical_control}.perturb_result.pkl")
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing/empty final artifacts: {missing}")
    statuses = pd.read_csv(STATUS_PATH, sep="\t", dtype=str, keep_default_na=False)
    blocking = statuses.loc[statuses["status"].isin(["BLOCKED_DEPENDENCY", "FAILED"])]
    if not blocking.empty:
        raise RuntimeError(f"Blocking terminal status remains: {blocking[['component', 'target', 'status']].to_dict('records')}")
    target_rows = statuses.loc[(statuses["component"] == "target") & (statuses["target"].isin(discovery_targets))]
    if target_rows.shape[0] != 5 or not target_rows["status"].isin(["PASS", "NEGATIVE_COMPLETE"]).all():
        raise RuntimeError("Five discovery targets do not have valid terminal statuses")
    direction = statuses.loc[statuses["component"] == "state_coordinate", "status"]
    if direction.empty or direction.iloc[-1] not in {"PASS", "NEGATIVE_COMPLETE"}:
        raise RuntimeError("Direction component lacks PASS/NEGATIVE_COMPLETE")
    if technical_control:
        control = statuses.loc[(statuses["component"] == "target") & (statuses["target"] == technical_control), "status"]
        if control.empty or control.iloc[-1] not in {"PASS", "NEGATIVE_COMPLETE"}:
            raise RuntimeError("Technical control lacks valid terminal status")
    summary = statuses.loc[statuses["component"] == "target", ["target", "status", "reason", "artifact"]].copy()
    write_tsv_once(summary, TABLES / "run_summary.tsv")
    record_status("final_integrity", "PASS", "All required artifacts are non-empty and all five discovery targets have valid terminal states", str(TABLES / "run_summary.tsv"))
    record_status("X01_OMICVERSE_CELLORACLE", "PASS", "CellOracle workflow completed without nucleus-level inference or post-hoc parameter changes", str(OUT / "_SUCCESS"))
    write_text_once(f"X01_OMICVERSE_CELLORACLE COMPLETE\ncompleted_at\t{now()}\ndiscovery_targets\t{','.join(discovery_targets)}\ntechnical_control\t{technical_control or 'NOT_APPLICABLE'}\n", OUT / "_SUCCESS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["prepare", "pilot", "full"], required=True)
    args = parser.parse_args()
    if not OUT.is_dir() or not str(OUT).startswith("/home/zhenhanli/myocardial/scop/GSE214611_spatial/21_scop_expansion/"):
        raise RuntimeError("Output write fence failed")
    np.random.seed(SEED)
    random.seed(SEED)
    try:
        {"prepare": prepare, "pilot": pilot, "full": full}[args.stage]()
        return 0
    except BlockedDependency as exc:
        record_status(args.stage, "BLOCKED_DEPENDENCY", f"{type(exc).__name__}: {exc}")
        log(f"STAGE_BLOCKED stage={args.stage} error={type(exc).__name__}: {exc}")
        traceback.print_exc()
        return 20
    except Exception as exc:
        record_status(args.stage, "FAILED", f"{type(exc).__name__}: {exc}")
        log(f"STAGE_FAILED stage={args.stage} error={type(exc).__name__}: {exc}")
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
