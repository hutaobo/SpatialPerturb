"""Core analysis tools for spatial perturbation inference."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import mannwhitneyu, pearsonr, spearmanr, ttest_ind
import statsmodels.api as sm

from ._utils import benjamini_hochberg, extract_matrix, safe_log2_fold_change
from .gr import collect_neighbors
from .schema import ensure_spatialperturb_schema
from .signatures import compare_program_concordance, derive_perturbation_programs

_FALLBACK_LR_NETWORK = pd.DataFrame(
    [
        {"ligand": "CXCL10", "receptor": "CXCR3", "resource": "fallback"},
        {"ligand": "GAS6", "receptor": "AXL", "resource": "fallback"},
        {"ligand": "NRG1", "receptor": "ERBB4", "resource": "fallback"},
        {"ligand": "TGFB1", "receptor": "TGFBR1", "resource": "fallback"},
        {"ligand": "LIG1", "receptor": "REC1", "resource": "fallback"},
    ]
)


def _default_feature_names(adata: AnnData, features: Sequence[str] | None) -> list[str]:
    if features is not None:
        return [str(feature) for feature in features if str(feature) in adata.var_names]
    barcode_columns = set(adata.uns.get("spatialperturb", {}).get("barcode_columns", []))
    return [str(feature) for feature in adata.var_names if str(feature) not in barcode_columns]


def _restrict_mask(
    adata: AnnData,
    *,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
) -> np.ndarray:
    mask = np.ones(adata.n_obs, dtype=bool)
    if cell_type is not None:
        cell_types = {cell_type} if isinstance(cell_type, str) else set(map(str, cell_type))
        mask &= adata.obs["cell_type"].astype(str).isin(cell_types).to_numpy()
    if roi is not None:
        rois = {roi} if isinstance(roi, str) else set(map(str, roi))
        mask &= adata.obs["roi"].astype(str).isin(rois).to_numpy()
    return mask


def _group_mask(
    adata: AnnData,
    *,
    perturbation: str,
    groupby: str,
    status_col: str,
) -> np.ndarray:
    include = (
        adata.obs["include_for_inference"].astype(bool).to_numpy()
        if "include_for_inference" in adata.obs
        else np.ones(adata.n_obs, dtype=bool)
    )
    return (
        (adata.obs[groupby].astype(str) == str(perturbation))
        & (adata.obs[status_col].astype(str) == "single")
        & include
    ).to_numpy()


def _mode_or_first(series: pd.Series) -> object:
    mode = series.mode(dropna=True)
    if not mode.empty:
        return mode.iloc[0]
    return series.iloc[0] if len(series) else np.nan


def _build_de_output(
    *,
    perturbation: str,
    control: str,
    effect_type: str,
    features: list[str],
    case_mean: np.ndarray,
    control_mean: np.ndarray,
    statistics: np.ndarray,
    pvalues: np.ndarray,
    case_n: int,
    control_n: int,
    method: str,
    case_sample_n: int | None = None,
    control_sample_n: int | None = None,
) -> pd.DataFrame:
    fdr = benjamini_hochberg(pvalues)
    output = pd.DataFrame(
        {
            "perturbation": perturbation,
            "control": control,
            "effect_type": effect_type,
            "method": method,
            "gene": features,
            "case_n": case_n,
            "control_n": control_n,
            "case_sample_n": case_sample_n,
            "control_sample_n": control_sample_n,
            "mean_case": case_mean,
            "mean_control": control_mean,
            "log2fc": safe_log2_fold_change(case_mean, control_mean),
            "statistic": statistics,
            "pvalue": pvalues,
            "fdr": fdr,
        }
    )
    return output.sort_values(["fdr", "pvalue", "log2fc"], ascending=[True, True, False]).reset_index(drop=True)


def _run_simple_de_matrix(
    *,
    case_matrix: np.ndarray,
    control_matrix: np.ndarray,
    features: list[str],
    perturbation: str,
    control: str,
    effect_type: str,
    effect_size_only: bool = False,
) -> pd.DataFrame:
    case_mean = case_matrix.mean(axis=0)
    control_mean = control_matrix.mean(axis=0)
    statistics = np.zeros(len(features), dtype=float)
    pvalues = np.ones(len(features), dtype=float)
    if not effect_size_only:
        for idx in range(len(features)):
            case_gene = case_matrix[:, idx]
            control_gene = control_matrix[:, idx]
            if case_gene.shape == control_gene.shape and np.allclose(case_gene, control_gene):
                continue
            result = mannwhitneyu(case_gene, control_gene, alternative="two-sided")
            statistics[idx] = float(result.statistic)
            pvalues[idx] = float(result.pvalue)

    return _build_de_output(
        perturbation=perturbation,
        control=control,
        effect_type=effect_type,
        features=features,
        case_mean=case_mean,
        control_mean=control_mean,
        statistics=statistics,
        pvalues=pvalues,
        case_n=case_matrix.shape[0],
        control_n=control_matrix.shape[0],
        method="simple",
    )


def _aggregate_pseudobulk(
    matrix: np.ndarray,
    obs: pd.DataFrame,
    *,
    sample_col: str,
    features: list[str],
    label: str,
    covariates: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = obs[sample_col].astype(str).rename("sample")
    frame = pd.DataFrame(matrix, index=obs.index, columns=features)
    aggregated = frame.groupby(sample_ids).sum()
    aggregated.index = [f"{label}:{sample}" for sample in aggregated.index]

    meta_cols = [sample_col, *(covariates or [])]
    meta = obs.loc[:, [col for col in meta_cols if col in obs.columns]].copy()
    meta["_sample"] = sample_ids.values
    sample_meta = meta.groupby("_sample").agg(_mode_or_first)
    sample_meta.index = aggregated.index
    sample_meta["sample"] = [idx.split(":", 1)[1] for idx in aggregated.index]
    sample_meta["label"] = label
    return aggregated, sample_meta


def _run_pseudobulk_de_matrix(
    *,
    case_matrix: np.ndarray,
    control_matrix: np.ndarray,
    case_obs: pd.DataFrame,
    control_obs: pd.DataFrame,
    features: list[str],
    perturbation: str,
    control: str,
    effect_type: str,
    sample_col: str,
    covariates: Sequence[str] | None = None,
    effect_size_only: bool = False,
) -> pd.DataFrame:
    case_counts, case_meta = _aggregate_pseudobulk(
        case_matrix,
        case_obs,
        sample_col=sample_col,
        features=features,
        label="case",
        covariates=covariates,
    )
    control_counts, control_meta = _aggregate_pseudobulk(
        control_matrix,
        control_obs,
        sample_col=sample_col,
        features=features,
        label="control",
        covariates=covariates,
    )

    combined_counts = pd.concat([case_counts, control_counts], axis=0)
    lib_size = combined_counts.sum(axis=1).replace(0, 1.0)
    normalized = np.log1p(combined_counts.div(lib_size, axis=0) * 1_000_000.0)

    design_meta = pd.concat([case_meta, control_meta], axis=0)
    design_meta["condition"] = [1] * len(case_meta) + [0] * len(control_meta)
    design = pd.DataFrame({"condition": design_meta["condition"].astype(float)}, index=normalized.index)
    if covariates:
        encoded = pd.get_dummies(design_meta.loc[:, list(covariates)], drop_first=True, dtype=float)
        design = pd.concat([design, encoded], axis=1)
    design = sm.add_constant(design, has_constant="add")

    statistics = np.zeros(len(features), dtype=float)
    pvalues = np.ones(len(features), dtype=float)
    if not effect_size_only:
        for idx, gene in enumerate(features):
            y = normalized[gene].astype(float)
            try:
                fit = sm.OLS(y, design).fit()
                statistics[idx] = float(fit.tvalues.get("condition", np.nan))
                pvalues[idx] = float(fit.pvalues.get("condition", np.nan))
            except Exception:
                fallback = ttest_ind(
                    normalized.loc[case_counts.index, gene],
                    normalized.loc[control_counts.index, gene],
                    equal_var=False,
                    nan_policy="omit",
                )
                statistics[idx] = float(fallback.statistic)
                pvalues[idx] = float(fallback.pvalue)

    return _build_de_output(
        perturbation=perturbation,
        control=control,
        effect_type=effect_type,
        features=features,
        case_mean=case_counts.mean(axis=0).to_numpy(),
        control_mean=control_counts.mean(axis=0).to_numpy(),
        statistics=statistics,
        pvalues=pvalues,
        case_n=case_matrix.shape[0],
        control_n=control_matrix.shape[0],
        case_sample_n=len(case_counts),
        control_sample_n=len(control_counts),
        method="pseudobulk",
    )


def _run_de(
    adata: AnnData,
    *,
    case_cells: Sequence[str] | None = None,
    control_cells: Sequence[str] | None = None,
    case_matrix: np.ndarray | None = None,
    control_matrix: np.ndarray | None = None,
    case_obs: pd.DataFrame | None = None,
    control_obs: pd.DataFrame | None = None,
    perturbation: str,
    control: str,
    effect_type: str,
    features: Sequence[str] | None = None,
    layer: str | None = None,
    method: str = "simple",
    sample_col: str | None = None,
    covariates: Sequence[str] | None = None,
    effect_size_only: bool = False,
) -> pd.DataFrame:
    resolved_features = _default_feature_names(adata, features)
    if not resolved_features:
        raise ValueError("No features available for differential expression.")

    if case_matrix is None or control_matrix is None:
        if case_cells is None or control_cells is None:
            raise ValueError("Either cells or explicit matrices must be provided.")
        case_matrix = extract_matrix(adata, layer=layer, obs_names=case_cells, var_names=resolved_features)
        control_matrix = extract_matrix(adata, layer=layer, obs_names=control_cells, var_names=resolved_features)
        case_obs = adata.obs.loc[list(case_cells)].copy()
        control_obs = adata.obs.loc[list(control_cells)].copy()
    else:
        if case_obs is None or control_obs is None:
            raise ValueError("case_obs and control_obs are required when explicit matrices are provided.")
        case_obs = case_obs.copy()
        control_obs = control_obs.copy()

    if method == "simple":
        return _run_simple_de_matrix(
            case_matrix=case_matrix,
            control_matrix=control_matrix,
            features=resolved_features,
            perturbation=perturbation,
            control=control,
            effect_type=effect_type,
            effect_size_only=effect_size_only,
        )
    if method == "pseudobulk":
        if sample_col is None:
            raise ValueError("sample_col is required when method='pseudobulk'.")
        return _run_pseudobulk_de_matrix(
            case_matrix=case_matrix,
            control_matrix=control_matrix,
            case_obs=case_obs,
            control_obs=control_obs,
            features=resolved_features,
            perturbation=perturbation,
            control=control,
            effect_type=effect_type,
            sample_col=sample_col,
            covariates=covariates,
            effect_size_only=effect_size_only,
        )
    raise ValueError("method must be 'simple' or 'pseudobulk'.")


def intrinsic_de(
    adata: AnnData,
    *,
    perturbation: str,
    control: str,
    groupby: str = "perturbation",
    status_col: str = "perturbation_status",
    features: Sequence[str] | None = None,
    layer: str | None = None,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
    min_cells_per_group: int = 2,
    min_samples_per_group: int = 2,
    method: str = "simple",
    sample_col: str | None = None,
    covariates: Sequence[str] | None = None,
    effect_size_only: bool = False,
    min_cells: int | None = None,
) -> pd.DataFrame:
    """Compare perturbed cells against control cells in the same ROI/cell type context."""
    ensure_spatialperturb_schema(adata)
    if min_cells is not None:
        min_cells_per_group = int(min_cells)
    base_mask = _restrict_mask(adata, cell_type=cell_type, roi=roi)
    case_mask = _group_mask(adata, perturbation=perturbation, groupby=groupby, status_col=status_col) & base_mask
    control_mask = _group_mask(adata, perturbation=control, groupby=groupby, status_col=status_col) & base_mask

    if case_mask.sum() < min_cells_per_group or control_mask.sum() < min_cells_per_group:
        raise ValueError(
            "Need at least "
            f"{min_cells_per_group} cells in both case and control groups, found "
            f"{case_mask.sum()} and {control_mask.sum()}."
        )
    if method == "pseudobulk":
        if sample_col is None:
            raise ValueError("sample_col is required when method='pseudobulk'.")
        case_sample_n = adata.obs.loc[case_mask, sample_col].astype(str).nunique()
        control_sample_n = adata.obs.loc[control_mask, sample_col].astype(str).nunique()
        if case_sample_n < min_samples_per_group or control_sample_n < min_samples_per_group:
            raise ValueError(
                "Need at least "
                f"{min_samples_per_group} samples in both case and control groups for pseudobulk, found "
                f"{case_sample_n} and {control_sample_n}."
            )

    return _run_de(
        adata,
        case_cells=adata.obs_names[case_mask],
        control_cells=adata.obs_names[control_mask],
        perturbation=perturbation,
        control=control,
        effect_type="intrinsic",
        features=features,
        layer=layer,
        method=method,
        sample_col=sample_col,
        covariates=covariates,
        effect_size_only=effect_size_only,
    )


def neighbor_de(
    adata: AnnData,
    *,
    perturbation: str,
    control: str,
    graph_key: str | None = None,
    groupby: str = "perturbation",
    status_col: str = "perturbation_status",
    exclude_perturbed: bool = True,
    drop_shared_neighbors: bool = False,
    weight_by_distance: bool = False,
    aggregate: str = "mean",
    features: Sequence[str] | None = None,
    layer: str | None = None,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
    min_cells_per_group: int = 2,
    min_samples_per_group: int = 2,
    method: str = "simple",
    sample_col: str | None = None,
    covariates: Sequence[str] | None = None,
    effect_size_only: bool = False,
    min_cells: int | None = None,
) -> pd.DataFrame:
    """Compare neighbors around perturbed cells against neighbors around control cells."""
    ensure_spatialperturb_schema(adata)
    if min_cells is not None:
        min_cells_per_group = int(min_cells)
    if aggregate not in {"mean", "sum", "pseudobulk"}:
        raise ValueError("aggregate must be one of {'mean', 'sum', 'pseudobulk'}.")
    base_mask = _restrict_mask(adata, cell_type=cell_type, roi=roi)
    neighbor_mask = _restrict_mask(adata, roi=roi) if roi is not None else np.ones(adata.n_obs, dtype=bool)
    allowed_neighbors = set(adata.obs_names[neighbor_mask])
    case_sources = adata.obs_names[
        _group_mask(adata, perturbation=perturbation, groupby=groupby, status_col=status_col) & base_mask
    ]
    control_sources = adata.obs_names[
        _group_mask(adata, perturbation=control, groupby=groupby, status_col=status_col) & base_mask
    ]

    case_neighbors = collect_neighbors(
        adata,
        cells=case_sources,
        graph_key=graph_key,
        exclude_perturbed=exclude_perturbed,
        as_frame=True,
    )
    control_neighbors = collect_neighbors(
        adata,
        cells=control_sources,
        graph_key=graph_key,
        exclude_perturbed=exclude_perturbed,
        as_frame=True,
    )

    if not case_neighbors.empty:
        case_neighbors = case_neighbors[case_neighbors["neighbor"].isin(allowed_neighbors)].copy()
    if not control_neighbors.empty:
        control_neighbors = control_neighbors[control_neighbors["neighbor"].isin(allowed_neighbors)].copy()

    case_cells = set(case_neighbors["neighbor"]) if not case_neighbors.empty else set()
    control_cells = set(control_neighbors["neighbor"]) if not control_neighbors.empty else set()
    shared = case_cells & control_cells
    if drop_shared_neighbors:
        case_neighbors = case_neighbors[~case_neighbors["neighbor"].isin(shared)].copy()
        control_neighbors = control_neighbors[~control_neighbors["neighbor"].isin(shared)].copy()

    resolved_features = _default_feature_names(adata, features)
    if not resolved_features:
        raise ValueError("No features available for differential expression.")

    def _aggregate_frame(frame: pd.DataFrame, sources: Sequence[str]) -> tuple[np.ndarray, pd.DataFrame]:
        if frame.empty:
            return np.empty((0, len(resolved_features)), dtype=float), adata.obs.iloc[0:0].copy()

        if aggregate == "pseudobulk":
            if sample_col is None:
                raise ValueError("sample_col is required when aggregate='pseudobulk'.")
            records: list[np.ndarray] = []
            obs_records: list[pd.Series] = []
            source_obs = adata.obs.loc[list(sources)].copy()
            for sample_id, sample_sources in source_obs.groupby(sample_col):
                source_names = set(sample_sources.index.astype(str))
                sample_edges = frame[frame["source"].isin(source_names)]
                if sample_edges.empty:
                    continue
                neighbor_names = list(dict.fromkeys(sample_edges["neighbor"].astype(str)))
                expr = extract_matrix(adata, layer=layer, obs_names=neighbor_names, var_names=resolved_features)
                if weight_by_distance:
                    min_distance = sample_edges.groupby("neighbor")["distance"].min().reindex(neighbor_names)
                    weights = 1.0 / np.maximum(min_distance.to_numpy(dtype=float), 1e-6)
                    row = (expr * weights[:, None]).sum(axis=0)
                else:
                    row = expr.sum(axis=0)
                records.append(np.asarray(row, dtype=float))
                meta = sample_sources.iloc[0].copy()
                meta[sample_col] = str(sample_id)
                meta["n_neighbors"] = int(len(neighbor_names))
                meta["n_sources"] = int(len(source_names))
                obs_records.append(meta)
            if not records:
                return np.empty((0, len(resolved_features)), dtype=float), adata.obs.iloc[0:0].copy()
            return np.vstack(records), pd.DataFrame(obs_records)

        records = []
        obs_records = []
        for source, edges in frame.groupby("source"):
            neighbor_names = edges["neighbor"].astype(str).tolist()
            if not neighbor_names:
                continue
            expr = extract_matrix(adata, layer=layer, obs_names=neighbor_names, var_names=resolved_features)
            if weight_by_distance:
                weights = 1.0 / np.maximum(edges["distance"].to_numpy(dtype=float), 1e-6)
                if aggregate == "mean":
                    row = np.average(expr, axis=0, weights=weights)
                else:
                    row = (expr * weights[:, None]).sum(axis=0)
            else:
                row = expr.mean(axis=0) if aggregate == "mean" else expr.sum(axis=0)
            meta = adata.obs.loc[str(source)].copy()
            meta["n_neighbors"] = int(len(neighbor_names))
            obs_records.append(meta)
            records.append(np.asarray(row, dtype=float))

        if not records:
            return np.empty((0, len(resolved_features)), dtype=float), adata.obs.iloc[0:0].copy()
        return np.vstack(records), pd.DataFrame(obs_records)

    case_matrix, case_obs = _aggregate_frame(case_neighbors, case_sources)
    control_matrix, control_obs = _aggregate_frame(control_neighbors, control_sources)
    if case_matrix.shape[0] < min_cells_per_group or control_matrix.shape[0] < min_cells_per_group:
        raise ValueError(
            "Need at least "
            f"{min_cells_per_group} aggregated neighbor profiles in both case and control groups, found "
            f"{case_matrix.shape[0]} and {control_matrix.shape[0]}."
        )
    if method == "pseudobulk":
        if sample_col is None:
            raise ValueError("sample_col is required when method='pseudobulk'.")
        case_sample_n = case_obs[sample_col].astype(str).nunique()
        control_sample_n = control_obs[sample_col].astype(str).nunique()
        if case_sample_n < min_samples_per_group or control_sample_n < min_samples_per_group:
            raise ValueError(
                "Need at least "
                f"{min_samples_per_group} samples in both case and control neighbor groups for pseudobulk, found "
                f"{case_sample_n} and {control_sample_n}."
            )

    output = _run_de(
        adata,
        case_matrix=case_matrix,
        control_matrix=control_matrix,
        case_obs=case_obs,
        control_obs=control_obs,
        perturbation=perturbation,
        control=control,
        effect_type="neighbor",
        features=resolved_features,
        layer=layer,
        method=method,
        sample_col=sample_col,
        covariates=covariates,
        effect_size_only=effect_size_only,
    )
    output["aggregate"] = aggregate
    output["weight_by_distance"] = weight_by_distance
    output["case_source_n"] = len(case_sources)
    output["control_source_n"] = len(control_sources)
    output["case_neighbor_n"] = len(set(case_neighbors["neighbor"])) if not case_neighbors.empty else 0
    output["control_neighbor_n"] = len(set(control_neighbors["neighbor"])) if not control_neighbors.empty else 0
    output["shared_neighbors_removed"] = len(shared) if drop_shared_neighbors else 0
    return output


def _load_lr_network(lr_network: str | pd.DataFrame | None = None) -> tuple[pd.DataFrame, str]:
    if isinstance(lr_network, pd.DataFrame):
        required = {"ligand", "receptor"}
        missing = required.difference(lr_network.columns)
        if missing:
            raise ValueError(f"lr_network must contain columns {required}, missing {missing}")
        source = "custom"
        if "resource" in lr_network.columns:
            resources = sorted({str(value) for value in lr_network["resource"].dropna().astype(str)})
            if resources:
                source = ",".join(resources)
        return lr_network.copy(), source
    if lr_network is None or lr_network == "fallback":
        return _FALLBACK_LR_NETWORK.copy(), "fallback"
    if lr_network == "custom":
        raise ValueError("When lr_network='custom', pass a DataFrame or a file path.")

    path = Path(str(lr_network)).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Ligand-receptor network file not found: {path}")
    if path.suffix.lower() in {".tsv", ".txt"}:
        table = pd.read_csv(path, sep="\t")
    else:
        table = pd.read_csv(path)
    required = {"ligand", "receptor"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"lr_network must contain columns {required}, missing {missing}")
    return table, str(path)


def differential_lr(
    adata: AnnData,
    *,
    perturbation: str,
    control: str,
    graph_key: str | None = None,
    lr_network: str | pd.DataFrame | None = None,
    groupby: str = "perturbation",
    status_col: str = "perturbation_status",
    source_groupby: str | None = None,
    target_groupby: str | None = None,
    layer: str | None = None,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Score ligand-receptor interactions between perturbed cells and their neighbors."""
    ensure_spatialperturb_schema(adata)
    base_mask = _restrict_mask(adata, cell_type=cell_type, roi=roi)
    neighbor_mask = _restrict_mask(adata, roi=roi) if roi is not None else np.ones(adata.n_obs, dtype=bool)
    case_sources = adata.obs_names[
        _group_mask(adata, perturbation=perturbation, groupby=groupby, status_col=status_col) & base_mask
    ]
    control_sources = adata.obs_names[
        _group_mask(adata, perturbation=control, groupby=groupby, status_col=status_col) & base_mask
    ]

    case_neighbors = collect_neighbors(adata, cells=case_sources, graph_key=graph_key, exclude_perturbed=True, as_frame=True)
    control_neighbors = collect_neighbors(
        adata, cells=control_sources, graph_key=graph_key, exclude_perturbed=True, as_frame=True
    )

    allowed_neighbors = set(adata.obs_names[neighbor_mask])
    if not case_neighbors.empty:
        case_neighbors = case_neighbors[case_neighbors["neighbor"].isin(allowed_neighbors)].copy()
    if not control_neighbors.empty:
        control_neighbors = control_neighbors[control_neighbors["neighbor"].isin(allowed_neighbors)].copy()

    for column in [source_groupby, target_groupby]:
        if column is not None and column not in adata.obs.columns:
            raise KeyError(f"{column!r} not found in adata.obs.")

    network, network_source = _load_lr_network(lr_network)
    ligands = sorted(set(network["ligand"]).intersection(adata.var_names))
    receptors = sorted(set(network["receptor"]).intersection(adata.var_names))
    if not ligands or not receptors:
        return pd.DataFrame(
            columns=[
                "perturbation",
                "control",
                "network_source",
                "source_group",
                "target_group",
                "ligand",
                "receptor",
                "case_score",
                "control_score",
                "diff_score",
                "rank",
            ]
        )

    records = []
    source_groups = [None]
    if source_groupby is not None:
        source_groups = sorted(
            {
                *adata.obs.loc[list(case_sources), source_groupby].astype(str),
                *adata.obs.loc[list(control_sources), source_groupby].astype(str),
            }
        )
    target_groups = [None]
    if target_groupby is not None:
        case_neighbor_cells = case_neighbors["neighbor"].astype(str).tolist() if not case_neighbors.empty else []
        control_neighbor_cells = control_neighbors["neighbor"].astype(str).tolist() if not control_neighbors.empty else []
        target_groups = sorted(
            {
                *adata.obs.loc[list(dict.fromkeys(case_neighbor_cells)), target_groupby].astype(str),
                *adata.obs.loc[list(dict.fromkeys(control_neighbor_cells)), target_groupby].astype(str),
            }
        )

    for source_group in source_groups:
        case_source_subset = list(case_sources)
        control_source_subset = list(control_sources)
        if source_groupby is not None:
            case_source_subset = adata.obs.loc[list(case_sources)]
            case_source_subset = case_source_subset[case_source_subset[source_groupby].astype(str) == str(source_group)].index.tolist()
            control_source_subset = adata.obs.loc[list(control_sources)]
            control_source_subset = control_source_subset[
                control_source_subset[source_groupby].astype(str) == str(source_group)
            ].index.tolist()

        if not case_source_subset or not control_source_subset:
            continue

        case_edge_subset = case_neighbors[case_neighbors["source"].isin(case_source_subset)].copy()
        control_edge_subset = control_neighbors[control_neighbors["source"].isin(control_source_subset)].copy()

        for target_group in target_groups:
            case_neighbor_subset = (
                list(dict.fromkeys(case_edge_subset["neighbor"].astype(str)))
                if not case_edge_subset.empty
                else []
            )
            control_neighbor_subset = (
                list(dict.fromkeys(control_edge_subset["neighbor"].astype(str)))
                if not control_edge_subset.empty
                else []
            )
            if target_groupby is not None:
                case_neighbor_subset = (
                    adata.obs.loc[case_neighbor_subset]
                    .loc[lambda df: df[target_groupby].astype(str) == str(target_group)]
                    .index.tolist()
                )
                control_neighbor_subset = (
                    adata.obs.loc[control_neighbor_subset]
                    .loc[lambda df: df[target_groupby].astype(str) == str(target_group)]
                    .index.tolist()
                )

            if not case_neighbor_subset or not control_neighbor_subset:
                continue

            case_source_expr = pd.Series(
                extract_matrix(adata, layer=layer, obs_names=case_source_subset, var_names=ligands).mean(axis=0),
                index=ligands,
            )
            control_source_expr = pd.Series(
                extract_matrix(adata, layer=layer, obs_names=control_source_subset, var_names=ligands).mean(axis=0),
                index=ligands,
            )
            case_neighbor_expr = pd.Series(
                extract_matrix(adata, layer=layer, obs_names=case_neighbor_subset, var_names=receptors).mean(axis=0),
                index=receptors,
            )
            control_neighbor_expr = pd.Series(
                extract_matrix(adata, layer=layer, obs_names=control_neighbor_subset, var_names=receptors).mean(axis=0),
                index=receptors,
            )

            for row in network[["ligand", "receptor"]].drop_duplicates().itertuples(index=False):
                if row.ligand not in case_source_expr.index or row.receptor not in case_neighbor_expr.index:
                    continue
                case_score = float(case_source_expr[row.ligand] * case_neighbor_expr[row.receptor])
                control_score = float(control_source_expr[row.ligand] * control_neighbor_expr[row.receptor])
                records.append(
                    {
                        "perturbation": perturbation,
                        "control": control,
                        "network_source": network_source,
                        "source_group": "all" if source_group is None else str(source_group),
                        "target_group": "all" if target_group is None else str(target_group),
                        "ligand": row.ligand,
                        "receptor": row.receptor,
                        "case_source_n": len(case_source_subset),
                        "control_source_n": len(control_source_subset),
                        "case_neighbor_n": len(case_neighbor_subset),
                        "control_neighbor_n": len(control_neighbor_subset),
                        "case_score": case_score,
                        "control_score": control_score,
                        "diff_score": case_score - control_score,
                    }
                )

    output = pd.DataFrame.from_records(records)
    if output.empty:
        return output
    output = output.sort_values(
        ["source_group", "target_group", "diff_score"],
        ascending=[True, True, False],
    ).reset_index(drop=True)
    output["rank"] = np.arange(1, len(output) + 1)
    return output


def platform_concordance(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    group_col: str = "perturbation",
    gene_col: str = "gene",
    score_col: str = "log2fc",
    top_n: int = 50,
    level: str = "both",
) -> pd.DataFrame:
    """Compare per-perturbation gene effects between two result tables."""
    if level not in {"gene", "program", "both"}:
        raise ValueError("level must be one of {'gene', 'program', 'both'}.")
    required = {group_col, gene_col, score_col}
    for label, table in {"left": left, "right": right}.items():
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"{label} results are missing required columns: {missing}")

    merged = left[[group_col, gene_col, score_col]].merge(
        right[[group_col, gene_col, score_col]],
        on=[group_col, gene_col],
        suffixes=("_left", "_right"),
    )
    left_programs = derive_perturbation_programs(left, group_col=group_col, gene_col=gene_col, score_col=score_col, top_n=top_n, direction="both")
    right_programs = derive_perturbation_programs(
        right,
        group_col=group_col,
        gene_col=gene_col,
        score_col=score_col,
        top_n=top_n,
        direction="both",
    )
    program_scores = compare_program_concordance(left_programs, right_programs)
    program_scores = (
        program_scores.rename(columns={"program": group_col, "score": "program_jaccard"})
        if not program_scores.empty
        else pd.DataFrame(columns=[group_col, "program_jaccard", "left_size", "right_size"])
    )

    records = []
    for perturbation, group in merged.groupby(group_col):
        x = group[f"{score_col}_left"].to_numpy(dtype=float)
        y = group[f"{score_col}_right"].to_numpy(dtype=float)
        if len(group) > 1:
            pearson = pearsonr(x, y).statistic
            spearman = spearmanr(x, y).statistic
        else:
            pearson = np.nan
            spearman = np.nan

        left_top = set(
            left[left[group_col] == perturbation]
            .assign(_abs=lambda df: df[score_col].abs())
            .sort_values("_abs", ascending=False)
            .head(top_n)[gene_col]
            .astype(str)
        )
        right_top = set(
            right[right[group_col] == perturbation]
            .assign(_abs=lambda df: df[score_col].abs())
            .sort_values("_abs", ascending=False)
            .head(top_n)[gene_col]
            .astype(str)
        )
        union = left_top | right_top
        jaccard = len(left_top & right_top) / len(union) if union else np.nan
        program_row = program_scores[program_scores[group_col].astype(str) == str(perturbation)]
        program_jaccard = float(program_row.iloc[0]["program_jaccard"]) if not program_row.empty else np.nan
        program_left_size = int(program_row.iloc[0]["left_size"]) if not program_row.empty else 0
        program_right_size = int(program_row.iloc[0]["right_size"]) if not program_row.empty else 0
        records.append(
            {
                "perturbation": perturbation,
                "n_shared_genes": len(group),
                "pearson": pearson,
                "spearman": spearman,
                "top_gene_jaccard": jaccard,
                "program_jaccard": program_jaccard,
                "program_left_size": program_left_size,
                "program_right_size": program_right_size,
            }
        )

    output = pd.DataFrame.from_records(records).sort_values("spearman", ascending=False)
    if level == "gene":
        return output.loc[:, ["perturbation", "n_shared_genes", "pearson", "spearman", "top_gene_jaccard"]]
    if level == "program":
        return output.loc[:, ["perturbation", "program_jaccard", "program_left_size", "program_right_size"]]
    return output


def power_curve(
    adata: AnnData,
    *,
    perturbation: str,
    control: str,
    mode: str = "intrinsic",
    feature: str | None = None,
    sample_sizes: Sequence[int] = (5, 10, 20, 30),
    n_boot: int = 100,
    alpha: float = 0.05,
    groupby: str = "perturbation",
    status_col: str = "perturbation_status",
    graph_key: str | None = None,
    layer: str | None = None,
    method: str = "simple",
    sample_col: str | None = None,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
    random_state: int = 0,
) -> pd.DataFrame:
    """Estimate detection power across sample sizes by bootstrap resampling."""
    ensure_spatialperturb_schema(adata)
    rng = np.random.default_rng(random_state)
    base_mask = _restrict_mask(adata, cell_type=cell_type, roi=roi)
    neighbor_mask = _restrict_mask(adata, roi=roi) if roi is not None else np.ones(adata.n_obs, dtype=bool)

    if mode == "intrinsic":
        case_pool = adata.obs_names[
            _group_mask(adata, perturbation=perturbation, groupby=groupby, status_col=status_col) & base_mask
        ]
        control_pool = adata.obs_names[
            _group_mask(adata, perturbation=control, groupby=groupby, status_col=status_col) & base_mask
        ]
    elif mode == "neighbor":
        case_sources = adata.obs_names[
            _group_mask(adata, perturbation=perturbation, groupby=groupby, status_col=status_col) & base_mask
        ]
        control_sources = adata.obs_names[
            _group_mask(adata, perturbation=control, groupby=groupby, status_col=status_col) & base_mask
        ]
        case_edges = collect_neighbors(adata, cells=case_sources, graph_key=graph_key, exclude_perturbed=True, as_frame=True)
        control_edges = collect_neighbors(adata, cells=control_sources, graph_key=graph_key, exclude_perturbed=True, as_frame=True)
        allowed_neighbors = set(adata.obs_names[neighbor_mask])
        if not case_edges.empty:
            case_edges = case_edges[case_edges["neighbor"].isin(allowed_neighbors)].copy()
        if not control_edges.empty:
            control_edges = control_edges[control_edges["neighbor"].isin(allowed_neighbors)].copy()
        case_pool = pd.Index(sorted(set(case_edges["neighbor"])))
        control_pool = pd.Index(sorted(set(control_edges["neighbor"])))
    else:
        raise ValueError("mode must be 'intrinsic' or 'neighbor'.")

    if feature is None:
        reference = intrinsic_de(
            adata,
            perturbation=perturbation,
            control=control,
            groupby=groupby,
            status_col=status_col,
            method=method,
            sample_col=sample_col,
            cell_type=cell_type,
            roi=roi,
            effect_size_only=True,
        )
        feature = str(reference.iloc[0]["gene"])

    if feature not in adata.var_names:
        raise KeyError(f"Feature {feature!r} not found in AnnData.var_names.")

    records = []
    if method == "pseudobulk":
        if sample_col is None:
            raise ValueError("sample_col is required when method='pseudobulk'.")
        case_values = (
            pd.DataFrame(
                {
                    "sample": adata.obs.loc[list(case_pool), sample_col].astype(str).to_numpy(),
                    "value": extract_matrix(adata, layer=layer, obs_names=case_pool, var_names=[feature]).ravel(),
                }
            )
            .groupby("sample")["value"]
            .mean()
            .to_numpy(dtype=float)
        )
        control_values = (
            pd.DataFrame(
                {
                    "sample": adata.obs.loc[list(control_pool), sample_col].astype(str).to_numpy(),
                    "value": extract_matrix(adata, layer=layer, obs_names=control_pool, var_names=[feature]).ravel(),
                }
            )
            .groupby("sample")["value"]
            .mean()
            .to_numpy(dtype=float)
        )
    else:
        case_values = extract_matrix(adata, layer=layer, obs_names=case_pool, var_names=[feature]).ravel()
        control_values = extract_matrix(adata, layer=layer, obs_names=control_pool, var_names=[feature]).ravel()

    max_case = len(case_values)
    max_control = len(control_values)

    for size in sample_sizes:
        effective_n = min(int(size), max_case, max_control)
        if effective_n < 2:
            records.append(
                {
                    "perturbation": perturbation,
                    "control": control,
                    "mode": mode,
                    "feature": feature,
                    "sample_size": int(size),
                    "effective_n": effective_n,
                    "method": method,
                    "power": np.nan,
                }
            )
            continue

        hits = 0
        for _ in range(int(n_boot)):
            case_sample = rng.choice(case_values, size=effective_n, replace=False)
            control_sample = rng.choice(control_values, size=effective_n, replace=False)
            if method == "pseudobulk":
                pvalue = ttest_ind(case_sample, control_sample, equal_var=False, nan_policy="omit").pvalue
            else:
                pvalue = mannwhitneyu(case_sample, control_sample, alternative="two-sided").pvalue
            if pvalue < alpha:
                hits += 1
        records.append(
            {
                "perturbation": perturbation,
                "control": control,
                "mode": mode,
                "feature": feature,
                "sample_size": int(size),
                "effective_n": effective_n,
                "method": method,
                "power": hits / float(n_boot),
            }
        )

    return pd.DataFrame.from_records(records)
