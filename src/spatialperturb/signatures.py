"""Program and signature helpers for SpatialPerturb."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import spearmanr

from ._utils import extract_matrix
from .gr import collect_neighbors


def build_signature_matrix(programs: Mapping[str, Sequence[str]]) -> pd.DataFrame:
    """Convert a program-to-gene mapping into a binary membership matrix."""
    normalized = {str(name): sorted({str(gene) for gene in genes}) for name, genes in programs.items()}
    if not normalized:
        return pd.DataFrame(dtype=int)

    genes = sorted({gene for values in normalized.values() for gene in values})
    matrix = pd.DataFrame(0, index=list(normalized.keys()), columns=genes, dtype=int)
    for program, members in normalized.items():
        matrix.loc[program, members] = 1
    return matrix


def derive_perturbation_programs(
    de_results: pd.DataFrame,
    *,
    group_col: str = "perturbation",
    gene_col: str = "gene",
    score_col: str = "log2fc",
    top_n: int = 50,
    direction: str = "both",
) -> dict[str, list[str]]:
    """Derive per-perturbation programs from a tidy DE table."""
    required = {group_col, gene_col, score_col}
    missing = required.difference(de_results.columns)
    if missing:
        raise ValueError(f"de_results is missing required columns: {missing}")
    if direction not in {"up", "down", "both"}:
        raise ValueError("direction must be one of {'up', 'down', 'both'}.")

    programs: dict[str, list[str]] = {}
    for group, frame in de_results.groupby(group_col):
        subset = frame.loc[:, [gene_col, score_col]].copy()
        subset[gene_col] = subset[gene_col].astype(str)
        subset = subset.dropna(subset=[score_col])
        if direction == "up":
            subset = subset[subset[score_col] > 0].sort_values(score_col, ascending=False)
        elif direction == "down":
            subset = subset[subset[score_col] < 0].sort_values(score_col, ascending=True)
        else:
            subset = subset.assign(_abs_score=subset[score_col].abs()).sort_values("_abs_score", ascending=False)
        genes = subset[gene_col].drop_duplicates().head(int(top_n)).tolist()
        programs[str(group)] = genes
    return programs


def score_programs(
    adata: AnnData,
    programs: Mapping[str, Sequence[str]],
    *,
    layer: str | None = None,
) -> pd.DataFrame:
    """Score each program as the mean expression of its member genes per cell."""
    scores: dict[str, np.ndarray] = {}
    for program, genes in programs.items():
        valid_genes = [str(gene) for gene in genes if str(gene) in adata.var_names]
        if not valid_genes:
            scores[str(program)] = np.zeros(adata.n_obs, dtype=float)
            continue
        matrix = extract_matrix(adata, layer=layer, var_names=valid_genes)
        scores[str(program)] = matrix.mean(axis=1)
    return pd.DataFrame(scores, index=adata.obs_names.astype(str))


def _resolve_score_frame(adata: AnnData, score_key: str | pd.DataFrame) -> pd.DataFrame:
    if isinstance(score_key, pd.DataFrame):
        scores = score_key.copy()
    else:
        if score_key not in adata.obsm:
            raise KeyError(f"{score_key!r} was not found in adata.obsm.")
        raw_scores = adata.obsm[score_key]
        if isinstance(raw_scores, pd.DataFrame):
            scores = raw_scores.copy()
        else:
            scores = pd.DataFrame(raw_scores, index=adata.obs_names.astype(str))
    scores.index = adata.obs_names.astype(str)
    scores.columns = scores.columns.astype(str)
    return scores


def neighbor_program_scores(
    adata: AnnData,
    *,
    score_key: str | pd.DataFrame = "program_scores",
    graph_key: str | None = None,
    key_added: str = "neighbor_program_scores",
    exclude_perturbed: bool = True,
) -> pd.DataFrame:
    """Average program scores across each cell's neighborhood."""
    scores = _resolve_score_frame(adata, score_key)
    neighbors = collect_neighbors(adata, graph_key=graph_key, exclude_perturbed=exclude_perturbed)
    output = pd.DataFrame(0.0, index=adata.obs_names.astype(str), columns=scores.columns.astype(str))
    for cell, neighbor_names in neighbors.items():
        if not neighbor_names:
            continue
        unique_neighbors = list(dict.fromkeys(map(str, neighbor_names)))
        output.loc[str(cell)] = scores.loc[unique_neighbors].mean(axis=0).to_numpy(dtype=float)
    adata.obsm[key_added] = output
    return output


def aggregate_program_scores(
    adata: AnnData,
    score_key: str | pd.DataFrame,
    *,
    groupby: str | Sequence[str],
) -> pd.DataFrame:
    """Aggregate per-cell program scores into a tidy long table."""
    scores = _resolve_score_frame(adata, score_key)
    group_cols = [groupby] if isinstance(groupby, str) else [str(column) for column in groupby]
    missing = [column for column in group_cols if column not in adata.obs.columns]
    if missing:
        raise KeyError(f"groupby columns not found in adata.obs: {missing}")

    obs = adata.obs.loc[scores.index, group_cols].copy()
    for column in group_cols:
        obs[column] = obs[column].astype(str)
    joined = pd.concat([obs, scores], axis=1)

    records: list[dict[str, object]] = []
    for group_values, frame in joined.groupby(group_cols, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        group_label = " | ".join(f"{column}={value}" for column, value in zip(group_cols, group_values, strict=False))
        n_cells = int(len(frame))
        means = frame.loc[:, scores.columns].mean(axis=0)
        for program, mean_score in means.items():
            record: dict[str, object] = {
                "grouping": " | ".join(group_cols),
                "group": group_label,
                "program": str(program),
                "mean_score": float(mean_score),
                "n_cells": n_cells,
            }
            for column, value in zip(group_cols, group_values, strict=False):
                record[str(column)] = str(value)
            records.append(record)

    return pd.DataFrame.from_records(records).sort_values(["grouping", "group", "program"]).reset_index(drop=True)


def _context_groups(adata: AnnData, groupby: str | Sequence[str] | None) -> list[tuple[dict[str, str], AnnData]]:
    if groupby is None:
        return [({}, adata)]

    group_cols = [groupby] if isinstance(groupby, str) else [str(column) for column in groupby]
    missing = [column for column in group_cols if column not in adata.obs.columns]
    if missing:
        raise KeyError(f"groupby columns not found in adata.obs: {missing}")

    obs = adata.obs.loc[:, group_cols].copy()
    for column in group_cols:
        obs[column] = obs[column].astype(str)

    groups: list[tuple[dict[str, str], AnnData]] = []
    for group_values, frame in obs.groupby(group_cols, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        metadata = {column: str(value) for column, value in zip(group_cols, group_values, strict=False)}
        subset = adata[frame.index.astype(str)].copy()
        groups.append((metadata, subset))
    return groups


def _resolve_reference_method(
    adata: AnnData,
    *,
    perturbation: str,
    control: str,
    perturbation_col: str,
    method: str,
    sample_col: str | None,
    status_col: str,
    min_samples_per_group: int,
) -> str:
    if method in {"simple", "pseudobulk"}:
        return method
    if sample_col is None or sample_col not in adata.obs.columns:
        return "simple"

    status = adata.obs[status_col].astype(str)
    perturbations = adata.obs[perturbation_col].astype(str)
    case_samples = adata.obs.loc[(status == "single") & (perturbations == str(perturbation)), sample_col].astype(str).nunique()
    control_samples = adata.obs.loc[(status == "single") & (perturbations == str(control)), sample_col].astype(str).nunique()
    if case_samples >= min_samples_per_group and control_samples >= min_samples_per_group:
        return "pseudobulk"
    return "simple"


def _format_program_name(perturbation: str, context: Mapping[str, str], *, append_context: bool) -> str:
    if not append_context or not context:
        return str(perturbation)
    suffix = ", ".join(f"{key}={value}" for key, value in context.items())
    return f"{perturbation} | {suffix}"


def build_reference_programs(
    adata: AnnData,
    *,
    control: str = "control",
    groupby: str | Sequence[str] | None = None,
    perturbation_col: str = "perturbation",
    status_col: str = "perturbation_status",
    method: str = "auto",
    sample_col: str | None = None,
    covariates: Sequence[str] | None = None,
    top_n: int = 50,
    direction: str = "both",
    min_cells_per_group: int = 2,
    min_samples_per_group: int = 2,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
    return_de_results: bool = False,
) -> dict[str, list[str]] | tuple[dict[str, list[str]], pd.DataFrame]:
    """Build per-perturbation reference programs from a Perturb-seq AnnData object."""
    from .tl import intrinsic_de

    if perturbation_col not in adata.obs.columns:
        raise KeyError(f"{perturbation_col!r} not found in adata.obs.")
    if status_col not in adata.obs.columns:
        raise KeyError(f"{status_col!r} not found in adata.obs.")

    context_subsets = _context_groups(adata, groupby)
    append_context = len(context_subsets) > 1
    de_tables: list[pd.DataFrame] = []

    for context, subset in context_subsets:
        status = subset.obs[status_col].astype(str)
        perturbations = subset.obs.loc[status == "single", perturbation_col].astype(str)
        groups = [
            perturbation
            for perturbation in sorted(perturbations.unique())
            if perturbation not in {str(control), "unassigned", "multiple"}
        ]
        for perturbation in groups:
            resolved_method = _resolve_reference_method(
                subset,
                perturbation=perturbation,
                control=control,
                perturbation_col=perturbation_col,
                method=method,
                sample_col=sample_col,
                status_col=status_col,
                min_samples_per_group=min_samples_per_group,
            )
            try:
                result = intrinsic_de(
                    subset,
                    perturbation=perturbation,
                    control=control,
                    groupby=perturbation_col,
                    status_col=status_col,
                    method=resolved_method,
                    sample_col=sample_col if resolved_method == "pseudobulk" else None,
                    covariates=covariates,
                    min_cells_per_group=min_cells_per_group,
                    min_samples_per_group=min_samples_per_group,
                    cell_type=cell_type,
                    roi=roi,
                )
            except ValueError:
                continue
            if result.empty:
                continue
            result = result.copy()
            result["program"] = _format_program_name(perturbation, context, append_context=append_context)
            result["reference_method"] = resolved_method
            for column, value in context.items():
                result[column] = value
            de_tables.append(result)

    de_results = pd.concat(de_tables, ignore_index=True) if de_tables else pd.DataFrame()
    programs = (
        derive_perturbation_programs(
            de_results,
            group_col="program",
            gene_col="gene",
            score_col="log2fc",
            top_n=top_n,
            direction=direction,
        )
        if not de_results.empty
        else {}
    )
    if return_de_results:
        return programs, de_results
    return programs


def compare_program_concordance(
    left: Mapping[str, Sequence[str]] | pd.DataFrame,
    right: Mapping[str, Sequence[str]] | pd.DataFrame,
) -> pd.DataFrame:
    """Compare program definitions or program scores across two inputs."""
    if isinstance(left, pd.DataFrame) and isinstance(right, pd.DataFrame):
        shared = [column for column in left.columns if column in right.columns]
        records: list[dict[str, object]] = []
        for program in shared:
            x = left[program].to_numpy(dtype=float)
            y = right[program].to_numpy(dtype=float)
            score = 1.0 if np.array_equal(x, y) else spearmanr(x, y).statistic
            records.append(
                {
                    "program": str(program),
                    "score": float(score),
                    "left_size": int(np.isfinite(x).sum()),
                    "right_size": int(np.isfinite(y).sum()),
                }
            )
        return pd.DataFrame.from_records(records)

    left_sets = {str(name): {str(gene) for gene in genes} for name, genes in dict(left).items()}
    right_sets = {str(name): {str(gene) for gene in genes} for name, genes in dict(right).items()}
    shared = sorted(set(left_sets) & set(right_sets))
    records = []
    for program in shared:
        left_genes = left_sets[program]
        right_genes = right_sets[program]
        union = left_genes | right_genes
        score = len(left_genes & right_genes) / len(union) if union else 1.0
        records.append(
            {
                "program": program,
                "score": float(score),
                "left_size": len(left_genes),
                "right_size": len(right_genes),
            }
        )
    return pd.DataFrame.from_records(records)
