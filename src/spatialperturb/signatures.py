"""Program and signature helpers for SpatialPerturb."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from ._utils import benjamini_hochberg, extract_matrix


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
    key = graph_key or ("sp_knn" if "sp_knn" in adata.obsp else "sp_radius")
    if key not in adata.obsp:
        raise KeyError(f"Spatial graph {key!r} not found in adata.obsp.")

    graph = adata.obsp[key].tocsr(copy=True)
    graph.data = np.ones(graph.nnz, dtype=float)
    score_values = scores.to_numpy(dtype=float)
    if exclude_perturbed:
        allowed_mask = adata.obs["perturbation_status"].astype(str).to_numpy() == "unassigned"
        graph = graph[:, allowed_mask]
        score_values = score_values[allowed_mask, :]

    row_counts = np.asarray(graph.sum(axis=1)).ravel()
    values = graph @ score_values
    valid_rows = row_counts > 0
    values[valid_rows, :] = values[valid_rows, :] / row_counts[valid_rows, None]
    values[~valid_rows, :] = 0.0

    output = pd.DataFrame(values, index=adata.obs_names.astype(str), columns=scores.columns.astype(str))
    adata.obsm[key_added] = output
    return output


def _aggregate_scores_from_obs(obs: pd.DataFrame, scores: pd.DataFrame, group_cols: Sequence[str]) -> pd.DataFrame:
    obs = obs.loc[scores.index, list(group_cols)].copy()
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

    return _aggregate_scores_from_obs(adata.obs, scores, group_cols)


def _feature_mean_expression(adata: AnnData, *, layer: str | None = None) -> pd.Series:
    matrix = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(matrix):
        means = np.asarray(matrix.mean(axis=0)).ravel()
    else:
        means = np.asarray(matrix, dtype=float).mean(axis=0)
    return pd.Series(means, index=adata.var_names.astype(str), dtype=float)


def _expression_bins(means: pd.Series, *, n_bins: int) -> dict[str, np.ndarray]:
    ranks = means.rank(method="first")
    bins = pd.qcut(ranks, q=min(int(n_bins), len(means)), labels=False, duplicates="drop")
    frame = pd.DataFrame({"gene": means.index.astype(str), "bin": bins.astype(int).to_numpy()})
    return {gene: frame.loc[frame["bin"] == bin_id, "gene"].to_numpy(dtype=str) for gene, bin_id in zip(frame["gene"], frame["bin"], strict=False)}


def _sample_expression_matched_program(
    genes: Sequence[str],
    *,
    bins_by_gene: Mapping[str, np.ndarray],
    all_genes: np.ndarray,
    rng: np.random.Generator,
) -> list[str]:
    sampled: list[str] = []
    original = {str(gene) for gene in genes}
    for gene in genes:
        candidates = bins_by_gene.get(str(gene), all_genes)
        candidates = np.asarray([candidate for candidate in candidates if candidate not in original], dtype=str)
        if candidates.size == 0:
            candidates = all_genes
        sampled.append(str(rng.choice(candidates)))
    return sampled


def _calibrate_observed_against_nulls(observed: pd.DataFrame, nulls: pd.DataFrame, *, min_cells: int) -> pd.DataFrame:
    if observed.empty:
        return observed.copy()
    required = {"group", "program", "mean_score", "n_cells"}
    missing = required.difference(observed.columns)
    if missing:
        raise ValueError(f"observed score table is missing required columns: {sorted(missing)}")

    output = observed.copy()
    output["mean_score"] = pd.to_numeric(output["mean_score"], errors="coerce")
    output["n_cells"] = pd.to_numeric(output["n_cells"], errors="coerce").fillna(0).astype(int)
    output["is_claim_level"] = output["n_cells"] >= int(min_cells)
    output["claim_status"] = np.where(output["is_claim_level"], "claim_ready", "exploratory_small_group")

    null_stats = (
        nulls.groupby(["group", "program"], dropna=False)["mean_score"]
        .agg(null_mean="mean", null_sd="std", null_n="count")
        .reset_index()
        if not nulls.empty
        else pd.DataFrame(columns=["group", "program", "null_mean", "null_sd", "null_n"])
    )
    output = output.merge(null_stats, on=["group", "program"], how="left")
    program_stats = (
        nulls.groupby("program", dropna=False)["mean_score"].agg(program_null_mean="mean", program_null_sd="std", program_null_n="count").reset_index()
        if not nulls.empty
        else pd.DataFrame(columns=["program", "program_null_mean", "program_null_sd", "program_null_n"])
    )
    output = output.merge(program_stats, on="program", how="left")
    for target, fallback in [("null_mean", "program_null_mean"), ("null_sd", "program_null_sd"), ("null_n", "program_null_n")]:
        output[target] = output[target].fillna(output[fallback])
    output["null_sd"] = output["null_sd"].fillna(0.0).replace(0.0, np.nan)
    output["z_score"] = (output["mean_score"] - output["null_mean"]) / output["null_sd"]
    output["z_score"] = output["z_score"].replace([np.inf, -np.inf], np.nan).fillna(0.0)

    pvalues: list[float] = []
    for row in output.itertuples(index=False):
        if nulls.empty:
            pvalues.append(1.0)
            continue
        frame = nulls[(nulls["program"].astype(str) == str(row.program)) & (nulls["group"].astype(str) == str(row.group))]
        if frame.empty:
            frame = nulls[nulls["program"].astype(str) == str(row.program)]
        if frame.empty:
            pvalues.append(1.0)
            continue
        null_values = pd.to_numeric(frame["mean_score"], errors="coerce").dropna().to_numpy(dtype=float)
        pvalues.append(float((np.sum(null_values >= float(row.mean_score)) + 1.0) / (null_values.size + 1.0)))
    output["empirical_pvalue"] = pvalues
    output["fdr"] = benjamini_hochberg(output["empirical_pvalue"])
    return output.drop(columns=[col for col in ["program_null_mean", "program_null_sd", "program_null_n"] if col in output.columns])


def calibrate_program_scores(
    adata: AnnData,
    programs: Mapping[str, Sequence[str]] | None = None,
    *,
    score_key: str | pd.DataFrame = "program_scores",
    groupby: str | Sequence[str] = ("cell_type", "roi"),
    n_random: int = 100,
    n_label_shuffles: int | None = None,
    n_bins: int = 20,
    min_cells: int = 50,
    seed: int = 0,
    layer: str | None = None,
    return_nulls: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Calibrate grouped program scores with expression-matched and label-shuffle nulls.

    The returned table keeps every group, but marks groups with fewer than
    ``min_cells`` as exploratory so downstream reports do not overclaim tiny ROIs.
    """
    scores = score_programs(adata, programs, layer=layer) if programs is not None else _resolve_score_frame(adata, score_key)
    observed = aggregate_program_scores(adata, scores, groupby=groupby)
    rng = np.random.default_rng(seed)
    null_tables: list[pd.DataFrame] = []

    if programs and n_random > 0:
        means = _feature_mean_expression(adata, layer=layer)
        all_genes = means.index.to_numpy(dtype=str)
        bins_by_gene = _expression_bins(means, n_bins=n_bins)
        normalized_programs = {str(name): [str(gene) for gene in genes if str(gene) in means.index] for name, genes in programs.items()}
        for iteration in range(int(n_random)):
            random_programs = {
                program: _sample_expression_matched_program(genes, bins_by_gene=bins_by_gene, all_genes=all_genes, rng=rng)
                for program, genes in normalized_programs.items()
            }
            random_scores = score_programs(adata, random_programs, layer=layer)
            random_grouped = aggregate_program_scores(adata, random_scores, groupby=groupby)
            random_grouped["null_type"] = "expression_matched"
            random_grouped["iteration"] = iteration
            null_tables.append(random_grouped)

    label_shuffle_count = int(n_random if n_label_shuffles is None else n_label_shuffles)
    if label_shuffle_count > 0:
        group_cols = [groupby] if isinstance(groupby, str) else [str(column) for column in groupby]
        obs_group_values = adata.obs.loc[:, group_cols].copy()
        for iteration in range(label_shuffle_count):
            permutation = rng.permutation(adata.n_obs)
            shuffled_obs = obs_group_values.copy()
            for column in group_cols:
                shuffled_obs[column] = obs_group_values[column].to_numpy()[permutation]
            shuffled_grouped = _aggregate_scores_from_obs(shuffled_obs, scores, group_cols)
            shuffled_grouped["null_type"] = "group_label_shuffle"
            shuffled_grouped["iteration"] = iteration
            null_tables.append(shuffled_grouped)

    nulls = pd.concat(null_tables, ignore_index=True) if null_tables else pd.DataFrame()
    calibrated = _calibrate_observed_against_nulls(observed, nulls, min_cells=min_cells)
    return (calibrated, nulls) if return_nulls else calibrated


def bootstrap_program_score_intervals(
    adata: AnnData,
    score_key: str | pd.DataFrame = "program_scores",
    *,
    groupby: str | Sequence[str] = ("cell_type", "roi"),
    n_bootstrap: int = 200,
    min_cells: int = 50,
    seed: int = 0,
) -> pd.DataFrame:
    """Estimate bootstrap confidence intervals for grouped program means."""
    scores = _resolve_score_frame(adata, score_key)
    group_cols = [groupby] if isinstance(groupby, str) else [str(column) for column in groupby]
    obs = adata.obs.loc[scores.index, group_cols].astype(str)
    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    joined = pd.concat([obs, scores], axis=1)
    for group_values, frame in joined.groupby(group_cols, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        n_cells = int(len(frame))
        if n_cells < int(min_cells):
            continue
        group_label = " | ".join(f"{column}={value}" for column, value in zip(group_cols, group_values, strict=False))
        values = frame.loc[:, scores.columns].to_numpy(dtype=float)
        boot = np.empty((int(n_bootstrap), values.shape[1]), dtype=float)
        for iteration in range(int(n_bootstrap)):
            sampled = rng.integers(0, n_cells, size=n_cells)
            boot[iteration, :] = values[sampled, :].mean(axis=0)
        lower = np.percentile(boot, 2.5, axis=0)
        upper = np.percentile(boot, 97.5, axis=0)
        mean = values.mean(axis=0)
        for idx, program in enumerate(scores.columns):
            record: dict[str, object] = {
                "grouping": " | ".join(group_cols),
                "group": group_label,
                "program": str(program),
                "mean_score": float(mean[idx]),
                "ci_low": float(lower[idx]),
                "ci_high": float(upper[idx]),
                "n_cells": n_cells,
                "n_bootstrap": int(n_bootstrap),
            }
            for column, value in zip(group_cols, group_values, strict=False):
                record[str(column)] = str(value)
            records.append(record)
    return pd.DataFrame.from_records(records)


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
    effect_size_only: bool = False,
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
                    effect_size_only=effect_size_only,
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


def program_redundancy_table(
    programs: Mapping[str, Sequence[str]],
    *,
    threshold: float = 0.5,
) -> pd.DataFrame:
    """Return pairwise Jaccard overlap and coarse redundancy groups for programs."""
    names = list(map(str, programs.keys()))
    gene_sets = {str(name): {str(gene) for gene in genes} for name, genes in programs.items()}
    parent = {name: name for name in names}

    def find(name: str) -> str:
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    def union(left: str, right: str) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    records: list[dict[str, object]] = []
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            union_genes = gene_sets[left] | gene_sets[right]
            score = len(gene_sets[left] & gene_sets[right]) / len(union_genes) if union_genes else 1.0
            if score >= threshold:
                union(left, right)
            records.append(
                {
                    "program_a": left,
                    "program_b": right,
                    "jaccard": float(score),
                    "program_a_size": len(gene_sets[left]),
                    "program_b_size": len(gene_sets[right]),
                    "redundant": bool(score >= threshold),
                }
            )
    group_ids = {root: idx + 1 for idx, root in enumerate(sorted({find(name) for name in names}))}
    for record in records:
        record["redundancy_group_a"] = group_ids[find(str(record["program_a"]))]
        record["redundancy_group_b"] = group_ids[find(str(record["program_b"]))]
    return pd.DataFrame.from_records(records)


def _split_reference_cells(
    adata: AnnData,
    *,
    perturbation_col: str,
    status_col: str,
    control: str,
    split_strategy: str,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, str]:
    rng = np.random.default_rng(seed)
    status = adata.obs[status_col].astype(str)
    perturbations = adata.obs[perturbation_col].astype(str)
    single_mask = status == "single"

    if split_strategy == "auto":
        if "guide_id" in adata.obs.columns and adata.obs["guide_id"].astype(str).nunique() >= 4:
            split_strategy = "guide"
        elif "cell_line" in adata.obs.columns and adata.obs["cell_line"].astype(str).nunique() >= 2:
            split_strategy = "cell_line"
        else:
            split_strategy = "cell"

    test_mask = np.zeros(adata.n_obs, dtype=bool)
    if split_strategy in {"guide", "cell_line"}:
        split_col = "guide_id" if split_strategy == "guide" else "cell_line"
        for perturbation in sorted(perturbations[single_mask].unique()):
            candidates = adata.obs.loc[single_mask & (perturbations == perturbation), split_col].astype(str).dropna().unique()
            if len(candidates) >= 2 or perturbation == str(control):
                test_value = str(rng.choice(candidates))
                test_mask |= ((adata.obs[split_col].astype(str) == test_value) & (perturbations == perturbation) & single_mask).to_numpy()
    elif split_strategy == "cell":
        for perturbation in sorted(perturbations[single_mask].unique()):
            idx = np.where(((perturbations == perturbation) & single_mask).to_numpy())[0]
            if idx.size < 2:
                continue
            n_test = max(1, int(np.ceil(idx.size * 0.25)))
            test_mask[rng.choice(idx, size=n_test, replace=False)] = True
    else:
        raise ValueError("split_strategy must be one of {'auto', 'guide', 'cell_line', 'cell'}.")

    train_mask = single_mask.to_numpy() & ~test_mask
    if not test_mask.any() or not train_mask.any():
        raise ValueError("Could not create a non-empty train/test split for reference validation.")
    return train_mask, test_mask, split_strategy


def validate_reference_programs(
    reference_adata: AnnData,
    *,
    programs: Mapping[str, Sequence[str]] | None = None,
    query_adata: AnnData | None = None,
    control: str = "control",
    perturbation_col: str = "perturbation",
    status_col: str = "perturbation_status",
    split_strategy: str = "auto",
    top_n: int = 50,
    direction: str = "both",
    method: str = "simple",
    seed: int = 0,
    min_cells_per_group: int = 2,
) -> pd.DataFrame:
    """Validate reference programs by held-out recovery of perturbation labels."""
    if perturbation_col not in reference_adata.obs.columns:
        raise KeyError(f"{perturbation_col!r} not found in reference_adata.obs.")
    if status_col not in reference_adata.obs.columns:
        raise KeyError(f"{status_col!r} not found in reference_adata.obs.")

    train_mask, test_mask, resolved_strategy = _split_reference_cells(
        reference_adata,
        perturbation_col=perturbation_col,
        status_col=status_col,
        control=control,
        split_strategy=split_strategy,
        seed=seed,
    )
    if programs is None:
        programs = build_reference_programs(
            reference_adata[train_mask].copy(),
            control=control,
            perturbation_col=perturbation_col,
            status_col=status_col,
            method=method,
            top_n=top_n,
            direction=direction,
            min_cells_per_group=min_cells_per_group,
            effect_size_only=True,
        )
    programs = {str(name): [str(gene) for gene in genes if str(gene) in reference_adata.var_names] for name, genes in programs.items()}
    test = reference_adata[test_mask].copy()
    scores = score_programs(test, programs)
    labels = test.obs[perturbation_col].astype(str)
    rng = np.random.default_rng(seed + 1)

    records: list[dict[str, object]] = []
    query_genes = set(map(str, query_adata.var_names)) if query_adata is not None else None
    for program, genes in programs.items():
        plain_program = str(program).split(":", 1)[-1].split(" | ", 1)[0]
        y_true = (labels == plain_program).to_numpy(dtype=int)
        score_values = scores[str(program)].to_numpy(dtype=float)
        positive_n = int(y_true.sum())
        negative_n = int((1 - y_true).sum())
        if positive_n > 0 and negative_n > 0:
            auroc = float(roc_auc_score(y_true, score_values))
            auprc = float(average_precision_score(y_true, score_values))
            shuffled = rng.permutation(y_true)
            shuffled_auroc = float(roc_auc_score(shuffled, score_values)) if len(np.unique(shuffled)) == 2 else np.nan
        else:
            auroc = np.nan
            auprc = np.nan
            shuffled_auroc = np.nan
        covered = sum(1 for gene in genes if query_genes is not None and gene in query_genes)
        records.append(
            {
                "program": str(program),
                "split_strategy": resolved_strategy,
                "train_n": int(train_mask.sum()),
                "test_n": int(test_mask.sum()),
                "positive_n": positive_n,
                "negative_n": negative_n,
                "auroc": auroc,
                "auprc": auprc,
                "shuffled_auroc": shuffled_auroc,
                "delta_auroc": auroc - shuffled_auroc if np.isfinite(auroc) and np.isfinite(shuffled_auroc) else np.nan,
                "program_gene_count": len(genes),
                "query_gene_coverage": covered / len(genes) if query_genes is not None and genes else np.nan,
            }
        )
    return pd.DataFrame.from_records(records).sort_values("auroc", ascending=False, na_position="last").reset_index(drop=True)


def spatial_autocorrelation_scores(
    adata: AnnData,
    *,
    score_key: str | pd.DataFrame = "program_scores",
    graph_key: str | None = None,
    n_permutations: int = 100,
    seed: int = 0,
) -> pd.DataFrame:
    """Compute graph Moran-style spatial autocorrelation for program scores."""
    scores = _resolve_score_frame(adata, score_key)
    key = graph_key or ("sp_knn" if "sp_knn" in adata.obsp else "sp_radius")
    if key not in adata.obsp:
        raise KeyError(f"Spatial graph {key!r} not found in adata.obsp.")
    graph = adata.obsp[key].tocsr(copy=True)
    graph.data = np.ones(graph.nnz, dtype=float)
    graph.setdiag(0)
    graph.eliminate_zeros()
    n = graph.shape[0]
    weight_sum = float(graph.sum())
    if weight_sum <= 0:
        raise ValueError(f"Spatial graph {key!r} has no edges.")

    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    for program in scores.columns:
        x = scores[str(program)].to_numpy(dtype=float)
        centered = x - np.nanmean(x)
        denom = float(np.dot(centered, centered))
        observed = 0.0 if denom == 0 else float((n / weight_sum) * (centered @ (graph @ centered)) / denom)
        null_values = np.zeros(int(n_permutations), dtype=float)
        for iteration in range(int(n_permutations)):
            permuted = rng.permutation(centered)
            perm_denom = float(np.dot(permuted, permuted))
            null_values[iteration] = 0.0 if perm_denom == 0 else float((n / weight_sum) * (permuted @ (graph @ permuted)) / perm_denom)
        null_mean = float(null_values.mean()) if null_values.size else 0.0
        null_sd = float(null_values.std(ddof=1)) if null_values.size > 1 else 0.0
        pvalue = float((np.sum(np.abs(null_values) >= abs(observed)) + 1.0) / (null_values.size + 1.0)) if null_values.size else 1.0
        records.append(
            {
                "program": str(program),
                "graph_key": key,
                "moran_i": observed,
                "null_mean": null_mean,
                "null_sd": null_sd,
                "z_score": (observed - null_mean) / null_sd if null_sd > 0 else 0.0,
                "pvalue": pvalue,
                "n_permutations": int(n_permutations),
            }
        )
    output = pd.DataFrame.from_records(records)
    output["fdr"] = benjamini_hochberg(output["pvalue"]) if not output.empty else []
    return output.sort_values("moran_i", ascending=False).reset_index(drop=True)
