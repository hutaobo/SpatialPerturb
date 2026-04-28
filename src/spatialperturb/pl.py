"""Plotting helpers for SpatialPerturb result tables."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData

from .gr import collect_neighbors


def barcode_spread(
    adata: AnnData,
    *,
    status_col: str = "perturbation_status",
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot perturbation assignment status counts."""
    data = adata.obs[status_col].astype(str).value_counts().rename_axis("status").reset_index(name="n_cells")
    ax = ax or plt.subplots(figsize=(6, 4))[1]
    sns.barplot(data=data, x="status", y="n_cells", hue="status", dodge=False, legend=False, ax=ax, palette="deep")
    ax.set_title("Perturbation assignment QC")
    ax.set_xlabel("")
    ax.set_ylabel("Cells")
    return ax


def own_vs_neighbor(
    intrinsic_results: pd.DataFrame,
    neighbor_results: pd.DataFrame,
    *,
    perturbation: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot intrinsic versus neighbor log fold changes for shared genes."""
    left = intrinsic_results.copy()
    right = neighbor_results.copy()
    if perturbation is not None:
        left = left[left["perturbation"] == perturbation]
        right = right[right["perturbation"] == perturbation]
    merged = left[["perturbation", "gene", "log2fc"]].merge(
        right[["perturbation", "gene", "log2fc"]],
        on=["perturbation", "gene"],
        suffixes=("_intrinsic", "_neighbor"),
    )
    ax = ax or plt.subplots(figsize=(5, 5))[1]
    sns.scatterplot(
        data=merged,
        x="log2fc_intrinsic",
        y="log2fc_neighbor",
        hue="perturbation",
        ax=ax,
        s=40,
    )
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_title("Cell-intrinsic vs neighbor effects")
    ax.set_xlabel("Intrinsic log2FC")
    ax.set_ylabel("Neighbor log2FC")
    return ax


def lr_pairs(
    lr_results: pd.DataFrame,
    *,
    top_n: int = 20,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot top ligand-receptor differential scores."""
    data = lr_results.copy().head(top_n)
    data["pair"] = data["ligand"].astype(str) + "->" + data["receptor"].astype(str)
    ax = ax or plt.subplots(figsize=(8, 4))[1]
    sns.barplot(data=data, x="pair", y="diff_score", hue="pair", dodge=False, legend=False, ax=ax, palette="viridis")
    ax.tick_params(axis="x", rotation=90)
    ax.set_title("Differential ligand-receptor pairs")
    ax.set_xlabel("")
    ax.set_ylabel("Case - control score")
    return ax


def lr_map(
    adata: AnnData,
    lr_results: pd.DataFrame,
    *,
    perturbation: str,
    graph_key: str | None = None,
    ligand: str | None = None,
    receptor: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot perturbed cells and their neighbors colored by ligand and receptor expression."""
    pair_table = lr_results[lr_results["perturbation"] == perturbation]
    if ligand is None or receptor is None:
        top_pair = pair_table.iloc[0]
        ligand = str(top_pair["ligand"])
        receptor = str(top_pair["receptor"])

    source_cells = adata.obs_names[
        (adata.obs["perturbation"].astype(str) == perturbation)
        & (adata.obs["perturbation_status"].astype(str) == "single")
    ]
    neighbor_edges = collect_neighbors(adata, cells=source_cells, graph_key=graph_key, exclude_perturbed=True, as_frame=True)
    neighbor_cells = pd.Index(sorted(set(neighbor_edges["neighbor"]))) if not neighbor_edges.empty else pd.Index([])

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    ax = ax or plt.subplots(figsize=(6, 5))[1]
    ax.scatter(coords[:, 0], coords[:, 1], color="lightgrey", s=20, alpha=0.5, label="all cells")

    if len(source_cells) > 0 and ligand in adata.var_names:
        source_idx = adata.obs_names.get_indexer(source_cells)
        ligand_values = np.asarray(adata[source_cells, [ligand]].X).ravel()
        ax.scatter(coords[source_idx, 0], coords[source_idx, 1], c=ligand_values, cmap="Reds", s=50, label=f"{ligand} in source")

    if len(neighbor_cells) > 0 and receptor in adata.var_names:
        neighbor_idx = adata.obs_names.get_indexer(neighbor_cells)
        receptor_values = np.asarray(adata[neighbor_cells, [receptor]].X).ravel()
        ax.scatter(
            coords[neighbor_idx, 0],
            coords[neighbor_idx, 1],
            c=receptor_values,
            cmap="Blues",
            s=50,
            marker="s",
            label=f"{receptor} in neighbors",
        )

    ax.set_title(f"{perturbation}: {ligand}->{receptor}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc="best")
    return ax


def platform_concordance(
    concordance_results: pd.DataFrame,
    *,
    metric: str = "spearman",
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot per-perturbation cross-platform concordance."""
    ax = ax or plt.subplots(figsize=(6, 4))[1]
    sns.barplot(
        data=concordance_results,
        x="perturbation",
        y=metric,
        hue="perturbation",
        dodge=False,
        legend=False,
        ax=ax,
        palette="crest",
    )
    ax.set_ylim(-1, 1)
    ax.set_title("Cross-platform concordance")
    ax.set_xlabel("")
    ax.set_ylabel(metric)
    return ax


def power_curve(
    power_results: pd.DataFrame,
    *,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot empirical power as a function of sample size."""
    ax = ax or plt.subplots(figsize=(6, 4))[1]
    sns.lineplot(
        data=power_results,
        x="sample_size",
        y="power",
        hue="perturbation",
        style="mode",
        marker="o",
        ax=ax,
    )
    ax.set_ylim(0, 1.05)
    ax.set_title("Power and sensitivity")
    ax.set_ylabel("Estimated power")
    ax.set_xlabel("Sample size")
    return ax
