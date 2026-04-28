"""Paper-grade reporting and figure rendering helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

from . import pl


def _workflow_schema_figure(ax: plt.Axes) -> plt.Axes:
    steps = [
        ("Input", "Prepared AnnData\nperturbations + spatial graph"),
        ("QC", "Assignment QC\nsingle vs multiple"),
        ("Effects", "Intrinsic + neighbor\nDE and power"),
        ("Communication", "Ligand-receptor\ndifferential scoring"),
        ("Concordance", "Cross-platform\ngene/program agreement"),
        ("Output", "Main figures\nand tidy tables"),
    ]
    ax.axis("off")
    y = 0.5
    x_positions = [0.08, 0.24, 0.40, 0.58, 0.76, 0.92]
    for idx, ((title, body), x) in enumerate(zip(steps, x_positions, strict=False)):
        ax.text(
            x,
            y,
            f"{title}\n{body}",
            ha="center",
            va="center",
            fontsize=11,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "#f6f7fb", "edgecolor": "#1f2937", "linewidth": 1.0},
            transform=ax.transAxes,
        )
        if idx < len(steps) - 1:
            ax.annotate(
                "",
                xy=(x_positions[idx + 1] - 0.06, y),
                xytext=(x + 0.06, y),
                xycoords=ax.transAxes,
                arrowprops={"arrowstyle": "->", "color": "#334155", "linewidth": 1.5},
            )
    ax.set_title("SpatialPerturb workflow and schema")
    return ax


def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def render_paper_figures(results: dict[str, Any], *, output_dir: str | Path) -> dict[str, str]:
    """Render the canonical paper figures from benchmark result tables."""
    output_root = Path(output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    figure_paths: dict[str, str] = {}

    fig, ax = plt.subplots(figsize=(14, 3.5))
    _workflow_schema_figure(ax)
    path = output_root / "workflow_schema.png"
    _save_figure(fig, path)
    figure_paths["workflow_schema"] = str(path)

    adata = results.get("adata")
    if adata is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        pl.barcode_spread(adata, ax=ax)
        path = output_root / "assignment_qc.png"
        _save_figure(fig, path)
        figure_paths["assignment_qc"] = str(path)

    intrinsic = results.get("intrinsic_de")
    neighbor = results.get("neighbor_de")
    if isinstance(intrinsic, pd.DataFrame) and isinstance(neighbor, pd.DataFrame) and not intrinsic.empty and not neighbor.empty:
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        perturbation = str(intrinsic.iloc[0]["perturbation"]) if "perturbation" in intrinsic.columns else None
        pl.own_vs_neighbor(intrinsic, neighbor, perturbation=perturbation, ax=ax)
        path = output_root / "own_vs_neighbor.png"
        _save_figure(fig, path)
        figure_paths["own_vs_neighbor"] = str(path)

    lr = results.get("differential_lr")
    if isinstance(lr, pd.DataFrame) and not lr.empty:
        fig, ax = plt.subplots(figsize=(9, 4.5))
        pl.lr_pairs(lr, ax=ax)
        path = output_root / "lr_differential.png"
        _save_figure(fig, path)
        figure_paths["lr_differential"] = str(path)

    concordance = results.get("platform_concordance")
    if isinstance(concordance, pd.DataFrame) and not concordance.empty:
        metric = "spearman" if "spearman" in concordance.columns else concordance.columns[-1]
        fig, ax = plt.subplots(figsize=(6.5, 4))
        pl.platform_concordance(concordance, metric=metric, ax=ax)
        path = output_root / "platform_concordance.png"
        _save_figure(fig, path)
        figure_paths["platform_concordance"] = str(path)

    power = results.get("power_curve")
    if isinstance(power, pd.DataFrame) and not power.empty:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        pl.power_curve(power, ax=ax)
        path = output_root / "power_curve.png"
        _save_figure(fig, path)
        figure_paths["power_curve"] = str(path)

    manifest_path = output_root / "figures_manifest.json"
    manifest_path.write_text(json.dumps(figure_paths, indent=2), encoding="utf-8")
    figure_paths["manifest"] = str(manifest_path)
    return figure_paths
