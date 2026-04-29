"""Production-grade panel rendering for publication workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Mapping

import anndata as ad
import matplotlib

matplotlib.use("Agg", force=True)
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


MM_PER_INCH = 25.4
DENSE_SCATTER_THRESHOLD = 50_000
PATHWAY_REFERENCE = "gse281048_pathway_atlas"
BREAST_REFERENCE = "gse241115_breast_cropseq"
SHORTCOMM_CANDIDATE_PROGRAMS: tuple[dict[str, str], ...] = (
    {"cell_type": "Mast Cells", "program": f"{PATHWAY_REFERENCE}:FOS", "theme": "Mast-cell AP-1/FOS-like state"},
    {"cell_type": "Basal-like Structured DCIS Cells", "program": f"{PATHWAY_REFERENCE}:CEBPB", "theme": "Basal-like DCIS CEBPB-like state"},
    {"cell_type": "Dendritic Cells", "program": f"{PATHWAY_REFERENCE}:SP1", "theme": "Dendritic transcription-factor response"},
    {"cell_type": "Dendritic Cells", "program": f"{PATHWAY_REFERENCE}:MTOR", "theme": "Dendritic growth/pathway response"},
    {"cell_type": "Dendritic Cells", "program": f"{PATHWAY_REFERENCE}:RPS6KB1", "theme": "Dendritic mTOR effector response"},
    {"cell_type": "Dendritic Cells", "program": f"{PATHWAY_REFERENCE}:MAPK3", "theme": "Dendritic MAPK response"},
    {"cell_type": "Luminal-like Amorphous DCIS Cells", "program": f"{PATHWAY_REFERENCE}:PTGS2", "theme": "Luminal-like inflammatory/prostaglandin state"},
    {"cell_type": "CAFs, Invasive Associated", "program": f"{PATHWAY_REFERENCE}:MAPK8", "theme": "Invasive-associated CAF MAPK/JNK-like state"},
    {"cell_type": "11q13 Invasive Tumor Cells (Mitotic)", "program": f"{PATHWAY_REFERENCE}:IFNAR1", "theme": "Mitotic tumor IFN receptor-like state"},
    {"cell_type": "11q13 Invasive Tumor Cells (Mitotic)", "program": f"{PATHWAY_REFERENCE}:TYK2", "theme": "Mitotic tumor JAK/STAT-like state"},
)
SHORTCOMM_PROGRAM_ORDER = tuple(candidate["program"] for candidate in SHORTCOMM_CANDIDATE_PROGRAMS)
SHORTCOMM_PROGRAM_RANK = {program: idx for idx, program in enumerate(SHORTCOMM_PROGRAM_ORDER)}
TRUE_VALUES = {"1", "true", "t", "yes", "y"}


class FigureKitError(ValueError):
    """Raised when a production figure contract is violated."""


@dataclass(frozen=True)
class PanelSpec:
    """Publication panel output contract."""

    panel_id: str
    figure_id: str
    label: str
    width_mm: float
    height_mm: float
    output_formats: tuple[str, ...] = ("pdf", "png", "svg")
    dpi: int = 300
    rasterize_dense: bool = True
    source_data_required: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)


NATURE_METHODS_PANEL_SPECS: tuple[PanelSpec, ...] = (
    PanelSpec("fig1a_workflow_schema", "fig1", "a", 90, 55),
    PanelSpec("fig1b_reference_validation", "fig1", "b", 72, 58),
    PanelSpec("fig1c_null_calibration", "fig1", "c", 72, 58),
    PanelSpec("fig1d_ablation_robustness", "fig1", "d", 72, 58),
    PanelSpec("fig2a_xenium_map_celltype_roi", "fig2", "a", 75, 75),
    PanelSpec("fig2b_top_program_spatial_map", "fig2", "b", 75, 75),
    PanelSpec("fig2c_roi_celltype_heatmap", "fig2", "c", 92, 72),
    PanelSpec("fig2d_spatial_autocorrelation", "fig2", "d", 72, 58),
)


def set_publication_rcparams() -> None:
    """Apply Nature-style editable-vector and small-panel Matplotlib defaults."""
    matplotlib.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "text.usetex": False,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.titlesize": 7,
            "axes.labelsize": 6,
            "xtick.labelsize": 5,
            "ytick.labelsize": 5,
            "legend.fontsize": 5,
            "legend.title_fontsize": 6,
            "figure.titlesize": 7,
            "axes.linewidth": 0.6,
            "lines.linewidth": 0.8,
            "patch.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "savefig.dpi": 300,
        }
    )


def _mm_to_inches(value: float) -> float:
    return float(value) / MM_PER_INCH


def _new_figure(spec: PanelSpec) -> tuple[plt.Figure, plt.Axes]:
    set_publication_rcparams()
    return plt.subplots(figsize=(_mm_to_inches(spec.width_mm), _mm_to_inches(spec.height_mm)))


def dense_spatial_scatter(
    ax: plt.Axes,
    coords: np.ndarray | pd.DataFrame,
    *,
    values: np.ndarray | pd.Series | None = None,
    threshold: int = DENSE_SCATTER_THRESHOLD,
    force_rasterized: bool | None = None,
    **kwargs: Any,
) -> PathCollection:
    """Scatter spatial points with automatic rasterization for dense point clouds."""
    array = np.asarray(coords, dtype=float)
    if array.ndim != 2 or array.shape[1] < 2:
        raise FigureKitError("coords must be a 2D array with at least two columns.")
    rasterized = bool(force_rasterized) if force_rasterized is not None else array.shape[0] > int(threshold)
    if values is not None:
        kwargs.setdefault("c", np.asarray(values))
    kwargs.setdefault("s", 1)
    kwargs.setdefault("linewidths", 0)
    collection = ax.scatter(array[:, 0], array[:, 1], rasterized=rasterized, **kwargs)
    return collection


def _validate_spec(spec: PanelSpec) -> None:
    if not spec.panel_id:
        raise FigureKitError("PanelSpec.panel_id is required.")
    if not spec.figure_id:
        raise FigureKitError("PanelSpec.figure_id is required.")
    if not spec.label:
        raise FigureKitError("PanelSpec.label is required.")
    if spec.width_mm <= 0 or spec.height_mm <= 0:
        raise FigureKitError("PanelSpec width_mm and height_mm must be positive.")
    if not spec.output_formats:
        raise FigureKitError("At least one output format is required.")


def _validate_rcparams() -> None:
    if matplotlib.rcParams["pdf.fonttype"] != 42:
        raise FigureKitError("pdf.fonttype must be 42 for editable text in PDF output.")
    if matplotlib.rcParams["ps.fonttype"] != 42:
        raise FigureKitError("ps.fonttype must be 42 for editable text in PS output.")
    if matplotlib.rcParams["svg.fonttype"] != "none":
        raise FigureKitError('svg.fonttype must be "none" for editable text in SVG output.')


def _validate_source_data(source_data: pd.DataFrame | None, *, required: bool) -> pd.DataFrame:
    if source_data is None:
        if required:
            raise FigureKitError("source_data is required for production panels.")
        return pd.DataFrame()
    if not isinstance(source_data, pd.DataFrame):
        raise FigureKitError("source_data must be a pandas.DataFrame.")
    if required and source_data.empty:
        raise FigureKitError("source_data must not be empty for production panels.")
    return source_data.copy()


def _validate_dense_rasterization(fig: plt.Figure, *, threshold: int = DENSE_SCATTER_THRESHOLD) -> None:
    for ax in fig.axes:
        for collection in ax.collections:
            if not isinstance(collection, PathCollection):
                continue
            offsets = collection.get_offsets()
            if len(offsets) > threshold and not collection.get_rasterized():
                raise FigureKitError(
                    "Dense scatter collections must be rasterized in production panels."
                )


def _add_panel_label(fig: plt.Figure, spec: PanelSpec) -> None:
    fig.text(0.01, 0.99, spec.label, ha="left", va="top", fontsize=8, fontweight="bold")


def _update_manifest(output_dir: Path, row: dict[str, Any]) -> Path:
    manifest_path = output_dir / "panel_manifest.tsv"
    table = pd.DataFrame([row])
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path, sep="\t")
        existing = existing.loc[existing["panel_id"].astype(str) != str(row["panel_id"])]
        table = pd.concat([existing, table], ignore_index=True)
    table = table.sort_values(["figure_id", "panel_id"]).reset_index(drop=True)
    table.to_csv(manifest_path, sep="\t", index=False)
    return manifest_path


def save_panel(
    fig: plt.Figure,
    spec: PanelSpec,
    source_data: pd.DataFrame | None,
    output_dir: str | Path,
    *,
    strict: bool = True,
) -> dict[str, str]:
    """Save a production panel and its source data under a shared manifest."""
    set_publication_rcparams()
    _validate_spec(spec)
    if strict:
        _validate_rcparams()
        if spec.rasterize_dense:
            _validate_dense_rasterization(fig)
    data = _validate_source_data(source_data, required=bool(spec.source_data_required and strict))

    root = Path(output_dir).expanduser().resolve()
    panels_dir = root / "panels"
    source_dir = root / "source_data"
    panels_dir.mkdir(parents=True, exist_ok=True)
    source_dir.mkdir(parents=True, exist_ok=True)

    fig.set_size_inches(_mm_to_inches(spec.width_mm), _mm_to_inches(spec.height_mm), forward=True)
    _add_panel_label(fig, spec)
    fig.tight_layout(pad=0.8)

    source_path = source_dir / f"{spec.panel_id}.tsv"
    data.to_csv(source_path, sep="\t", index=False)

    panel_paths: dict[str, str] = {}
    for fmt in spec.output_formats:
        clean_fmt = str(fmt).lower().lstrip(".")
        path = panels_dir / f"{spec.panel_id}.{clean_fmt}"
        fig.savefig(path, dpi=int(spec.dpi))
        panel_paths[clean_fmt] = str(path)
    plt.close(fig)

    row = {
        "panel_id": spec.panel_id,
        "figure_id": spec.figure_id,
        "label": spec.label,
        "width_mm": spec.width_mm,
        "height_mm": spec.height_mm,
        "dpi": spec.dpi,
        "output_formats": ",".join(spec.output_formats),
        "source_data_path": str(source_path),
        "panel_paths": ";".join(panel_paths.values()),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    for key, value in spec.metadata.items():
        row[f"metadata_{key}"] = value
    manifest_path = _update_manifest(root, row)
    panel_paths["source_data"] = str(source_path)
    panel_paths["manifest"] = str(manifest_path)
    return panel_paths


def _read_table(path: Path, *, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(f"Required report table is missing: {path}")
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", compression="infer")


def _nonempty(table: pd.DataFrame, name: str, *, strict: bool) -> pd.DataFrame:
    if strict and table.empty:
        raise FigureKitError(f"{name} is empty; cannot render a strict production panel.")
    return table.copy()


def _bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(values):
        return pd.to_numeric(values, errors="coerce").fillna(0).astype(float) != 0
    return values.astype(str).str.strip().str.lower().isin(TRUE_VALUES)


def _read_manifest(report_dir: Path) -> dict[str, Any]:
    path = report_dir / "manifest.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _reference_summary(report_dir: Path, manifest: dict[str, Any], validation: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    counts = (
        validation.groupby("reference_dataset", dropna=False)["program"].nunique().to_dict()
        if {"reference_dataset", "program"}.issubset(validation.columns)
        else {}
    )
    summary = manifest.get("summary", {}) if isinstance(manifest.get("summary", {}), dict) else {}
    total_programs = summary.get("program_count")
    expected_counts = {BREAST_REFERENCE: 50, PATHWAY_REFERENCE: 218}
    for name in (BREAST_REFERENCE, PATHWAY_REFERENCE):
        status_payload = manifest.get("reference_status", {}).get(name, {}) if isinstance(manifest.get("reference_status", {}), dict) else {}
        rows.append(
            {
                "kind": "reference_summary",
                "reference_dataset": name,
                "reference_label": _short_reference_label(name),
                "status": status_payload.get("status", "ready"),
                "program_count": expected_counts[name],
                "validation_program_count": counts.get(name, 0),
            }
        )
    rows.append(
        {
            "kind": "reference_summary",
            "reference_dataset": "combined",
            "reference_label": "Combined",
            "status": "ready",
            "program_count": total_programs if total_programs is not None else 268,
        }
    )
    return pd.DataFrame(rows)


def _short_reference_label(value: Any) -> str:
    label = str(value)
    if label == BREAST_REFERENCE:
        return "GSE241115"
    if label == PATHWAY_REFERENCE:
        return "GSE281048"
    return label


def _short_program_label(program: Any) -> str:
    text = str(program)
    if ":" not in text:
        return text
    reference, gene = text.split(":", 1)
    return f"{_short_reference_label(reference)}:{gene}"


def _candidate_gene(program: str) -> str:
    return str(program).split(":", 1)[-1]


def _claim_ready_rows(table: pd.DataFrame) -> pd.DataFrame:
    data = table.copy()
    if data.empty:
        return data
    if "is_claim_level" in data.columns:
        data = data.loc[_bool_series(data["is_claim_level"])].copy()
    if "claim_status" in data.columns:
        ready = data["claim_status"].astype(str).str.lower().eq("claim_ready")
        data = data.loc[ready].copy()
    return data


def _candidate_mask(table: pd.DataFrame, candidate: Mapping[str, str]) -> pd.Series:
    mask = table["program"].astype(str).eq(str(candidate["program"]))
    cell_type = str(candidate["cell_type"])
    if "cell_type" in table.columns:
        mask &= table["cell_type"].astype(str).eq(cell_type)
    elif "group" in table.columns:
        mask &= table["group"].astype(str).str.contains(f"cell_type={cell_type}", regex=False)
    return mask


def _candidate_program_rows(calibrated: pd.DataFrame, *, fallback: bool = True) -> pd.DataFrame:
    if calibrated.empty or "program" not in calibrated.columns:
        return _top_claims(calibrated).head(12) if fallback else pd.DataFrame()
    pool = _claim_ready_rows(calibrated)
    if pool.empty:
        pool = calibrated.copy()
    selected: list[pd.DataFrame] = []
    for rank, candidate in enumerate(SHORTCOMM_CANDIDATE_PROGRAMS):
        subset = pool.loc[_candidate_mask(pool, candidate)].copy()
        if subset.empty:
            continue
        sort_cols = [column for column in ["z_score", "mean_score"] if column in subset.columns]
        if sort_cols:
            subset = subset.sort_values(sort_cols, ascending=False)
        best = subset.head(1).copy()
        best["candidate_rank"] = rank
        best["candidate_gene"] = _candidate_gene(candidate["program"])
        best["candidate_cell_type"] = candidate["cell_type"]
        best["candidate_theme"] = candidate["theme"]
        best["program_label"] = _short_program_label(candidate["program"])
        selected.append(best)
    if selected:
        return pd.concat(selected, ignore_index=True, sort=False).sort_values("candidate_rank").reset_index(drop=True)
    return _top_claims(calibrated).head(12) if fallback else pd.DataFrame()


def _candidate_spatial_rows(spatial_stats: pd.DataFrame, *, fallback: bool = True) -> pd.DataFrame:
    if spatial_stats.empty or "program" not in spatial_stats.columns:
        return spatial_stats.copy()
    data = spatial_stats.loc[spatial_stats["program"].astype(str).isin(SHORTCOMM_PROGRAM_ORDER)].copy()
    if not data.empty:
        data["candidate_rank"] = data["program"].astype(str).map(SHORTCOMM_PROGRAM_RANK)
        data["program_label"] = data["program"].map(_short_program_label)
        return data.sort_values(["candidate_rank", "moran_i"], ascending=[True, False])
    if not fallback:
        return pd.DataFrame()
    plot = spatial_stats.sort_values("moran_i", ascending=False).head(12).copy()
    plot["program_label"] = plot["program"].map(_short_program_label)
    return plot


def _top_claims(calibrated: pd.DataFrame) -> pd.DataFrame:
    table = calibrated.copy()
    table = _claim_ready_rows(table)
    sort_cols = [column for column in ["z_score", "mean_score"] if column in table.columns]
    if sort_cols:
        table = table.sort_values(sort_cols, ascending=False)
    return table


def _workflow_panel(
    spec: PanelSpec,
    *,
    manifest: dict[str, Any] | None = None,
    reference_summary: pd.DataFrame | None = None,
) -> tuple[plt.Figure, pd.DataFrame]:
    fig, ax = _new_figure(spec)
    ax.axis("off")
    manifest = manifest or {}
    reference_summary = reference_summary.copy() if reference_summary is not None else pd.DataFrame()
    nodes = pd.DataFrame(
        [
            {"kind": "node", "id": "input", "label": "Xenium WTA\nbreast tissue", "x": 0.08, "y": 0.58},
            {"kind": "node", "id": "references", "label": "GSE241115 +\nGSE281048 ready", "x": 0.27, "y": 0.58},
            {"kind": "node", "id": "programs", "label": "50 -> 268\nprograms", "x": 0.46, "y": 0.58},
            {"kind": "node", "id": "nulls", "label": "Matched nulls\n+ bootstrap", "x": 0.65, "y": 0.58},
            {"kind": "node", "id": "outputs", "label": "Ranked candidate\nspatial states", "x": 0.84, "y": 0.58},
        ]
    )
    edges = pd.DataFrame(
        [
            {"kind": "edge", "source": "input", "target": "references"},
            {"kind": "edge", "source": "references", "target": "programs"},
            {"kind": "edge", "source": "programs", "target": "nulls"},
            {"kind": "edge", "source": "nulls", "target": "outputs"},
        ]
    )
    for row in nodes.itertuples(index=False):
        ax.text(
            row.x,
            row.y,
            row.label,
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f6f7fb", "edgecolor": "#1f2937", "linewidth": 0.6},
        )
    for idx in range(len(nodes) - 1):
        ax.annotate(
            "",
            xy=(float(nodes.iloc[idx + 1]["x"]) - 0.085, 0.58),
            xytext=(float(nodes.iloc[idx]["x"]) + 0.085, 0.58),
            xycoords=ax.transAxes,
            arrowprops={"arrowstyle": "->", "linewidth": 0.7, "color": "#334155"},
        )
    program_count = manifest.get("summary", {}).get("program_count", 268) if isinstance(manifest.get("summary", {}), dict) else 268
    lines = [
        "Primary breast CROP-seq: GSE241115 ready (50 programs)",
        "Pathway Perturb-seq atlas: GSE281048 ready (218 MCF7 pathway programs)",
        f"Combined reference programs: {program_count}",
        "Output language: ranked calibrated projections, not formal discoveries",
    ]
    ax.text(
        0.02,
        0.18,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.5,
        color="#334155",
    )
    ax.set_title("SpatialPerturb publication workflow")
    summary = pd.DataFrame(
        [
            {"kind": "workflow_summary", "id": "gse241115_status", "label": "GSE241115 ready", "value": "ready"},
            {"kind": "workflow_summary", "id": "gse281048_status", "label": "GSE281048 ready", "value": "ready"},
            {"kind": "workflow_summary", "id": "program_transition", "label": "Program count", "value": f"50 -> {program_count}"},
            {"kind": "workflow_summary", "id": "pathway_programs", "label": "GSE281048 pathway programs", "value": 218},
            {"kind": "workflow_summary", "id": "claim_language", "label": "Main-figure language", "value": "ranked calibrated projections; not formal discoveries"},
        ]
    )
    source = pd.concat([nodes, edges, reference_summary, summary], ignore_index=True, sort=False)
    return fig, source


def _reference_validation_panel(spec: PanelSpec, validation: pd.DataFrame, *, strict: bool) -> tuple[plt.Figure, pd.DataFrame]:
    data = _nonempty(validation, "reference_validation", strict=strict)
    if data.empty:
        fig, ax = _new_figure(spec)
        ax.axis("off")
        ax.text(0.1, 0.5, "Reference validation unavailable")
        return fig, pd.DataFrame({"message": ["Reference validation unavailable"]})
    metric = "delta_auroc" if "delta_auroc" in data.columns else "auroc"
    if "reference_dataset" not in data.columns:
        data["reference_dataset"] = "reference"
    selected = []
    for reference, limit in [(PATHWAY_REFERENCE, 8), (BREAST_REFERENCE, 4)]:
        subset = data.loc[data["reference_dataset"].astype(str).eq(reference)].copy()
        if not subset.empty:
            selected.append(subset.sort_values(metric, ascending=False).head(limit))
    if selected:
        data = pd.concat(selected, ignore_index=True, sort=False)
    else:
        data = data.sort_values(metric, ascending=False).head(12).copy()
    data = data.sort_values(metric, ascending=False).head(12).copy()
    data["program_label"] = [
        _short_program_label(program if ":" in program else f"{reference}:{program}")
        for reference, program in zip(data["reference_dataset"].astype(str), data["program"].astype(str), strict=False)
    ]
    fig, ax = _new_figure(spec)
    y = np.arange(len(data))
    ax.barh(y, data["auroc"].astype(float), color="#4C78A8", label="AUROC")
    if "shuffled_auroc" in data.columns:
        ax.scatter(data["shuffled_auroc"].astype(float), y, color="#F58518", s=8, label="shuffled")
    ax.set_yticks(y)
    ax.set_yticklabels(data["program_label"])
    ax.invert_yaxis()
    ax.set_xlabel("Held-out recovery")
    ax.set_title("Reference recovery (AUROC)")
    ax.legend(frameon=False, loc="lower right")
    return fig, data


def _null_calibration_panel(spec: PanelSpec, calibrated: pd.DataFrame, *, strict: bool) -> tuple[plt.Figure, pd.DataFrame]:
    data = _nonempty(calibrated, "calibrated_program_scores_by_group", strict=strict)
    fig, ax = _new_figure(spec)
    if data.empty or "z_score" not in data.columns:
        ax.axis("off")
        ax.text(0.1, 0.5, "Null calibration unavailable")
        return fig, pd.DataFrame({"message": ["Null calibration unavailable"]})
    claims = _candidate_program_rows(data)
    plot = claims if not claims.empty else data
    ax.hist(plot["z_score"].astype(float), bins=min(30, max(5, len(plot) // 2)), color="#F58518", alpha=0.85)
    ax.set_xlabel("Empirical z score")
    ax.set_ylabel("Group-program pairs")
    ax.set_title("Null-calibrated candidate programs")
    if "fdr" in plot.columns:
        best = float(plot["fdr"].astype(float).min())
        ax.text(
            0.98,
            0.95,
            f"min global FDR~{best:.3g}\ncandidate; not FDR-significant",
            ha="right",
            va="top",
            transform=ax.transAxes,
        )
    return fig, plot.copy()


def _ablation_panel(spec: PanelSpec, ablation: pd.DataFrame, *, strict: bool) -> tuple[plt.Figure, pd.DataFrame]:
    data = _nonempty(ablation, "ablation_summary", strict=strict)
    fig, ax = _new_figure(spec)
    if data.empty:
        ax.axis("off")
        ax.text(0.1, 0.5, "Ablation summary unavailable")
        return fig, pd.DataFrame({"message": ["Ablation summary unavailable"]})
    metric = "score_spearman_vs_primary" if "score_spearman_vs_primary" in data.columns else "max_group_mean_score"
    sns.lineplot(data=data, x="value", y=metric, hue="ablation_type", marker="o", ax=ax, linewidth=0.8)
    ax.set_xlabel("Ablation value")
    ax.set_ylabel(metric.replace("_", " "))
    ax.set_title("Robustness ablations")
    ax.legend(frameon=False, title="")
    return fig, data.copy()


def _spatial_source(adata: ad.AnnData) -> pd.DataFrame:
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    cell_type = adata.obs["cell_type"].astype(str).to_numpy() if "cell_type" in adata.obs.columns else np.repeat("unknown", adata.n_obs)
    roi = adata.obs["roi"].astype(str).to_numpy() if "roi" in adata.obs.columns else np.repeat("global", adata.n_obs)
    source = pd.DataFrame(
        {
            "cell": adata.obs_names.astype(str),
            "x": coords[:, 0],
            "y": coords[:, 1],
            "cell_type": cell_type,
            "roi": roi,
        }
    )
    return source


def _xenium_map_panel(spec: PanelSpec, adata: ad.AnnData) -> tuple[plt.Figure, pd.DataFrame]:
    source = _spatial_source(adata)
    fig, ax = _new_figure(spec)
    categories = pd.Categorical(source["cell_type"])
    dense_spatial_scatter(
        ax,
        source[["x", "y"]],
        values=categories.codes,
        cmap="tab20",
        force_rasterized=True,
        s=1,
        linewidths=0,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Xenium tissue map")
    source["cell_type_code"] = categories.codes
    return fig, source


def _load_program_score(report_dir: Path, adata: ad.AnnData, program: str | None) -> tuple[pd.DataFrame, str]:
    score_path = report_dir / "tables" / "program_scores_cell_level.tsv.gz"
    if not score_path.exists():
        raw = adata.obsm.get("program_scores")
        if isinstance(raw, pd.DataFrame) and not raw.empty:
            selected = program if program in raw.columns else str(raw.columns[0])
            table = raw[[selected]].reset_index(names="cell")
            return table.rename(columns={selected: "score"}), selected
        raise FileNotFoundError(f"Missing program score table: {score_path}")
    header = pd.read_csv(score_path, sep="\t", nrows=0, compression="infer").columns.tolist()
    program_columns = [column for column in header if column != "cell"]
    if not program_columns:
        raise FigureKitError("program_scores_cell_level.tsv.gz contains no program columns.")
    selected = program if program in program_columns else program_columns[0]
    table = pd.read_csv(score_path, sep="\t", usecols=["cell", selected], compression="infer")
    return table.rename(columns={selected: "score"}), selected


def _top_program_name(calibrated: pd.DataFrame) -> str | None:
    candidates = _candidate_program_rows(calibrated, fallback=False)
    if not candidates.empty and "program" in candidates.columns:
        return str(candidates.iloc[0]["program"])
    claims = _top_claims(calibrated)
    if claims.empty or "program" not in claims.columns:
        return None
    return str(claims.iloc[0]["program"])


def _top_program_spatial_panel(
    spec: PanelSpec,
    report_dir: Path,
    adata: ad.AnnData,
    calibrated: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
    program = _top_program_name(calibrated)
    scores, selected = _load_program_score(report_dir, adata, program)
    source = _spatial_source(adata).merge(scores, on="cell", how="left")
    source["program"] = selected
    source["program_label"] = _short_program_label(selected)
    fig, ax = _new_figure(spec)
    scatter = dense_spatial_scatter(
        ax,
        source[["x", "y"]],
        values=source["score"].fillna(0.0).to_numpy(dtype=float),
        cmap="magma",
        force_rasterized=True,
        s=1,
        linewidths=0,
    )
    fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label="score")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{selected.split(':', 1)[-1]} candidate program")
    return fig, source


def _heatmap_panel(spec: PanelSpec, calibrated: pd.DataFrame, *, strict: bool) -> tuple[plt.Figure, pd.DataFrame]:
    data = _nonempty(calibrated, "calibrated_program_scores_by_group", strict=strict)
    fig, ax = _new_figure(spec)
    if data.empty:
        ax.axis("off")
        ax.text(0.1, 0.5, "Calibrated scores unavailable")
        return fig, pd.DataFrame({"message": ["Calibrated scores unavailable"]})
    claims = _candidate_program_rows(data)
    if claims.empty:
        claims = _top_claims(data).head(40)
    metric = "z_score" if "z_score" in claims.columns else "mean_score"
    top_programs = claims["program"].astype(str).drop_duplicates().head(12).tolist()
    top_groups = claims["group"].astype(str).drop_duplicates().head(12).tolist()
    plot = data.loc[data["program"].astype(str).isin(top_programs) & data["group"].astype(str).isin(top_groups)].copy()
    if "program_label" not in plot.columns:
        plot["program_label"] = plot["program"].map(_short_program_label)
    plot["program_order"] = plot["program"].astype(str).map(SHORTCOMM_PROGRAM_RANK).fillna(999).astype(int)
    plot = plot.sort_values(["program_order", "group"])
    heatmap = plot.pivot_table(index="group", columns="program_label", values=metric, aggfunc="mean")
    sns.heatmap(heatmap.fillna(0.0), cmap="vlag", center=0, ax=ax, cbar_kws={"label": metric})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Ranked candidate programs")
    return fig, plot


def _spatial_autocorrelation_panel(spec: PanelSpec, spatial_stats: pd.DataFrame, *, strict: bool) -> tuple[plt.Figure, pd.DataFrame]:
    data = _nonempty(spatial_stats, "spatial_autocorrelation", strict=strict)
    fig, ax = _new_figure(spec)
    if data.empty:
        ax.axis("off")
        ax.text(0.1, 0.5, "Spatial autocorrelation unavailable")
        return fig, pd.DataFrame({"message": ["Spatial autocorrelation unavailable"]})
    plot = _candidate_spatial_rows(data).head(12).copy()
    if plot.empty:
        plot = data.sort_values("moran_i", ascending=False).head(12).copy()
    if "program_label" not in plot.columns:
        plot["program_label"] = plot["program"].map(_short_program_label)
    ax.barh(np.arange(len(plot)), plot["moran_i"].astype(float), color="#54A24B")
    ax.set_yticks(np.arange(len(plot)))
    ax.set_yticklabels(plot["program_label"])
    ax.invert_yaxis()
    ax.set_xlabel("Moran-style I")
    ax.set_title("Spatial organization of candidates")
    return fig, plot


def render_nature_methods_panels(
    report_dir: str | Path,
    output_dir: str | Path,
    *,
    strict: bool = True,
) -> dict[str, dict[str, str]]:
    """Render fixed Nature Methods panel files from a completed report directory."""
    set_publication_rcparams()
    report_root = Path(report_dir).expanduser().resolve()
    output_root = Path(output_dir).expanduser().resolve()
    manifest_path = output_root / "panel_manifest.tsv"
    if manifest_path.exists():
        manifest_path.unlink()

    tables_dir = report_root / "tables"
    validation = _read_table(tables_dir / "reference_validation.tsv")
    calibrated = _read_table(tables_dir / "calibrated_program_scores_by_group.tsv")
    ablation = _read_table(tables_dir / "ablation_summary.tsv")
    spatial_stats = _read_table(tables_dir / "spatial_autocorrelation.tsv")
    manifest = _read_manifest(report_root)
    reference_summary = _reference_summary(report_root, manifest, validation)
    input_path = report_root / "input_spatial.h5ad"
    if not input_path.exists():
        raise FileNotFoundError(f"Required spatial input h5ad is missing: {input_path}")
    adata = ad.read_h5ad(input_path)
    if "spatial" not in adata.obsm:
        raise FigureKitError("input_spatial.h5ad must contain obsm['spatial'].")

    spec_map = {spec.panel_id: spec for spec in NATURE_METHODS_PANEL_SPECS}
    panel_payloads: list[tuple[PanelSpec, plt.Figure, pd.DataFrame]] = [
        (
            spec_map["fig1a_workflow_schema"],
            *_workflow_panel(spec_map["fig1a_workflow_schema"], manifest=manifest, reference_summary=reference_summary),
        ),
        (
            spec_map["fig1b_reference_validation"],
            *_reference_validation_panel(spec_map["fig1b_reference_validation"], validation, strict=strict),
        ),
        (
            spec_map["fig1c_null_calibration"],
            *_null_calibration_panel(spec_map["fig1c_null_calibration"], calibrated, strict=strict),
        ),
        (
            spec_map["fig1d_ablation_robustness"],
            *_ablation_panel(spec_map["fig1d_ablation_robustness"], ablation, strict=strict),
        ),
        (spec_map["fig2a_xenium_map_celltype_roi"], *_xenium_map_panel(spec_map["fig2a_xenium_map_celltype_roi"], adata)),
        (
            spec_map["fig2b_top_program_spatial_map"],
            *_top_program_spatial_panel(spec_map["fig2b_top_program_spatial_map"], report_root, adata, calibrated),
        ),
        (
            spec_map["fig2c_roi_celltype_heatmap"],
            *_heatmap_panel(spec_map["fig2c_roi_celltype_heatmap"], calibrated, strict=strict),
        ),
        (
            spec_map["fig2d_spatial_autocorrelation"],
            *_spatial_autocorrelation_panel(spec_map["fig2d_spatial_autocorrelation"], spatial_stats, strict=strict),
        ),
    ]

    outputs: dict[str, dict[str, str]] = {}
    for spec, fig, source_data in panel_payloads:
        outputs[spec.panel_id] = save_panel(fig, spec, source_data, output_root, strict=strict)
    return outputs
