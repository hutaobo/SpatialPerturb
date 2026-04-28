"""Benchmark orchestration for SpatialPerturb analyses."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy.stats import spearmanr

from .datasets import available_datasets, load_public_dataset
from .gr import build_spatial_graph
from .io import read_xenium
from .reports import render_paper_figures
from .schema import ensure_spatialperturb_schema
from .signatures import (
    aggregate_program_scores,
    bootstrap_program_score_intervals,
    build_reference_programs,
    build_signature_matrix,
    calibrate_program_scores,
    neighbor_program_scores,
    program_redundancy_table,
    score_programs,
    spatial_autocorrelation_scores,
    validate_reference_programs,
)
from .tl import differential_lr, intrinsic_de, neighbor_de, platform_concordance, power_curve

_BENCHMARK_CATALOG = pd.DataFrame(
    [
        {
            "benchmark": "shen_2026_core",
            "description": "Reproduce intrinsic, neighbor, ligand-receptor, power, and figure outputs on a spatial perturbation dataset.",
            "required_inputs": "Prepared spatial AnnData with perturbation assignments and spatial coordinates.",
        },
        {
            "benchmark": "cross_platform_concordance",
            "description": "Compare perturbation signatures between spatial and dissociated reference datasets.",
            "required_inputs": "Two tidy DE result tables aligned on perturbation and gene.",
        },
        {
            "benchmark": "reference_projection",
            "description": "Project reference-derived programs onto a spatial or Xenium dataset and summarize neighborhood context.",
            "required_inputs": "Prepared spatial AnnData plus one or more reference AnnData objects or registered datasets.",
        },
        {
            "benchmark": "breast_reference_projection",
            "description": "Project breast Perturb-seq reference programs onto Xenium WTA breast tissue.",
            "required_inputs": "Xenium WTA AnnData plus GSE241115 and optionally GSE281048 prepared reference datasets.",
        },
        {
            "benchmark": "nature_methods_breast_shortcomm",
            "description": "Publication-grade breast Xenium reference projection with null calibration, validation, spatial statistics, and paper-ready figures.",
            "required_inputs": "Prepared Xenium WTA AnnData plus GSE241115 and, when available, GSE281048 MCF7 pathway reference.",
        },
    ]
)


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, AnnData):
        return {
            "type": "AnnData",
            "n_obs": int(value.n_obs),
            "n_vars": int(value.n_vars),
            "dataset": value.uns.get("spatialperturb", {}).get("dataset_name"),
        }
    if isinstance(value, pd.DataFrame):
        return {"type": "DataFrame", "rows": int(len(value)), "columns": list(map(str, value.columns))}
    if isinstance(value, pd.Series):
        return {"type": "Series", "length": int(len(value))}
    return value


def _concat_or_empty(tables: Sequence[pd.DataFrame]) -> pd.DataFrame:
    frames = [table for table in tables if table is not None and not table.empty]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _package_version() -> str:
    try:
        from importlib.metadata import version

        return version("SpatialPerturb")
    except Exception:
        return "local-dev"


def _infer_perturbations(adata: AnnData, control: str) -> list[str]:
    status = adata.obs["perturbation_status"].astype(str)
    perturbations = adata.obs.loc[status == "single", "perturbation"].astype(str)
    excluded = {str(control), "unassigned", "multiple"}
    return sorted([name for name in perturbations.unique() if name not in excluded])


def _write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def _write_frame(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def _write_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=_json_default), encoding="utf-8")


def _write_program_heatmap(grouped_scores: pd.DataFrame, path: Path, *, title: str) -> str | None:
    if grouped_scores.empty:
        return None
    heatmap_data = grouped_scores.pivot_table(index="group", columns="program", values="mean_score", aggfunc="mean").fillna(0.0)
    if heatmap_data.empty:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    fig_width = max(6.0, 0.45 * heatmap_data.shape[1] + 2.5)
    fig_height = max(4.0, 0.4 * heatmap_data.shape[0] + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(heatmap_data, cmap="viridis", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Program")
    ax.set_ylabel("Group")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def available_benchmarks() -> pd.DataFrame:
    """Return available benchmark tracks and their required inputs."""
    return _BENCHMARK_CATALOG.copy()


def run_cross_platform_benchmark(
    spatial_results: pd.DataFrame,
    reference_results: pd.DataFrame,
    *,
    config: Mapping[str, Any] | None = None,
    output_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Run the cross-platform concordance benchmark on paired DE result tables."""
    cfg = dict(config or {})
    concordance = platform_concordance(
        spatial_results,
        reference_results,
        top_n=int(cfg.get("top_n", 50)),
        level=str(cfg.get("level", "both")),
    )
    if output_dir is not None:
        output_root = Path(output_dir).expanduser().resolve()
        _write_table(concordance, output_root / "tables" / "platform_concordance.tsv")
        _write_manifest(
            output_root / "manifest.json",
            {
                "benchmark": "cross_platform_concordance",
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "package_version": _package_version(),
                "config": cfg,
                "summary": {"rows": len(concordance)},
            },
        )
    return concordance


def run_core_benchmark(
    dataset_or_adata: str | AnnData,
    *,
    perturbations: Sequence[str] | None = None,
    control: str | None = None,
    target_map: Mapping[str, str] | None = None,
    lr_network: str | pd.DataFrame | None = None,
    graph_key: str | None = None,
    sample_sizes: Sequence[int] | None = None,
    cell_type: str | Sequence[str] | None = None,
    roi: str | Sequence[str] | None = None,
    config: Mapping[str, Any] | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Run the core benchmark suite for a spatial perturbation dataset."""
    cfg = dict(config or {})
    cache_dir = Path(cfg.get("cache_dir", ".spatialperturb-cache")).expanduser().resolve()

    if isinstance(dataset_or_adata, AnnData):
        adata = ensure_spatialperturb_schema(dataset_or_adata.copy())
        dataset_name = str(adata.uns.get("spatialperturb", {}).get("dataset_name", "custom_dataset"))
    else:
        dataset_name = str(dataset_or_adata)
        adata = load_public_dataset(dataset_name, cache_dir=cache_dir)

    if perturbations is not None:
        cfg["perturbations"] = list(perturbations)
    if control is not None:
        cfg["control"] = control
    if target_map is not None:
        cfg["target_map"] = dict(target_map)
    if lr_network is not None:
        cfg["lr_network"] = lr_network
    if graph_key is not None:
        cfg["graph_key"] = graph_key
    if sample_sizes is not None:
        cfg["sample_sizes"] = list(sample_sizes)
    if cell_type is not None:
        cfg["cell_type"] = cell_type
    if roi is not None:
        cfg["roi"] = roi

    graph_key = cfg.get("graph_key")
    if graph_key is None:
        graph_key = "sp_knn" if "sp_knn" in adata.obsp else ("sp_radius" if "sp_radius" in adata.obsp else None)
    if graph_key is None and "spatial" in adata.obsm:
        build_spatial_graph(
            adata,
            mode=str(cfg.get("graph_mode", "knn")),
            k=int(cfg.get("k", 15)),
            radius=cfg.get("radius"),
        )
        graph_key = "sp_knn" if "sp_knn" in adata.obsp else "sp_radius"

    control = str(cfg.get("control", "control"))
    perturbations = list(cfg.get("perturbations") or _infer_perturbations(adata, control))
    target_map = dict(cfg.get("target_map") or {})
    sample_col = cfg.get("sample_col")
    if sample_col is None and "sample" in adata.obs.columns:
        sample_col = "sample"
    method = str(cfg.get("method", "pseudobulk" if sample_col is not None else "simple"))
    cell_type = cfg.get("cell_type")
    roi = cfg.get("roi")
    lr_network = cfg.get("lr_network", "fallback")
    sample_sizes = tuple(cfg.get("sample_sizes", (5, 10, 20)))
    covariates = cfg.get("covariates")
    min_cells_per_group = int(cfg.get("min_cells_per_group", 2))
    min_samples_per_group = int(cfg.get("min_samples_per_group", 2))

    intrinsic_tables: list[pd.DataFrame] = []
    neighbor_tables: list[pd.DataFrame] = []
    lr_tables: list[pd.DataFrame] = []
    power_tables: list[pd.DataFrame] = []

    for perturbation in perturbations:
        intrinsic_tables.append(
            intrinsic_de(
                adata,
                perturbation=perturbation,
                control=control,
                cell_type=cell_type,
                roi=roi,
                method=method,
                sample_col=sample_col,
                covariates=covariates,
                min_cells_per_group=min_cells_per_group,
                min_samples_per_group=min_samples_per_group,
            )
        )
        if graph_key is not None:
            neighbor_tables.append(
                neighbor_de(
                    adata,
                    perturbation=perturbation,
                    control=control,
                    graph_key=graph_key,
                    cell_type=cell_type,
                    roi=roi,
                    method=method,
                    sample_col=sample_col,
                    covariates=covariates,
                    aggregate=str(cfg.get("neighbor_aggregate", "mean")),
                    weight_by_distance=bool(cfg.get("weight_by_distance", False)),
                    drop_shared_neighbors=bool(cfg.get("drop_shared_neighbors", False)),
                    min_cells_per_group=min_cells_per_group,
                    min_samples_per_group=min_samples_per_group,
                )
            )
            lr_tables.append(
                differential_lr(
                    adata,
                    perturbation=perturbation,
                    control=control,
                    graph_key=graph_key,
                    lr_network=lr_network,
                    source_groupby=cfg.get("source_groupby"),
                    target_groupby=cfg.get("target_groupby"),
                    cell_type=cell_type,
                    roi=roi,
                )
            )
        feature = target_map.get(perturbation)
        power_tables.append(
            power_curve(
                adata,
                perturbation=perturbation,
                control=control,
                feature=feature,
                sample_sizes=sample_sizes,
                graph_key=graph_key,
                method=method,
                sample_col=sample_col,
                cell_type=cell_type,
                roi=roi,
                n_boot=int(cfg.get("n_boot", 100)),
                alpha=float(cfg.get("alpha", 0.05)),
            )
        )

    results: dict[str, Any] = {
        "intrinsic_de": _concat_or_empty(intrinsic_tables),
        "neighbor_de": _concat_or_empty(neighbor_tables),
        "differential_lr": _concat_or_empty(lr_tables),
        "power_curve": _concat_or_empty(power_tables),
        "dataset_catalog": available_datasets(),
    }

    reference_results = cfg.get("reference_results")
    reference_input = cfg.get("reference_dataset") or cfg.get("reference_adata")
    if reference_results is None and reference_input is not None:
        if isinstance(reference_input, AnnData):
            reference_adata = ensure_spatialperturb_schema(reference_input.copy())
        else:
            reference_adata = load_public_dataset(str(reference_input), cache_dir=cache_dir)
        reference_tables = [
            intrinsic_de(
                reference_adata,
                perturbation=perturbation,
                control=control,
                cell_type=cell_type,
                roi=roi,
                method=method if sample_col is not None and "sample" in reference_adata.obs.columns else "simple",
                sample_col=sample_col if sample_col in reference_adata.obs.columns else None,
                covariates=covariates,
                min_cells_per_group=min_cells_per_group,
                min_samples_per_group=min_samples_per_group,
            )
            for perturbation in perturbations
        ]
        reference_results = _concat_or_empty(reference_tables)

    if isinstance(reference_results, pd.DataFrame) and not reference_results.empty and not results["intrinsic_de"].empty:
        results["platform_concordance"] = run_cross_platform_benchmark(
            results["intrinsic_de"],
            reference_results,
            config={"top_n": cfg.get("top_n", 50), "level": cfg.get("concordance_level", "both")},
        )

    if output_dir is not None:
        output_root = Path(output_dir).expanduser().resolve()
        tables_dir = output_root / "tables"
        figures_dir = output_root / "figures"
        output_root.mkdir(parents=True, exist_ok=True)

        table_keys = ["intrinsic_de", "neighbor_de", "differential_lr", "power_curve", "platform_concordance"]
        for key in table_keys:
            table = results.get(key)
            if isinstance(table, pd.DataFrame) and not table.empty:
                _write_table(table, tables_dir / f"{key}.tsv")
        prepared_input = output_root / "input.h5ad"
        adata.write_h5ad(prepared_input)
        report_results = {"dataset": dataset_name, "adata": adata, **results}
        figure_paths = render_paper_figures(report_results, output_dir=figures_dir)

        manifest = {
            "benchmark": "shen_2026_core",
            "dataset": dataset_name,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "package_version": _package_version(),
            "config": cfg,
            "tables": {key: str(tables_dir / f"{key}.tsv") for key in table_keys if (tables_dir / f"{key}.tsv").exists()},
            "figures": figure_paths,
            "input_h5ad": str(prepared_input),
            "summary": {
                "perturbations": perturbations,
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "method": method,
                "sample_col": sample_col,
            },
        }
        _write_manifest(output_root / "manifest.json", manifest)
        _write_manifest(output_root / "config.json", cfg)
        results["report_dir"] = str(output_root)
        results["manifest"] = manifest
    elif bool(cfg.get("include_input", False)):
        results["dataset"] = dataset_name
        results["adata"] = adata

    return results


def run_reference_projection_benchmark(
    spatial_input: str | Path | AnnData,
    *,
    reference_datasets: Sequence[str],
    config: Mapping[str, Any] | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Project reference-derived programs onto a spatial or Xenium dataset."""
    cfg = dict(config or {})
    cache_dir = Path(cfg.get("cache_dir", ".spatialperturb-cache")).expanduser().resolve()

    if isinstance(spatial_input, AnnData):
        spatial_adata = ensure_spatialperturb_schema(spatial_input.copy())
        dataset_name = str(spatial_adata.uns.get("spatialperturb", {}).get("dataset_name", "custom_xenium"))
    else:
        spatial_path = Path(spatial_input).expanduser()
        if spatial_path.suffix == ".h5ad":
            spatial_adata = ad.read_h5ad(spatial_path)
            ensure_spatialperturb_schema(spatial_adata, metadata={"platform": "xenium", "source_path": str(spatial_path)})
        else:
            spatial_adata = read_xenium(
                spatial_path,
                cell_group_path=cfg.get("cell_group_path"),
                roi_geojson_path=cfg.get("roi_geojson_path"),
                sample_name=cfg.get("sample_name"),
                load_molecules=bool(cfg.get("load_molecules", False)),
            )
        dataset_name = spatial_path.stem

    if "sp_knn" not in spatial_adata.obsp and "sp_radius" not in spatial_adata.obsp:
        build_spatial_graph(spatial_adata, mode="knn", k=int(cfg.get("k", 15)))

    reference_adatas = dict(cfg.get("reference_adatas") or {})
    reference_objects: dict[str, AnnData] = {}
    reference_rows: list[dict[str, Any]] = []
    all_programs: dict[str, list[str]] = {}
    reference_de_tables: list[pd.DataFrame] = []

    for dataset_name_ref in reference_datasets:
        if str(dataset_name_ref) in reference_adatas:
            reference = ensure_spatialperturb_schema(reference_adatas[str(dataset_name_ref)].copy())
        else:
            reference = load_public_dataset(str(dataset_name_ref), cache_dir=cache_dir)
        if str(dataset_name_ref) == "gse281048_pathway_atlas":
            pathway_cell_line = str(cfg.get("pathway_cell_line", "MCF7"))
            if "cell_line" not in reference.obs.columns:
                raise KeyError("gse281048_pathway_atlas requires obs['cell_line'] for MCF7 filtering.")
            reference = reference[reference.obs["cell_line"].astype(str) == pathway_cell_line].copy()
            if reference.n_obs == 0:
                raise ValueError(f"No cells remain after filtering gse281048_pathway_atlas to {pathway_cell_line!r}.")
        reference_objects[str(dataset_name_ref)] = reference

        programs, de_results = build_reference_programs(
            reference,
            control=str(cfg.get("reference_control", "control")),
            groupby=cfg.get("default_reference_groupby"),
            method=str(cfg.get("reference_method", "auto")),
            sample_col=cfg.get("reference_sample_col", "sample" if "sample" in reference.obs.columns else None),
            covariates=cfg.get("reference_covariates"),
            top_n=int(cfg.get("top_n", 50)),
            direction=str(cfg.get("direction", "both")),
            effect_size_only=bool(cfg.get("reference_effect_size_only", False)),
            return_de_results=True,
        )
        all_programs.update({f"{dataset_name_ref}:{program}": genes for program, genes in programs.items()})
        if not de_results.empty:
            de_results = de_results.copy()
            de_results["reference_dataset"] = str(dataset_name_ref)
            reference_de_tables.append(de_results)
        reference_rows.append(
            {
                "dataset": str(dataset_name_ref),
                "n_obs": int(reference.n_obs),
                "n_vars": int(reference.n_vars),
                "platform": str(reference.uns.get("spatialperturb", {}).get("platform", "unknown")),
                "n_programs": int(len(programs)),
            }
        )

    reference_summary = pd.DataFrame(reference_rows)
    program_scores = score_programs(spatial_adata, all_programs) if all_programs else pd.DataFrame(index=spatial_adata.obs_names.astype(str))
    spatial_adata.obsm["program_scores"] = program_scores
    neighborhood_scores = neighbor_program_scores(spatial_adata, score_key="program_scores")
    grouped_scores = aggregate_program_scores(spatial_adata, "program_scores", groupby=cfg.get("groupby", ["cell_type", "roi"]))
    grouped_neighbor_scores = aggregate_program_scores(spatial_adata, neighborhood_scores, groupby=cfg.get("groupby", ["cell_type", "roi"]))
    reference_de = _concat_or_empty(reference_de_tables)
    program_membership = build_signature_matrix(all_programs) if all_programs else pd.DataFrame()

    results: dict[str, Any] = {
        "adata": spatial_adata,
        "program_scores": program_scores,
        "program_scores_by_group": grouped_scores,
        "neighbor_program_scores": neighborhood_scores,
        "neighbor_program_scores_by_group": grouped_neighbor_scores,
        "reference_de": reference_de,
        "reference_program_membership": program_membership,
        "reference_summary": reference_summary,
        "dataset_catalog": available_datasets(),
    }

    if output_dir is not None:
        output_root = Path(output_dir).expanduser().resolve()
        tables_dir = output_root / "tables"
        figures_dir = output_root / "figures"
        references_dir = output_root / "references"
        output_root.mkdir(parents=True, exist_ok=True)

        spatial_adata.write_h5ad(output_root / "input_spatial.h5ad")
        if not reference_summary.empty:
            _write_table(reference_summary, tables_dir / "reference_summary.tsv")
        if not grouped_scores.empty:
            _write_table(grouped_scores, tables_dir / "program_scores_by_group.tsv")
        if not grouped_neighbor_scores.empty:
            _write_table(grouped_neighbor_scores, tables_dir / "neighbor_program_scores_by_group.tsv")
        if not program_membership.empty:
            _write_frame(program_membership.reset_index(names="program"), tables_dir / "reference_program_membership.tsv")
        if not program_scores.empty:
            _write_frame(program_scores.reset_index(names="cell"), tables_dir / "program_scores_cell_level.tsv.gz")
        if not neighborhood_scores.empty:
            _write_frame(neighborhood_scores.reset_index(names="cell"), tables_dir / "neighbor_program_scores_cell_level.tsv.gz")
        if not reference_de.empty:
            _write_table(reference_de, tables_dir / "reference_de.tsv")

        references_dir.mkdir(parents=True, exist_ok=True)
        for reference_name, reference in reference_objects.items():
            reference.write_h5ad(references_dir / f"{reference_name}.h5ad")

        figure_paths: dict[str, str] = {}
        heatmap_path = _write_program_heatmap(grouped_scores, figures_dir / "program_scores_heatmap.png", title="Program scores by group")
        if heatmap_path is not None:
            figure_paths["program_scores_heatmap"] = heatmap_path
        neighbor_heatmap_path = _write_program_heatmap(
            grouped_neighbor_scores,
            figures_dir / "neighbor_program_scores_heatmap.png",
            title="Neighborhood program scores by group",
        )
        if neighbor_heatmap_path is not None:
            figure_paths["neighbor_program_scores_heatmap"] = neighbor_heatmap_path

        manifest = {
            "benchmark": "breast_reference_projection",
            "dataset": dataset_name,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "package_version": _package_version(),
            "config": cfg,
            "reference_datasets": list(map(str, reference_datasets)),
            "tables": {
                key: str(path)
                for key, path in {
                    "reference_summary": tables_dir / "reference_summary.tsv",
                    "program_scores_by_group": tables_dir / "program_scores_by_group.tsv",
                    "neighbor_program_scores_by_group": tables_dir / "neighbor_program_scores_by_group.tsv",
                    "reference_program_membership": tables_dir / "reference_program_membership.tsv",
                    "program_scores_cell_level": tables_dir / "program_scores_cell_level.tsv.gz",
                    "neighbor_program_scores_cell_level": tables_dir / "neighbor_program_scores_cell_level.tsv.gz",
                    "reference_de": tables_dir / "reference_de.tsv",
                }.items()
                if path.exists()
            },
            "figures": figure_paths,
            "summary": {
                "n_obs": int(spatial_adata.n_obs),
                "n_vars": int(spatial_adata.n_vars),
                "reference_count": int(len(reference_rows)),
                "program_count": int(program_scores.shape[1]),
            },
        }
        _write_manifest(output_root / "manifest.json", manifest)
        results["manifest"] = manifest
        results["report_dir"] = str(output_root)

    return results


def _load_spatial_for_publication(spatial_input: str | Path | AnnData, cfg: Mapping[str, Any]) -> tuple[AnnData, str]:
    if isinstance(spatial_input, AnnData):
        adata = ensure_spatialperturb_schema(spatial_input.copy())
        return adata, str(adata.uns.get("spatialperturb", {}).get("dataset_name", "custom_xenium"))
    spatial_path = Path(spatial_input).expanduser()
    if spatial_path.suffix == ".h5ad":
        adata = ad.read_h5ad(spatial_path)
        ensure_spatialperturb_schema(adata, metadata={"platform": "xenium", "source_path": str(spatial_path)})
        return adata, spatial_path.stem
    adata = read_xenium(
        spatial_path,
        cell_group_path=cfg.get("cell_group_path"),
        roi_geojson_path=cfg.get("roi_geojson_path"),
        sample_name=cfg.get("sample_name"),
        load_molecules=bool(cfg.get("load_molecules", False)),
    )
    return adata, spatial_path.stem


def _filter_reference_for_projection(name: str, reference: AnnData, cfg: Mapping[str, Any]) -> AnnData:
    if name != "gse281048_pathway_atlas":
        return reference
    pathway_cell_line = str(cfg.get("pathway_cell_line", "MCF7"))
    if "cell_line" not in reference.obs.columns:
        raise KeyError("gse281048_pathway_atlas requires obs['cell_line'] for MCF7 filtering.")
    filtered = reference[reference.obs["cell_line"].astype(str) == pathway_cell_line].copy()
    if filtered.n_obs == 0:
        raise ValueError(f"No cells remain after filtering gse281048_pathway_atlas to {pathway_cell_line!r}.")
    return filtered


def _load_references_for_publication(
    reference_datasets: Sequence[str],
    *,
    cache_dir: Path,
    cfg: Mapping[str, Any],
) -> tuple[dict[str, AnnData], dict[str, dict[str, Any]]]:
    configured = dict(cfg.get("reference_adatas") or {})
    references: dict[str, AnnData] = {}
    status: dict[str, dict[str, Any]] = {}
    strict = bool(cfg.get("strict_references", False))
    for name in map(str, reference_datasets):
        try:
            reference = ensure_spatialperturb_schema(configured[name].copy()) if name in configured else load_public_dataset(name, cache_dir=cache_dir)
            reference = _filter_reference_for_projection(name, reference, cfg)
            references[name] = reference
            status[name] = {
                "status": "ready",
                "n_obs": int(reference.n_obs),
                "n_vars": int(reference.n_vars),
            }
        except Exception as exc:
            if strict or name == "gse241115_breast_cropseq":
                raise
            status[name] = {
                "status": "blocked",
                "reason": f"{name.upper()}_UNAVAILABLE",
                "message": str(exc),
            }
    if not references:
        raise RuntimeError("No reference datasets were available for the Nature Methods breast analysis.")
    return references, status


def _programs_from_membership(membership: pd.DataFrame) -> dict[str, list[str]]:
    if membership.empty:
        return {}
    frame = membership.copy()
    if "program" in frame.columns:
        frame = frame.set_index("program")
    return {
        str(program): [str(gene) for gene, value in row.items() if int(value) == 1]
        for program, row in frame.iterrows()
    }


def _score_table_spearman(left: pd.DataFrame, right: pd.DataFrame) -> float:
    if left.empty or right.empty:
        return float("nan")
    merged = left.loc[:, ["group", "program", "mean_score"]].merge(
        right.loc[:, ["group", "program", "mean_score"]],
        on=["group", "program"],
        suffixes=("_left", "_right"),
    )
    if len(merged) < 2:
        return float("nan")
    return float(spearmanr(merged["mean_score_left"], merged["mean_score_right"]).statistic)


def _run_publication_ablations(
    spatial_adata: AnnData,
    references: Mapping[str, AnnData],
    baseline_grouped: pd.DataFrame,
    baseline_neighbor_grouped: pd.DataFrame,
    *,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    control = str(cfg.get("reference_control", "control"))
    groupby = cfg.get("groupby", ["cell_type", "roi"])
    for top_n in cfg.get("top_n_values", (25, 50, 100)):
        programs: dict[str, list[str]] = {}
        for name, reference in references.items():
            ref_programs = build_reference_programs(
                reference,
                control=control,
                method=str(cfg.get("reference_method", "auto")),
                sample_col=cfg.get("reference_sample_col", "sample" if "sample" in reference.obs.columns else None),
                top_n=int(top_n),
                direction=str(cfg.get("direction", "both")),
                effect_size_only=bool(cfg.get("reference_effect_size_only", True)),
            )
            programs.update({f"{name}:{program}": genes for program, genes in ref_programs.items()})
        scores = score_programs(spatial_adata, programs) if programs else pd.DataFrame(index=spatial_adata.obs_names.astype(str))
        grouped = aggregate_program_scores(spatial_adata, scores, groupby=groupby) if not scores.empty else pd.DataFrame()
        records.append(
            {
                "ablation_type": "top_n",
                "value": int(top_n),
                "program_count": int(scores.shape[1]) if not scores.empty else 0,
                "score_spearman_vs_primary": _score_table_spearman(baseline_grouped, grouped),
                "max_group_mean_score": float(grouped["mean_score"].max()) if not grouped.empty else np.nan,
            }
        )

    base_scores = _resolve_publication_scores(spatial_adata)
    for k in cfg.get("graph_k_values", (5, 15, 30)):
        copied = spatial_adata.copy()
        graph_key = f"sp_knn_k{int(k)}"
        build_spatial_graph(copied, mode="knn", k=int(k), graph_key=graph_key)
        neighbor = neighbor_program_scores(copied, score_key=base_scores, graph_key=graph_key, key_added=f"neighbor_program_scores_k{int(k)}")
        grouped_neighbor = aggregate_program_scores(copied, neighbor, groupby=groupby)
        records.append(
            {
                "ablation_type": "graph_k",
                "value": int(k),
                "program_count": int(base_scores.shape[1]),
                "score_spearman_vs_primary": _score_table_spearman(baseline_neighbor_grouped, grouped_neighbor),
                "max_group_mean_score": float(grouped_neighbor["mean_score"].max()) if not grouped_neighbor.empty else np.nan,
            }
        )
    return pd.DataFrame.from_records(records)


def _resolve_publication_scores(adata: AnnData) -> pd.DataFrame:
    raw = adata.obsm.get("program_scores")
    if isinstance(raw, pd.DataFrame):
        return raw.copy()
    if raw is None:
        return pd.DataFrame(index=adata.obs_names.astype(str))
    return pd.DataFrame(raw, index=adata.obs_names.astype(str))


def _write_nature_methods_summary(
    path: Path,
    *,
    manifest: Mapping[str, Any],
    top_claims: pd.DataFrame,
    reference_status: Mapping[str, Any],
) -> None:
    lines = [
        "# SpatialPerturb Nature Methods Brief Communication Summary",
        "",
        "## Format Target",
        "- Article type: Nature Methods Brief Communication.",
        "- Constraints: 70-word abstract, 1,200-word main text including abstract/references/legends, maximum two main display items, Online Methods with subheadings.",
        "",
        "## Run Summary",
        f"- Spatial cells: `{manifest.get('summary', {}).get('n_obs', 'NA')}`.",
        f"- Spatial genes/features: `{manifest.get('summary', {}).get('n_vars', 'NA')}`.",
        f"- Reference datasets: `{', '.join(map(str, manifest.get('reference_datasets', [])))}`.",
        f"- Program count: `{manifest.get('summary', {}).get('program_count', 'NA')}`.",
        "",
        "## Claim-Level Biological Readout",
    ]
    if top_claims.empty:
        lines.append("- No claim-level calibrated groups passed the minimum-cell filter.")
    else:
        for row in top_claims.head(10).itertuples(index=False):
            lines.append(
                f"- `{row.group}` resembles `{row.program}` "
                f"(mean={float(row.mean_score):.4g}, z={float(row.z_score):.3g}, FDR={float(row.fdr):.3g}, n={int(row.n_cells)})."
            )
    if reference_status:
        lines.extend(["", "## Reference Status"])
        for name, payload in reference_status.items():
            if isinstance(payload, Mapping):
                lines.append(f"- `{name}`: `{payload.get('status', 'unknown')}` {payload.get('reason', '')} {payload.get('message', '')}".strip())
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "- Projection scores indicate transcriptional similarity to a Perturb-seq-derived program, not a real perturbation in the tissue.",
            "- Cell-line references and FFPE Xenium tissue differ in context; claims should be framed as candidate regulatory-state hypotheses.",
            "- Small ROI/cell-type groups are retained in supplementary tables but excluded from claim-level summaries.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_publication_biological_interpretation(
    path: Path,
    *,
    calibrated: pd.DataFrame,
    spatial_stats: pd.DataFrame,
    reference_status: Mapping[str, Any],
) -> None:
    claims = calibrated.loc[calibrated.get("is_claim_level", False)].copy() if not calibrated.empty else pd.DataFrame()
    if not claims.empty:
        claims = claims.sort_values(["z_score", "mean_score"], ascending=[False, False])
    lines = [
        "# Biological Interpretation",
        "",
        "The strongest calibrated signals should be interpreted as spatial localization of Perturb-seq reference-like transcriptional states. They do not show that the Xenium tissue contains the corresponding knockout, CRISPRi perturbation, or drug treatment.",
        "",
        "## Main Finding",
    ]
    if claims.empty:
        lines.append("- No calibrated claim-level group was available after applying the minimum-cell filter.")
    else:
        for row in claims.head(12).itertuples(index=False):
            program = str(row.program)
            if any(token in program.upper() for token in ("PRPF6", "UXT", "SS18L2", "SMAD5", "PMF1", "PFDN5")):
                meaning = "luminal/secretory invasive tumor-state similarity, driven by epithelial genes such as TFF1/TFF3/AGR2/IGFBP5 when present in the program."
            elif any(token in program.upper() for token in ("IFN", "STAT", "IRF", "TNF", "TGFB")):
                meaning = "cytokine or pathway-response similarity that may localize inflammatory or stromal microenvironments."
            else:
                meaning = "reference perturbation-associated transcriptional-state similarity."
            lines.append(f"- `{row.group}`: `{program}` (z={float(row.z_score):.3g}, FDR={float(row.fdr):.3g}, n={int(row.n_cells)}), interpreted as {meaning}")
    if not spatial_stats.empty:
        lines.extend(["", "## Spatial Organization"])
        for row in spatial_stats.sort_values("moran_i", ascending=False).head(8).itertuples(index=False):
            lines.append(f"- `{row.program}` has Moran-style spatial autocorrelation `{float(row.moran_i):.3g}` (FDR={float(row.fdr):.3g}), supporting local organization of the state.")
    blocked = {name: payload for name, payload in reference_status.items() if isinstance(payload, Mapping) and payload.get("status") == "blocked"}
    if blocked:
        lines.extend(["", "## Blocked Optional References"])
        for name, payload in blocked.items():
            lines.append(f"- `{name}`: `{payload.get('reason', 'blocked')}`; {payload.get('message', '')}")
    lines.extend(
        [
            "",
            "## Caveats",
            "- GSE241115 and GSE281048 are cell-line references; the query is FFPE breast tissue.",
            "- Projection is an association-style score and should motivate validation rather than serve as causal proof.",
            "- ROI and cell-group annotation quality directly shapes all biological conclusions.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _render_nature_methods_figures(
    *,
    adata: AnnData,
    calibrated: pd.DataFrame,
    validation: pd.DataFrame,
    spatial_stats: pd.DataFrame,
    output_dir: Path,
    min_cells: int,
    seed: int,
) -> dict[str, str]:
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: dict[str, str] = {}
    claims = calibrated.loc[calibrated["n_cells"] >= int(min_cells)].copy() if not calibrated.empty else pd.DataFrame()
    if not claims.empty:
        claims = claims.sort_values(["z_score", "mean_score"], ascending=[False, False])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].axis("off")
    axes[0, 0].text(
        0.02,
        0.96,
        "SpatialPerturb\nPerturb-seq reference programs -> Xenium WTA tissue\nnull calibration + held-out validation + spatial graph statistics",
        va="top",
        ha="left",
        fontsize=12,
        weight="bold",
    )
    if not validation.empty and "auroc" in validation.columns:
        plot = validation.dropna(subset=["auroc"]).head(12)
        sns.barplot(data=plot, x="auroc", y="program", ax=axes[0, 1], color="#4C78A8")
        axes[0, 1].set_xlim(0, 1)
        axes[0, 1].set_title("Held-out reference recovery")
    else:
        axes[0, 1].axis("off")
        axes[0, 1].text(0.1, 0.5, "Reference validation unavailable")
    if not claims.empty:
        axes[1, 0].hist(claims["z_score"].astype(float), bins=30, color="#F58518", alpha=0.85)
        axes[1, 0].set_title("Null-calibrated claim-level z scores")
        axes[1, 0].set_xlabel("Empirical z")
    else:
        axes[1, 0].axis("off")
        axes[1, 0].text(0.1, 0.5, "No claim-level calibrated scores")
    if not spatial_stats.empty:
        plot = spatial_stats.sort_values("moran_i", ascending=False).head(12)
        sns.barplot(data=plot, x="moran_i", y="program", ax=axes[1, 1], color="#54A24B")
        axes[1, 1].set_title("Spatial autocorrelation")
    else:
        axes[1, 1].axis("off")
        axes[1, 1].text(0.1, 0.5, "Spatial autocorrelation unavailable")
    fig.tight_layout()
    path = figures_dir / "main_figure_1.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    figure_paths["main_figure_1"] = str(path)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    scores = _resolve_publication_scores(adata)
    if not claims.empty and "spatial" in adata.obsm and not scores.empty:
        top_program = str(claims.iloc[0]["program"])
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        rng = np.random.default_rng(seed)
        idx = np.arange(adata.n_obs)
        if idx.size > 20000:
            idx = rng.choice(idx, size=20000, replace=False)
        color = scores[top_program].to_numpy(dtype=float)[idx] if top_program in scores.columns else np.zeros(len(idx))
        scatter = axes[0, 0].scatter(coords[idx, 0], coords[idx, 1], c=color, s=1, cmap="magma", linewidths=0)
        axes[0, 0].set_title(f"Spatial map: {top_program}")
        axes[0, 0].set_aspect("equal", adjustable="box")
        fig.colorbar(scatter, ax=axes[0, 0], fraction=0.046, pad=0.04)
    else:
        axes[0, 0].axis("off")
        axes[0, 0].text(0.1, 0.5, "Spatial program map unavailable")

    if not claims.empty:
        top_programs = claims["program"].drop_duplicates().head(8).tolist()
        top_groups = claims["group"].drop_duplicates().head(10).tolist()
        heatmap = claims[claims["program"].isin(top_programs) & claims["group"].isin(top_groups)].pivot_table(
            index="group",
            columns="program",
            values="z_score",
            aggfunc="mean",
        )
        sns.heatmap(heatmap.fillna(0), cmap="vlag", center=0, ax=axes[0, 1])
        axes[0, 1].set_title("Claim-level calibrated enrichment")
    else:
        axes[0, 1].axis("off")
        axes[0, 1].text(0.1, 0.5, "No claim-level heatmap")

    if not spatial_stats.empty:
        plot = spatial_stats.sort_values("moran_i", ascending=False).head(10)
        sns.scatterplot(data=plot, x="moran_i", y="z_score", size="n_permutations", legend=False, ax=axes[1, 0], color="#B279A2")
        for row in plot.itertuples(index=False):
            axes[1, 0].text(float(row.moran_i), float(row.z_score), str(row.program).split(":", 1)[-1], fontsize=7)
        axes[1, 0].set_title("Spatially organized programs")
    else:
        axes[1, 0].axis("off")

    if not claims.empty:
        plot = claims.head(12).copy()
        plot["label"] = plot["program"].astype(str).str.split(":", n=1).str[-1]
        sns.barplot(data=plot, x="z_score", y="label", ax=axes[1, 1], color="#E45756")
        axes[1, 1].set_title("Top robust biological programs")
    else:
        axes[1, 1].axis("off")
    fig.tight_layout()
    path = figures_dir / "main_figure_2.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    figure_paths["main_figure_2"] = str(path)
    return figure_paths


def run_nature_methods_breast_analysis(
    spatial_input: str | Path | AnnData,
    *,
    reference_datasets: Sequence[str] | None = None,
    config: Mapping[str, Any] | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Run the publication-grade breast Xenium analysis for a Nature Methods Brief Communication."""
    cfg = {
        "cache_dir": ".spatialperturb-cache",
        "groupby": ["cell_type", "roi"],
        "k": 15,
        "top_n": 50,
        "top_n_values": (25, 50, 100),
        "graph_k_values": (5, 15, 30),
        "reference_control": "control",
        "reference_method": "auto",
        "reference_effect_size_only": True,
        "direction": "both",
        "pathway_cell_line": "MCF7",
        "n_random": 25,
        "n_label_shuffles": 25,
        "n_spatial_permutations": 25,
        "n_bootstrap": 100,
        "min_claim_cells": 50,
        "seed": 0,
    }
    cfg.update(dict(config or {}))
    cache_dir = Path(cfg.get("cache_dir", ".spatialperturb-cache")).expanduser().resolve()
    output_root = Path(output_dir or "reports/nature_methods_breast_shortcomm").expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    tables_dir = output_root / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    requested_references = list(reference_datasets or ["gse241115_breast_cropseq", "gse281048_pathway_atlas"])
    spatial_adata, dataset_name = _load_spatial_for_publication(spatial_input, cfg)
    if "sp_knn" not in spatial_adata.obsp and "sp_radius" not in spatial_adata.obsp:
        build_spatial_graph(spatial_adata, mode="knn", k=int(cfg.get("k", 15)))

    references, reference_status = _load_references_for_publication(requested_references, cache_dir=cache_dir, cfg=cfg)
    projection_cfg = dict(cfg)
    projection_cfg["reference_adatas"] = references
    projection_cfg["top_n"] = int(cfg.get("top_n", 50))
    projection = run_reference_projection_benchmark(
        spatial_adata,
        reference_datasets=list(references.keys()),
        config=projection_cfg,
        output_dir=output_root,
    )
    spatial_adata = projection["adata"]
    programs = _programs_from_membership(projection["reference_program_membership"])

    calibrated, calibration_nulls = calibrate_program_scores(
        spatial_adata,
        programs,
        score_key="program_scores",
        groupby=cfg.get("groupby", ["cell_type", "roi"]),
        n_random=int(cfg.get("n_random", 25)),
        n_label_shuffles=int(cfg.get("n_label_shuffles", 25)),
        min_cells=int(cfg.get("min_claim_cells", 50)),
        seed=int(cfg.get("seed", 0)),
        return_nulls=True,
    )
    bootstrap = bootstrap_program_score_intervals(
        spatial_adata,
        "program_scores",
        groupby=cfg.get("groupby", ["cell_type", "roi"]),
        n_bootstrap=int(cfg.get("n_bootstrap", 100)),
        min_cells=int(cfg.get("min_claim_cells", 50)),
        seed=int(cfg.get("seed", 0)) + 11,
    )
    if not bootstrap.empty:
        calibrated = calibrated.merge(
            bootstrap.loc[:, ["group", "program", "ci_low", "ci_high", "n_bootstrap"]],
            on=["group", "program"],
            how="left",
        )
    spatial_stats = spatial_autocorrelation_scores(
        spatial_adata,
        score_key="program_scores",
        n_permutations=int(cfg.get("n_spatial_permutations", 25)),
        seed=int(cfg.get("seed", 0)) + 23,
    )
    redundancy = program_redundancy_table(programs)

    validation_tables: list[pd.DataFrame] = []
    for name, reference in references.items():
        try:
            validation = validate_reference_programs(
                reference,
                query_adata=spatial_adata,
                control=str(cfg.get("reference_control", "control")),
                top_n=int(cfg.get("top_n", 50)),
                direction=str(cfg.get("direction", "both")),
                method="simple",
                seed=int(cfg.get("seed", 0)) + 37,
            )
            validation["reference_dataset"] = name
            validation_tables.append(validation)
        except Exception as exc:
            reference_status.setdefault(name, {})["validation_status"] = "blocked"
            reference_status.setdefault(name, {})["validation_message"] = str(exc)
    reference_validation = _concat_or_empty(validation_tables)

    ablation = _run_publication_ablations(
        spatial_adata,
        references,
        projection["program_scores_by_group"],
        projection["neighbor_program_scores_by_group"],
        cfg=cfg,
    )

    _write_table(calibrated, tables_dir / "calibrated_program_scores_by_group.tsv")
    if not calibration_nulls.empty:
        _write_frame(calibration_nulls, tables_dir / "calibration_nulls.tsv.gz")
    if not reference_validation.empty:
        _write_table(reference_validation, tables_dir / "reference_validation.tsv")
    _write_table(spatial_stats, tables_dir / "spatial_autocorrelation.tsv")
    _write_table(ablation, tables_dir / "ablation_summary.tsv")
    if not redundancy.empty:
        _write_table(redundancy, tables_dir / "program_redundancy.tsv")
    if not bootstrap.empty:
        _write_table(bootstrap, tables_dir / "bootstrap_program_score_intervals.tsv")

    figure_paths = _render_nature_methods_figures(
        adata=spatial_adata,
        calibrated=calibrated,
        validation=reference_validation,
        spatial_stats=spatial_stats,
        output_dir=output_root,
        min_cells=int(cfg.get("min_claim_cells", 50)),
        seed=int(cfg.get("seed", 0)),
    )
    top_claims = calibrated.loc[calibrated["n_cells"] >= int(cfg.get("min_claim_cells", 50))].sort_values(
        ["z_score", "mean_score"], ascending=[False, False]
    )

    manifest = {
        "benchmark": "nature_methods_breast_shortcomm",
        "dataset": dataset_name,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "package_version": _package_version(),
        "config": cfg,
        "reference_datasets": list(references.keys()),
        "requested_reference_datasets": requested_references,
        "reference_status": reference_status,
        "tables": {
            key: str(path)
            for key, path in {
                "calibrated_program_scores_by_group": tables_dir / "calibrated_program_scores_by_group.tsv",
                "calibration_nulls": tables_dir / "calibration_nulls.tsv.gz",
                "reference_validation": tables_dir / "reference_validation.tsv",
                "spatial_autocorrelation": tables_dir / "spatial_autocorrelation.tsv",
                "ablation_summary": tables_dir / "ablation_summary.tsv",
                "program_redundancy": tables_dir / "program_redundancy.tsv",
                "bootstrap_program_score_intervals": tables_dir / "bootstrap_program_score_intervals.tsv",
            }.items()
            if path.exists()
        },
        "figures": figure_paths,
        "summary": {
            "n_obs": int(spatial_adata.n_obs),
            "n_vars": int(spatial_adata.n_vars),
            "reference_count": int(len(references)),
            "program_count": int(len(programs)),
            "claim_level_rows": int(len(top_claims)),
        },
    }
    _write_manifest(output_root / "manifest.json", manifest)
    _write_nature_methods_summary(output_root / "nature_methods_summary.md", manifest=manifest, top_claims=top_claims, reference_status=reference_status)
    _write_publication_biological_interpretation(
        output_root / "biological_interpretation.md",
        calibrated=calibrated,
        spatial_stats=spatial_stats,
        reference_status=reference_status,
    )

    projection.update(
        {
            "manifest": manifest,
            "report_dir": str(output_root),
            "calibrated_program_scores_by_group": calibrated,
            "calibration_nulls": calibration_nulls,
            "reference_validation": reference_validation,
            "spatial_autocorrelation": spatial_stats,
            "ablation_summary": ablation,
            "program_redundancy": redundancy,
        }
    )
    return projection
