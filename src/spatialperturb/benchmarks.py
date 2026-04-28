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
import pandas as pd
import seaborn as sns
from anndata import AnnData

from .datasets import available_datasets, load_public_dataset
from .gr import build_spatial_graph
from .io import read_xenium
from .reports import render_paper_figures
from .schema import ensure_spatialperturb_schema
from .signatures import (
    aggregate_program_scores,
    build_reference_programs,
    build_signature_matrix,
    neighbor_program_scores,
    score_programs,
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
