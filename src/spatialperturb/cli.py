"""Command-line interface for SpatialPerturb."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import typer

from . import __version__
from .benchmarks import available_benchmarks, run_core_benchmark, run_nature_methods_breast_analysis, run_reference_projection_benchmark
from .datasets import available_datasets, fetch_dataset, prepare_dataset
from .io import read_xenium
from .schema import schema_summary, validate_spatialperturb_schema

app = typer.Typer(help="SpatialPerturb CLI", add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Command group for SpatialPerturb utilities."""


@app.command()
def version() -> None:
    """Print the installed package version."""
    typer.echo(__version__)


@app.command("datasets")
def list_datasets() -> None:
    """List built-in dataset cards."""
    typer.echo(available_datasets().to_string(index=False))


@app.command("benchmarks")
def list_benchmarks() -> None:
    """List benchmark tracks."""
    typer.echo(available_benchmarks().to_string(index=False))


@app.command("fetch-dataset")
def fetch_dataset_command(
    name: str,
    cache_dir: Path = Path(".spatialperturb-cache"),
    force: bool = False,
) -> None:
    """Download raw public files for a registered dataset."""
    result = fetch_dataset(name, cache_dir=cache_dir, force=force)
    typer.echo(f"dataset={result['dataset']}")
    typer.echo(f"status={result['status']}")
    typer.echo(f"raw_dir={result['raw_dir']}")


@app.command("prepare-dataset")
def prepare_dataset_command(
    name: str,
    cache_dir: Path = Path(".spatialperturb-cache"),
    output_dir: Path | None = None,
) -> None:
    """Prepare a dataset into SpatialPerturb-compatible outputs."""
    result = prepare_dataset(name, cache_dir=cache_dir, output_dir=output_dir)
    typer.echo(f"dataset={result['dataset']}")
    typer.echo(f"status={result['status']}")
    typer.echo(f"prepared_path={result.get('prepared_path')}")
    if result.get("note"):
        typer.echo(result["note"])


@app.command("prepare-xenium")
def prepare_xenium_command(
    path: Path,
    output_path: Path,
    cell_group_path: Path | None = None,
    roi_geojson_path: Path | None = None,
    sample_name: str | None = None,
    load_molecules: bool = False,
) -> None:
    """Prepare a Xenium directory into a schema-compliant h5ad file."""
    adata = read_xenium(
        path,
        cell_group_path=cell_group_path,
        roi_geojson_path=roi_geojson_path,
        sample_name=sample_name,
        load_molecules=load_molecules,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    typer.echo(f"prepared_path={output_path.resolve()}")
    typer.echo(f"n_obs={adata.n_obs}")
    typer.echo(f"n_vars={adata.n_vars}")


@app.command("run-benchmark")
def run_benchmark_command(
    dataset: str,
    cache_dir: Path = Path(".spatialperturb-cache"),
    output_dir: Path | None = None,
    reference_dataset: str | None = None,
    method: str = "auto",
) -> None:
    """Run the paper-grade benchmark workflow for a dataset."""
    config: dict[str, object] = {"cache_dir": cache_dir}
    if reference_dataset is not None:
        config["reference_dataset"] = reference_dataset
    if method != "auto":
        config["method"] = method
    report_root = output_dir or Path("reports") / dataset
    results = run_core_benchmark(dataset, config=config, output_dir=report_root)
    typer.echo(f"dataset={dataset}")
    typer.echo(f"report_dir={report_root.resolve()}")
    typer.echo(f"tables={sorted(key for key, value in results.items() if hasattr(value, 'to_csv'))}")


@app.command("render-paper-figures")
def render_paper_figures_command(
    dataset: str,
    cache_dir: Path = Path(".spatialperturb-cache"),
    output_dir: Path | None = None,
    reference_dataset: str | None = None,
) -> None:
    """Render the canonical paper figures for a dataset."""
    config: dict[str, object] = {"cache_dir": cache_dir}
    if reference_dataset is not None:
        config["reference_dataset"] = reference_dataset
    report_root = output_dir or Path("reports") / dataset
    run_core_benchmark(dataset, config=config, output_dir=report_root)
    typer.echo(f"dataset={dataset}")
    typer.echo(f"figures_dir={(report_root.resolve() / 'figures')}")


@app.command("run-reference-benchmark")
def run_reference_benchmark_command(
    spatial_input: Path,
    output_dir: Path,
    cache_dir: Path = Path(".spatialperturb-cache"),
    cell_group_path: Path | None = None,
    roi_geojson_path: Path | None = None,
    sample_name: str | None = None,
    reference_datasets: str = "gse241115_breast_cropseq,gse281048_pathway_atlas",
) -> None:
    """Run the Xenium + Perturb-seq reference projection workflow."""
    datasets = [name.strip() for name in reference_datasets.split(",") if name.strip()]
    config: dict[str, object] = {
        "cache_dir": cache_dir,
        "cell_group_path": cell_group_path,
        "roi_geojson_path": roi_geojson_path,
        "sample_name": sample_name,
    }
    results = run_reference_projection_benchmark(
        spatial_input,
        reference_datasets=datasets,
        config=config,
        output_dir=output_dir,
    )
    typer.echo(f"dataset={results['manifest']['dataset'] if 'manifest' in results else spatial_input}")
    typer.echo(f"report_dir={output_dir.resolve()}")
    typer.echo(f"references={','.join(datasets)}")


@app.command("run-nature-methods-breast-analysis")
def run_nature_methods_breast_analysis_command(
    spatial_input: Path,
    output_dir: Path,
    cache_dir: Path = Path(".spatialperturb-cache"),
    cell_group_path: Path | None = None,
    roi_geojson_path: Path | None = None,
    sample_name: str | None = None,
    reference_datasets: str = "gse241115_breast_cropseq,gse281048_pathway_atlas",
    n_random: int = 25,
    n_spatial_permutations: int = 25,
    n_bootstrap: int = 100,
    min_claim_cells: int = 50,
) -> None:
    """Run the publication-grade Nature Methods breast short-communication workflow."""
    datasets = [name.strip() for name in reference_datasets.split(",") if name.strip()]
    config: dict[str, object] = {
        "cache_dir": cache_dir,
        "cell_group_path": cell_group_path,
        "roi_geojson_path": roi_geojson_path,
        "sample_name": sample_name,
        "n_random": n_random,
        "n_label_shuffles": n_random,
        "n_spatial_permutations": n_spatial_permutations,
        "n_bootstrap": n_bootstrap,
        "min_claim_cells": min_claim_cells,
        "reference_effect_size_only": True,
    }
    results = run_nature_methods_breast_analysis(
        spatial_input,
        reference_datasets=datasets,
        config=config,
        output_dir=output_dir,
    )
    manifest = results.get("manifest", {})
    typer.echo(f"dataset={manifest.get('dataset', spatial_input)}")
    typer.echo(f"report_dir={output_dir.resolve()}")
    typer.echo(f"references={','.join(map(str, manifest.get('reference_datasets', datasets)))}")
    typer.echo(f"claim_level_rows={manifest.get('summary', {}).get('claim_level_rows', 'NA')}")


@app.command("validate")
def validate(path: Path) -> None:
    """Validate that an h5ad file follows the SpatialPerturb schema."""
    adata = ad.read_h5ad(path)
    validate_spatialperturb_schema(adata)
    typer.echo(schema_summary(adata).rename_axis("metric").reset_index(name="value").to_string(index=False))


if __name__ == "__main__":
    app()
