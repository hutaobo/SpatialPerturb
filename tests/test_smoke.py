from pathlib import Path

import pandas as pd
import spatialperturb as sp
from typer.testing import CliRunner

import spatialperturb.cli as cli_module
from spatialperturb.cli import app


def test_import():
    assert hasattr(sp, "__version__")
    assert hasattr(sp, "build_signature_matrix")
    assert hasattr(sp, "build_reference_programs")
    assert hasattr(sp, "intrinsic_de")
    assert hasattr(sp, "load_demo_dataset")


def test_cli_version():
    runner = CliRunner()
    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0
    assert sp.__version__ in result.stdout


def test_cli_lists_datasets_and_benchmarks():
    runner = CliRunner()
    datasets = runner.invoke(app, ["datasets"])
    benchmarks = runner.invoke(app, ["benchmarks"])

    assert datasets.exit_code == 0
    assert "shen_2026_stereoseq" in datasets.stdout
    assert benchmarks.exit_code == 0
    assert "cross_platform_concordance" in benchmarks.stdout


def test_dataset_lifecycle_and_cli_workflow(tmp_path):
    fetch = sp.fetch_dataset("demo_spatialperturb", cache_dir=tmp_path)
    prepare = sp.prepare_dataset("demo_spatialperturb", cache_dir=tmp_path)
    adata = sp.load_public_dataset("demo_spatialperturb", cache_dir=tmp_path)

    assert fetch["status"] == "built_in"
    assert prepare["status"] == "ready"
    assert adata.n_obs > 0

    runner = CliRunner()
    prepare_cli = runner.invoke(app, ["prepare-dataset", "demo_spatialperturb", "--cache-dir", str(tmp_path)])
    benchmark_cli = runner.invoke(
        app,
        [
            "run-benchmark",
            "demo_spatialperturb",
            "--cache-dir",
            str(tmp_path),
            "--output-dir",
            str(tmp_path / "reports"),
        ],
    )
    render_cli = runner.invoke(
        app,
        [
            "render-paper-figures",
            "demo_spatialperturb",
            "--cache-dir",
            str(tmp_path),
            "--output-dir",
            str(tmp_path / "figures-report"),
        ],
    )

    assert prepare_cli.exit_code == 0
    assert benchmark_cli.exit_code == 0
    assert render_cli.exit_code == 0
    assert (tmp_path / "reports" / "manifest.json").exists()
    assert (tmp_path / "figures-report" / "figures" / "workflow_schema.png").exists()


def test_prepare_xenium_and_run_reference_cli_commands(tmp_path, monkeypatch):
    expression = pd.DataFrame([[1, 0], [0, 2]], index=["cell_a", "cell_b"], columns=["GeneA", "GeneB"])
    cells = pd.DataFrame({"x": [0, 1], "y": [1, 0]}, index=expression.index)
    xenium_dir = tmp_path / "xenium"
    xenium_dir.mkdir()
    expression.to_csv(xenium_dir / "counts.csv")
    cells.to_csv(xenium_dir / "cells.csv")

    runner = CliRunner()
    prepare_result = runner.invoke(
        app,
        [
            "prepare-xenium",
            str(xenium_dir),
            str(tmp_path / "prepared_xenium.h5ad"),
        ],
    )

    def fake_run_reference_projection_benchmark(spatial_input, *, reference_datasets, config, output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        manifest_path = Path(output_dir) / "manifest.json"
        manifest_path.write_text('{"dataset":"demo"}', encoding="utf-8")
        return {"manifest": {"dataset": "demo"}}

    monkeypatch.setattr(cli_module, "run_reference_projection_benchmark", fake_run_reference_projection_benchmark)
    reference_result = runner.invoke(
        app,
        [
            "run-reference-benchmark",
            str(tmp_path / "prepared_xenium.h5ad"),
            str(tmp_path / "reference-report"),
        ],
    )

    assert prepare_result.exit_code == 0
    assert (tmp_path / "prepared_xenium.h5ad").exists()
    assert reference_result.exit_code == 0
    assert (tmp_path / "reference-report" / "manifest.json").exists()


def test_run_nature_methods_breast_analysis_cli_command(tmp_path, monkeypatch):
    runner = CliRunner()

    def fake_run_nature_methods_breast_analysis(spatial_input, *, reference_datasets, config, output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        (Path(output_dir) / "manifest.json").write_text('{"benchmark":"nature_methods_breast_shortcomm"}', encoding="utf-8")
        return {
            "manifest": {
                "dataset": "demo",
                "reference_datasets": reference_datasets,
                "summary": {"claim_level_rows": 3},
            }
        }

    monkeypatch.setattr(cli_module, "run_nature_methods_breast_analysis", fake_run_nature_methods_breast_analysis)
    result = runner.invoke(
        app,
        [
            "run-nature-methods-breast-analysis",
            str(tmp_path / "prepared_xenium.h5ad"),
            str(tmp_path / "nature-report"),
            "--n-random",
            "2",
            "--n-spatial-permutations",
            "2",
            "--n-bootstrap",
            "3",
        ],
    )

    assert result.exit_code == 0
    assert "claim_level_rows=3" in result.stdout
    assert (tmp_path / "nature-report" / "manifest.json").exists()
