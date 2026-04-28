"""Run the publication-grade breast Xenium analysis with A100 defaults."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import spatialperturb as sp


DATA_ROOT = Path("/data/taobo.hu/SpatialPerturb")
CACHE_DIR = DATA_ROOT / "cache"
INPUT_DIR = DATA_ROOT / "inputs" / "xenium_wta_breast"
PREPARED_DIR = DATA_ROOT / "prepared"
REPORT_DIR = DATA_ROOT / "reports" / "nature_methods_breast_shortcomm"
SPATIAL_OUTPUT = PREPARED_DIR / "xenium_wta_breast.h5ad"
REFERENCE_STATUS_PATH = REPORT_DIR / "reference_prepare_status.json"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def _r_seurat_available() -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False, "Rscript is not available on PATH."
    result = subprocess.run(
        [rscript, "-e", "suppressPackageStartupMessages({library(Seurat); library(Matrix)}); cat('Seurat OK\\n')"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return False, result.stderr.strip() or result.stdout.strip() or "Seurat could not be loaded."
    return True, result.stdout.strip()


def _prepare_xenium() -> Path:
    if SPATIAL_OUTPUT.exists():
        return SPATIAL_OUTPUT
    PREPARED_DIR.mkdir(parents=True, exist_ok=True)
    cell_group_path = INPUT_DIR / "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv"
    roi_geojson_path = INPUT_DIR / "xenium_explorer_annotations.geojson"
    required_inputs = [INPUT_DIR / "cell_feature_matrix.h5", INPUT_DIR / "cells.csv.gz"]
    missing_inputs = [path for path in required_inputs if not path.exists()]
    if missing_inputs:
        raise FileNotFoundError(
            "Missing required Xenium inputs. Run scripts/a100_sync_xenium_minimal.ps1 first. "
            f"Missing: {', '.join(map(str, missing_inputs))}"
        )
    adata = sp.read_xenium(
        INPUT_DIR,
        cell_group_path=cell_group_path if cell_group_path.exists() else None,
        roi_geojson_path=roi_geojson_path if roi_geojson_path.exists() else None,
        sample_name="xenium_wta_breast",
    )
    adata.write_h5ad(SPATIAL_OUTPUT)
    return SPATIAL_OUTPUT


def _prepare_reference_dataset(name: str) -> dict[str, Any]:
    fetch = sp.fetch_dataset(name, cache_dir=CACHE_DIR)
    prepare = sp.prepare_dataset(name, cache_dir=CACHE_DIR)
    return {"status": "ready", "fetch": fetch, "prepare": prepare}


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    spatial_path = _prepare_xenium()
    reference_status: dict[str, Any] = {}

    print("[1/5] Preparing GSE241115 breast CROP-seq reference...", flush=True)
    reference_status["gse241115_breast_cropseq"] = _prepare_reference_dataset("gse241115_breast_cropseq")

    print("[2/5] Preparing GSE281048 pathway atlas if R/Seurat is available...", flush=True)
    r_ready, r_message = _r_seurat_available()
    if r_ready:
        try:
            reference_status["gse281048_pathway_atlas"] = {
                "r_check": r_message,
                **_prepare_reference_dataset("gse281048_pathway_atlas"),
            }
        except Exception as exc:
            reference_status["gse281048_pathway_atlas"] = {
                "status": "blocked",
                "reason": "GSE281048_PREPARE_FAILED",
                "message": str(exc),
                "r_check": r_message,
            }
    else:
        reference_status["gse281048_pathway_atlas"] = {
            "status": "blocked",
            "reason": "GSE281048_BLOCKED_R_ENV_MISSING",
            "message": r_message,
        }
    _write_json(REFERENCE_STATUS_PATH, reference_status)

    print("[3/5] Running Nature Methods breast short-communication analysis...", flush=True)
    results = sp.run_nature_methods_breast_analysis(
        spatial_path,
        reference_datasets=["gse241115_breast_cropseq", "gse281048_pathway_atlas"],
        config={
            "cache_dir": CACHE_DIR,
            "cell_group_path": INPUT_DIR / "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv",
            "roi_geojson_path": INPUT_DIR / "xenium_explorer_annotations.geojson",
            "sample_name": "xenium_wta_breast",
            "reference_effect_size_only": True,
            "pathway_cell_line": "MCF7",
            "n_random": 25,
            "n_label_shuffles": 25,
            "n_spatial_permutations": 25,
            "n_bootstrap": 100,
            "min_claim_cells": 50,
            "seed": 20260428,
        },
        output_dir=REPORT_DIR,
    )

    print("[4/5] Packaging manuscript-facing summary...", flush=True)
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "package_nature_methods_shortcomm_report.py"),
            "--report-dir",
            str(REPORT_DIR),
        ],
        check=True,
    )

    print("[5/5] Complete.", flush=True)
    print(
        json.dumps(
            {
                "report_dir": str(REPORT_DIR),
                "manifest": str(REPORT_DIR / "manifest.json"),
                "nature_methods_summary": str(REPORT_DIR / "nature_methods_summary.md"),
                "biological_interpretation": str(REPORT_DIR / "biological_interpretation.md"),
                "reference_status": reference_status,
                "summary": results.get("manifest", {}).get("summary", {}),
            },
            indent=2,
            default=str,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
