"""Run the Xenium breast reference projection workflow with A100-friendly defaults."""

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
REPORT_DIR = DATA_ROOT / "reports" / "breast_reference_projection"
SPATIAL_OUTPUT = PREPARED_DIR / "xenium_wta_breast.h5ad"
REFERENCE_STATUS_PATH = REPORT_DIR / "reference_status.json"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def _r_seurat_available() -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False, "Rscript is not available on PATH."
    result = subprocess.run(
        [rscript, "-e", "suppressPackageStartupMessages(library(Seurat)); cat('Seurat OK\\n')"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return False, result.stderr.strip() or result.stdout.strip() or "Seurat could not be loaded."
    return True, result.stdout.strip()


def _prepare_reference_dataset(name: str) -> dict[str, Any]:
    print(f"Preparing reference dataset: {name}", flush=True)
    fetch = sp.fetch_dataset(name, cache_dir=CACHE_DIR)
    prepare = sp.prepare_dataset(name, cache_dir=CACHE_DIR)
    return {"fetch": fetch, "prepare": prepare}


def _inject_reference_status(manifest_path: Path, reference_status: dict[str, Any]) -> None:
    if not manifest_path.exists():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["reference_status"] = reference_status
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
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

    print("[1/5] Preparing Xenium breast WTA dataset...", flush=True)
    adata = sp.read_xenium(
        INPUT_DIR,
        cell_group_path=cell_group_path if cell_group_path.exists() else None,
        roi_geojson_path=roi_geojson_path if roi_geojson_path.exists() else None,
        sample_name="xenium_wta_breast",
    )
    adata.write_h5ad(SPATIAL_OUTPUT)

    reference_status: dict[str, Any] = {}
    reference_datasets = ["gse241115_breast_cropseq"]

    print("[2/5] Fetching and preparing breast-specific CROP-seq reference...", flush=True)
    reference_status["gse241115_breast_cropseq"] = {"status": "ready", **_prepare_reference_dataset("gse241115_breast_cropseq")}

    print("[3/5] Checking optional pathway Perturb-seq atlas prerequisites...", flush=True)
    r_ready, r_message = _r_seurat_available()
    if r_ready:
        try:
            reference_status["gse281048_pathway_atlas"] = {
                "status": "ready",
                "r_check": r_message,
                **_prepare_reference_dataset("gse281048_pathway_atlas"),
            }
            reference_datasets.append("gse281048_pathway_atlas")
        except Exception as exc:
            reference_status["gse281048_pathway_atlas"] = {
                "status": "blocked",
                "reason": "GSE281048_PREPARE_FAILED",
                "message": str(exc),
            }
    else:
        reference_status["gse281048_pathway_atlas"] = {
            "status": "blocked",
            "reason": "GSE281048_BLOCKED_RSCRIPT_MISSING",
            "message": r_message,
        }
    _write_json(REFERENCE_STATUS_PATH, reference_status)

    print("[4/5] Running reference projection benchmark...", flush=True)
    results = sp.run_reference_projection_benchmark(
        SPATIAL_OUTPUT,
        reference_datasets=reference_datasets,
        config={
            "cache_dir": CACHE_DIR,
            "cell_group_path": cell_group_path if cell_group_path.exists() else None,
            "roi_geojson_path": roi_geojson_path if roi_geojson_path.exists() else None,
            "sample_name": "xenium_wta_breast",
            "pathway_cell_line": "MCF7",
            "reference_effect_size_only": True,
            "k": 15,
        },
        output_dir=REPORT_DIR,
    )

    manifest_path = REPORT_DIR / "manifest.json"
    _inject_reference_status(manifest_path, reference_status)

    print("[5/5] Writing biological interpretation...", flush=True)
    from interpret_breast_reference_projection import interpret_report

    interpretation = interpret_report(REPORT_DIR)
    manifest = results.get("manifest", {})
    print(
        json.dumps(
            {
                "report_dir": str(REPORT_DIR),
                "manifest": str(manifest_path),
                "interpretation": str(interpretation["markdown_path"]),
                "reference_datasets": reference_datasets,
                "summary": manifest.get("summary", {}),
                "reference_status": reference_status,
            },
            indent=2,
            default=str,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
