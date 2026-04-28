"""Run the Xenium breast reference projection workflow with A100-friendly defaults."""

from __future__ import annotations

import json
from pathlib import Path
import sys

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


def main() -> None:
    input_path = INPUT_DIR
    cell_group_path = INPUT_DIR / "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv"
    roi_geojson_path = INPUT_DIR / "xenium_explorer_annotations.geojson"

    print("[1/4] Preparing Xenium breast WTA dataset...", flush=True)
    adata = sp.read_xenium(
        input_path,
        cell_group_path=cell_group_path if cell_group_path.exists() else None,
        roi_geojson_path=roi_geojson_path if roi_geojson_path.exists() else None,
        sample_name="xenium_wta_breast",
    )
    SPATIAL_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(SPATIAL_OUTPUT)

    print("[2/4] Fetching and preparing breast-specific CROP-seq reference...", flush=True)
    sp.fetch_dataset("gse241115_breast_cropseq", cache_dir=CACHE_DIR)
    sp.prepare_dataset("gse241115_breast_cropseq", cache_dir=CACHE_DIR)

    print("[3/4] Fetching and preparing pathway Perturb-seq atlas...", flush=True)
    sp.fetch_dataset("gse281048_pathway_atlas", cache_dir=CACHE_DIR)
    sp.prepare_dataset("gse281048_pathway_atlas", cache_dir=CACHE_DIR)

    print("[4/4] Running reference projection benchmark...", flush=True)
    results = sp.run_reference_projection_benchmark(
        SPATIAL_OUTPUT,
        reference_datasets=["gse241115_breast_cropseq", "gse281048_pathway_atlas"],
        config={
            "cache_dir": CACHE_DIR,
            "cell_group_path": cell_group_path if cell_group_path.exists() else None,
            "roi_geojson_path": roi_geojson_path if roi_geojson_path.exists() else None,
            "sample_name": "xenium_wta_breast",
            "pathway_cell_line": "MCF7",
            "k": 15,
        },
        output_dir=REPORT_DIR,
    )

    manifest_path = REPORT_DIR / "manifest.json"
    manifest = results.get("manifest", {})
    print(json.dumps({"report_dir": str(REPORT_DIR), "manifest": str(manifest_path), "summary": manifest.get("summary", {})}, indent=2))


if __name__ == "__main__":
    main()
