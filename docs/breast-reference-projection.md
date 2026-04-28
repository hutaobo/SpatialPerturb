# Breast Reference Projection

This page describes the breast Xenium WTA + Perturb-seq reference projection workflow implemented in SpatialPerturb.

## Goal

The workflow upgrades the analysis from “spatial data already has perturbation labels” to:

```text
unperturbed Xenium WTA tissue + Perturb-seq reference atlas -> spatial perturbation-like program scores
```

The main biological question is not whether a tissue cell was perturbed. Instead, it asks whether spatial cell states resemble transcriptional programs learned from Perturb-seq references.

## Reference Datasets

### Primary reference: `gse241115_breast_cropseq`

- GEO accession: `GSE241115`
- Biological role: breast cancer CROP-seq reference.
- Prepared output: standard `AnnData`.
- Perturbation metadata:
  - `guide_id`: original sgRNA/protospacer call.
  - `perturbation`: guide name after removing `_sgRNA\d+` or `_sg\d+`.
  - `perturbation_status`: `single`, `multiple`, or `unassigned`.
  - controls are normalized to `control`.
- Guide/intergenic features are recorded as `uns["spatialperturb"]["barcode_columns"]` and excluded from expression DE and program gene selection.

### Optional reference: `gse281048_pathway_atlas`

- GEO accession: `GSE281048`
- Biological role: pathway Perturb-seq atlas, with `MCF7` used as the default breast cancer subset.
- Preparation requires `Rscript` and Seurat because the raw files are Seurat `.rds.gz` objects.
- If R/Seurat is unavailable, A100 scripts record `GSE281048_BLOCKED_RSCRIPT_MISSING` and continue with the primary reference.

## Xenium Input

`read_xenium()` supports real 10x Xenium outs:

```python
import spatialperturb as sp

adata = sp.read_xenium(
    "/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast",
    cell_group_path="/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv",
    roi_geojson_path="/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/xenium_explorer_annotations.geojson",
    sample_name="xenium_wta_breast",
)
```

The reader uses:

- `cell_feature_matrix.h5` for cell-by-gene counts.
- `cells.csv.gz` for centroids and cell metadata.
- optional cell-group CSV for `obs["cell_type"]`.
- optional Xenium Explorer GeoJSON ROI polygons for `obs["roi"]`.

ROI assignment uses `matplotlib.path.Path.contains_points`, so no `shapely` or `geopandas` dependency is required.

## Python API

```python
results = sp.run_reference_projection_benchmark(
    adata,
    reference_datasets=["gse241115_breast_cropseq"],
    config={
        "cache_dir": "/data/taobo.hu/SpatialPerturb/cache",
        "k": 15,
        "groupby": ["cell_type", "roi"],
        "reference_effect_size_only": True,
    },
    output_dir="/data/taobo.hu/SpatialPerturb/reports/breast_reference_projection",
)
```

Default strategy:

- Xenium graph: `mode="knn"`, `k=15`.
- Spatial aggregation: `cell_type` and `roi`.
- Primary reference: `gse241115_breast_cropseq`.
- Optional pathway reference: `gse281048_pathway_atlas`, filtered to `cell_line == "MCF7"` when available.
- Reference DE: `pseudobulk` when sufficient samples exist; otherwise `simple`.
- Full-scale A100 run: `reference_effect_size_only=True` to rank programs by log2 fold-change without expensive per-gene significance tests.

## CLI

```bash
spatialperturb prepare-xenium \
  /data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast \
  /data/taobo.hu/SpatialPerturb/prepared/xenium_wta_breast.h5ad \
  --cell-group-path /data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv \
  --roi-geojson-path /data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/xenium_explorer_annotations.geojson \
  --sample-name xenium_wta_breast

spatialperturb run-reference-benchmark \
  /data/taobo.hu/SpatialPerturb/prepared/xenium_wta_breast.h5ad \
  /data/taobo.hu/SpatialPerturb/reports/breast_reference_projection \
  --cache-dir /data/taobo.hu/SpatialPerturb/cache \
  --reference-datasets gse241115_breast_cropseq
```

## A100 Run

The A100 run uses:

```text
/data/taobo.hu/SpatialPerturb/cache
/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast
/data/taobo.hu/SpatialPerturb/prepared
/data/taobo.hu/SpatialPerturb/reports/breast_reference_projection
```

Launch:

```bash
tmux new -d -s sp_breast_ref \
  "bash /data/taobo.hu/SpatialPerturb/code/SpatialPerturb/scripts/a100_run_breast_reference_projection.sh 2>&1 | tee /data/taobo.hu/SpatialPerturb/reports/breast_reference_projection/run.log"

bash /data/taobo.hu/SpatialPerturb/code/SpatialPerturb/scripts/a100_monitor_breast_reference_projection.sh --watch
```

The monitor writes:

- `status.json`
- `status.md`

Expected final states:

- `COMPLETE`
- `COMPLETE_WITH_BLOCKED_OPTIONAL_REFERENCES`
- `FAILED`

## Outputs

The report directory contains:

- `input_spatial.h5ad`
- `references/gse241115_breast_cropseq.h5ad`
- `tables/program_scores_cell_level.tsv.gz`
- `tables/program_scores_by_group.tsv`
- `tables/neighbor_program_scores_cell_level.tsv.gz`
- `tables/neighbor_program_scores_by_group.tsv`
- `tables/reference_de.tsv`
- `tables/reference_program_membership.tsv`
- `tables/top_programs_by_roi_cell_type.tsv`
- `tables/top_neighbor_programs.tsv`
- `figures/program_scores_heatmap.png`
- `figures/neighbor_program_scores_heatmap.png`
- `manifest.json`
- `biological_interpretation.md`

## Biological Interpretation

Interpretation defaults:

- Cell-level score: transcriptional similarity to a Perturb-seq-derived program.
- ROI/cell-type aggregation: where a reference-like state localizes in the tissue.
- Neighbor score: whether adjacent cells show similar program activity.

For the breast Xenium WTA run, the strongest cleaned programs localized to 11q13 invasive tumor cell states and were enriched for luminal/secretory breast epithelial genes such as `TFF1`, `TFF3`, `AGR2`, `MUCL1`, `IGFBP5`, `DSCAM-AS1`, `TIMP3`, and `SPTSSB`.

Important caveats:

- Projection is not causal perturbation evidence.
- Perturb-seq references are cell-line-derived, while Xenium is FFPE tissue.
- ROI and cell-group annotation quality directly affects interpretation.
- Full-scale effect-size-only runs rank genes by log2 fold-change; do not use placeholder p-values/FDR for statistical claims.
