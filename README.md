# SpatialPerturb

SpatialPerturb is a Python toolkit for combining **Spatial Transcriptomics** with **Perturb-seq** workflows. It is built around `AnnData`/`scanpy` objects and focuses on perturbation signatures, reference projection, spatial program scoring, neighborhood summaries, and benchmark-ready reports.

The current package supports two complementary analysis modes:

- Spatial perturbation analysis when the spatial assay already contains perturbation labels.
- Reference projection when the spatial assay is unperturbed tissue, such as Xenium WTA, and Perturb-seq is used as an external reference atlas.

## Highlights

- Stable SpatialPerturb `AnnData` schema with provenance in `uns["spatialperturb"]`.
- Xenium reader for real 10x outs using `cell_feature_matrix.h5` and `cells.csv.gz`.
- Optional Xenium cell-group CSV and Xenium Explorer GeoJSON ROI annotation.
- Public dataset lifecycle helpers: `fetch_dataset()`, `prepare_dataset()`, `load_public_dataset()`.
- Intrinsic DE, neighborhood DE, ligand-receptor scoring, power curves, and platform concordance.
- Perturb-seq reference program construction and spatial projection onto unperturbed tissue.
- A100-ready scripts for breast Xenium WTA + Perturb-seq reference projection.
- Reproducible report directories with tables, heatmaps, manifests, and biological interpretation notes.

## Install

```bash
pip install SpatialPerturb
```

For heavier ecosystem interop:

```bash
pip install "SpatialPerturb[interop]"
```

For development:

```bash
python -m pip install -e ".[interop]"
pytest -q
```

## Quick Start

```python
import spatialperturb as sp

adata = sp.load_demo_dataset()

intrinsic = sp.intrinsic_de(
    adata,
    perturbation="Lrrk2",
    control="control",
    method="pseudobulk",
    sample_col="sample",
    cell_type="neuron",
    roi="hippocampus",
)

neighbor = sp.neighbor_de(
    adata,
    perturbation="Lrrk2",
    control="control",
    method="pseudobulk",
    sample_col="sample",
    aggregate="pseudobulk",
    cell_type="neuron",
    roi="hippocampus",
)

programs = sp.derive_perturbation_programs(intrinsic, top_n=50, direction="both")
scores = sp.score_programs(adata, programs)
```

## Xenium WTA + Perturb-seq Reference Projection

SpatialPerturb can project Perturb-seq-derived programs onto an unperturbed Xenium tissue sample. This is useful when the tissue is not genetically perturbed but you want to ask which perturbation-like transcriptional programs appear in spatial context.

```python
import spatialperturb as sp

spatial = sp.read_xenium(
    "/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast",
    cell_group_path="/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv",
    roi_geojson_path="/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast/xenium_explorer_annotations.geojson",
    sample_name="xenium_wta_breast",
)

results = sp.run_reference_projection_benchmark(
    spatial,
    reference_datasets=["gse241115_breast_cropseq"],
    config={
        "cache_dir": "/data/taobo.hu/SpatialPerturb/cache",
        "k": 15,
        "reference_effect_size_only": True,
    },
    output_dir="/data/taobo.hu/SpatialPerturb/reports/breast_reference_projection",
)
```

The score means **transcriptional similarity to a Perturb-seq program**. It does not mean the tissue cell has the corresponding knockout, guide, or drug perturbation.

## Registered Datasets

```python
import spatialperturb as sp

sp.available_datasets()
```

Key dataset cards:

- `gse241115_breast_cropseq`: breast cancer stem-like state CROP-seq reference; automatic GEO `RAW.tar` parsing into `AnnData`.
- `gse281048_pathway_atlas`: pathway Perturb-seq atlas with MCF7 support; R/Seurat-backed preparation when `Rscript` and Seurat are available.
- `shen_2026_scrnaseq`: `GSE274058` scRNA-seq Perturb-seq-style reference; automatic fetch/prepare/load.
- `shen_2026_stereoseq`: `GSE274447` spatial track; automatic raw fetch/extraction with final conversion expected from GEF or tabular export.
- `demo_spatialperturb`: deterministic paired demo data.

## CLI

```bash
spatialperturb datasets
spatialperturb benchmarks
spatialperturb fetch-dataset gse241115_breast_cropseq
spatialperturb prepare-dataset gse241115_breast_cropseq
spatialperturb prepare-xenium /path/to/xenium_outs /path/to/xenium_wta_breast.h5ad
spatialperturb run-reference-benchmark /path/to/xenium_wta_breast.h5ad /path/to/report
spatialperturb run-benchmark demo_spatialperturb --output-dir reports/demo
spatialperturb validate path/to/data.h5ad
```

## A100 Breast Workflow

The repository includes scripts for the A100 environment used for the Xenium WTA breast reference projection run:

```bash
pwsh scripts/a100_sync_xenium_minimal.ps1
ssh sscb-a100.scilifelab.se
tmux new -d -s sp_breast_ref \
  "bash /data/taobo.hu/SpatialPerturb/code/SpatialPerturb/scripts/a100_run_breast_reference_projection.sh 2>&1 | tee /data/taobo.hu/SpatialPerturb/reports/breast_reference_projection/run.log"
bash /data/taobo.hu/SpatialPerturb/code/SpatialPerturb/scripts/a100_monitor_breast_reference_projection.sh --watch
```

Expected output directory:

```text
/data/taobo.hu/SpatialPerturb/reports/breast_reference_projection
```

Important outputs include `manifest.json`, `biological_interpretation.md`, `tables/program_scores_by_group.tsv`, `tables/neighbor_program_scores_by_group.tsv`, cell-level score tables, and heatmaps.

## Package Layout

- `spatialperturb.io`: AnnData ingestion helpers, including Xenium and Stereo-seq style readers.
- `spatialperturb.pp`: perturbation assignment and QC.
- `spatialperturb.gr`: spatial graph construction and neighbor collection.
- `spatialperturb.tl`: intrinsic DE, neighbor DE, ligand-receptor scoring, concordance, and power.
- `spatialperturb.signatures`: perturbation program derivation, scoring, aggregation, and neighborhood program scoring.
- `spatialperturb.datasets`: dataset registry plus public `fetch/prepare/load`.
- `spatialperturb.benchmarks`: benchmark orchestration and report manifests.
- `spatialperturb.reports`: fixed paper figure rendering.

## Documentation

ReadTheDocs builds from `docs/` using Sphinx/MyST. The main workflow pages are:

- `docs/workflow.md`
- `docs/benchmarks.md`
- `docs/breast-reference-projection.md`
- `docs/api.md`

## Citation

Please cite the package if you find it useful. See `CITATION.cff`.
