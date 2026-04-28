# SpatialPerturb

SpatialPerturb is an AnnData-native framework for spatial perturbation inference across sequencing-based and imaging-based platforms.

It now ships a benchmark-oriented workflow built around:

- a stable `AnnData` schema,
- `fetch -> prepare -> load` public dataset lifecycle helpers,
- intrinsic and neighborhood differential effects with `simple` and `pseudobulk` modes,
- ligand-receptor differential scoring with fixed fallback or custom LR resources,
- perturbation-program and cross-platform concordance metrics,
- paper-style figure rendering and report manifests.

## Install

```bash
pip install SpatialPerturb
```

For heavier ecosystem interop:

```bash
pip install "SpatialPerturb[interop]"
```

## Quick start

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

lr = sp.differential_lr(adata, perturbation="Lrrk2", control="control", lr_network="fallback")
power = sp.power_curve(adata, perturbation="Lrrk2", control="control", method="pseudobulk", sample_col="sample")
programs = sp.derive_perturbation_programs(intrinsic, top_n=10, direction="both")
```

## Public dataset lifecycle

```python
import spatialperturb as sp

sp.available_datasets()

sp.fetch_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
sp.prepare_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
adata = sp.load_public_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
```

Registered public tracks:

- `shen_2026_stereoseq` -> `GSE274447`
- `shen_2026_scrnaseq` -> `GSE274058`
- `demo_spatialperturb` -> deterministic paired demo data

Notes:

- `shen_2026_scrnaseq` supports automatic fetch and preparation from the GEO raw archive.
- `shen_2026_stereoseq` supports automatic fetch and extraction, but final preparation still requires a preconverted `.h5ad` or tabular export from the raw GEF files.

## Paper-grade benchmark workflow

```python
import spatialperturb as sp

results = sp.run_core_benchmark(
    "demo_spatialperturb",
    config={
        "cache_dir": ".spatialperturb-cache",
        "reference_dataset": "demo_spatialperturb",
        "method": "pseudobulk",
        "sample_col": "sample",
        "concordance_level": "both",
    },
    output_dir="reports/demo_spatialperturb",
)
```

This writes:

- tidy tables under `reports/.../tables/`
- fixed paper figures under `reports/.../figures/`
- a machine-readable `manifest.json`
- the exact `input.h5ad` used for the run

## CLI

```bash
spatialperturb datasets
spatialperturb fetch-dataset shen_2026_scrnaseq
spatialperturb prepare-dataset shen_2026_scrnaseq
spatialperturb run-benchmark demo_spatialperturb --output-dir reports/demo
spatialperturb render-paper-figures demo_spatialperturb --output-dir reports/demo-figs
spatialperturb validate path/to/data.h5ad
```

## Package layout

- `spatialperturb.io`: AnnData ingestion helpers.
- `spatialperturb.pp`: perturbation assignment and QC.
- `spatialperturb.gr`: spatial graph construction and neighbor collection.
- `spatialperturb.tl`: intrinsic DE, neighbor DE, ligand-receptor scoring, concordance, and power.
- `spatialperturb.pl`: plotting helpers for benchmark figures.
- `spatialperturb.signatures`: perturbation program derivation and scoring.
- `spatialperturb.datasets`: dataset registry plus public `fetch/prepare/load`.
- `spatialperturb.benchmarks`: benchmark orchestration and report manifests.
- `spatialperturb.reports`: fixed paper figure rendering.

## Development

```bash
python -m pip install --upgrade build pytest twine
python -m build
pytest -q
```

## Citation

Please cite the package if you find it useful. See `CITATION.cff`.
