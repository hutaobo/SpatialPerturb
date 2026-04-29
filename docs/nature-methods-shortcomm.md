# Nature Methods Brief Communication Workflow

This workflow turns the breast Xenium WTA reference projection into a submission-grade analysis bundle.

## What It Adds

- Held-out Perturb-seq reference validation with AUROC/AUPRC.
- Expression-matched random programs and group-label shuffles for empirical calibration.
- Ranked candidate filtering with `min_claim_cells=50`; smaller ROI/cell-type groups remain supplementary.
- Bootstrap confidence intervals for group-level program scores.
- Moran-style spatial autocorrelation on the Xenium spatial graph.
- Ablations for `top_n = 25, 50, 100` and graph `k = 5, 15, 30`.
- Two main figures plus a Nature Methods short-communication scaffold.

## Python API

```python
import spatialperturb as sp

results = sp.run_nature_methods_breast_analysis(
    "/data/taobo.hu/SpatialPerturb/prepared/xenium_wta_breast.h5ad",
    reference_datasets=["gse241115_breast_cropseq", "gse281048_pathway_atlas"],
    config={
        "cache_dir": "/data/taobo.hu/SpatialPerturb/cache",
        "pathway_cell_line": "MCF7",
        "reference_effect_size_only": True,
        "n_random": 25,
        "n_spatial_permutations": 25,
        "n_bootstrap": 100,
        "min_claim_cells": 50,
    },
    output_dir="/data/taobo.hu/SpatialPerturb/reports/nature_methods_breast_shortcomm",
)
```

## CLI

```bash
spatialperturb run-nature-methods-breast-analysis \
  /data/taobo.hu/SpatialPerturb/prepared/xenium_wta_breast.h5ad \
  /data/taobo.hu/SpatialPerturb/reports/nature_methods_breast_shortcomm \
  --cache-dir /data/taobo.hu/SpatialPerturb/cache
```

## A100 Run

```bash
tmux new -d -s sp_nm_breast \
  "bash /data/taobo.hu/SpatialPerturb/code/SpatialPerturb/scripts/a100_run_nature_methods_breast_analysis.sh 2>&1 | tee /data/taobo.hu/SpatialPerturb/reports/nature_methods_breast_shortcomm/run.log"

bash /data/taobo.hu/SpatialPerturb/code/SpatialPerturb/scripts/a100_monitor_nature_methods_breast_analysis.sh --watch
```

The A100 runner first tries to create `/data/taobo.hu/SpatialPerturb/envs/r-seurat` with micromamba, mamba, or conda. For the submission run, `gse281048_pathway_atlas` is expected to be prepared and filtered to MCF7 through that R/Seurat environment; if a development environment lacks R, the manifest records the secondary-reference failure explicitly without changing the primary GSE241115 analysis.

The current submission run treats both references as ready: GSE241115 contributes the breast CROP-seq baseline and GSE281048 contributes 218 MCF7 pathway programs, expanding the reference panel from 50 to 268 programs.

## Key Outputs

- `manifest.json`
- `nature_methods_summary.md`
- `nature_methods_shortcomm_scaffold.md`
- `biological_interpretation.md`
- `figures/main_figure_1.png`
- `figures/main_figure_2.png`
- `tables/reference_validation.tsv`
- `tables/calibrated_program_scores_by_group.tsv`
- `tables/spatial_autocorrelation.tsv`
- `tables/ablation_summary.tsv`

## Production Panel Output

The submission-facing panel layer lives in `spatialperturb.figurekit`. It is separate from exploratory plotting helpers in `spatialperturb.pl` and enforces production defaults for Nature-style figure assembly:

- PDF/PS fonts are exported as Type 42 and SVG text is left editable.
- Panel sizes are declared in millimeters through `PanelSpec`.
- Axis labels, ticks, legends and panel labels use fixed small-format typography.
- Dense spatial scatter layers are rasterized while axes, labels, legends and colorbars remain vector.
- Every saved panel exports `.pdf`, `.png` and `.svg`, writes a matching Source Data `.tsv` file and records an entry in `panel_manifest.tsv`.

```python
import spatialperturb as sp

outputs = sp.render_nature_methods_panels(
    "/data/taobo.hu/SpatialPerturb/reports/nature_methods_breast_shortcomm",
    "/data/taobo.hu/SpatialPerturb/manuscripts/nature_methods_shortcomm",
    strict=True,
)
```

This writes eight independent panel files under `panels/` and their matched source tables under `source_data/`:

- `fig1a_workflow_schema`
- `fig1b_reference_validation`
- `fig1c_null_calibration`
- `fig1d_ablation_robustness`
- `fig2a_xenium_map_celltype_roi`
- `fig2b_top_program_spatial_map`
- `fig2c_roi_celltype_heatmap`
- `fig2d_spatial_autocorrelation`

The main panels intentionally use a curated short-communication selection rather than all 268 programs. Figure 1 emphasizes reference readiness, `50 -> 268` program expansion, GSE281048 validation and null-calibration caveats. Figure 2 focuses on ranked candidate spatial programs with interpretable breast-cancer context: Mast-cell `FOS`, basal-like structured DCIS `CEBPB`, dendritic `SP1`/`MTOR`/`RPS6KB1`/`MAPK3`, luminal-like amorphous DCIS `PTGS2`, invasive-associated CAF `MAPK8`, and mitotic invasive tumor `IFNAR1`/`TYK2`. Full FDR values, bootstrap intervals, ablations, redundancy and all 268 programs belong in Source Data, supplementary tables or Extended Data.

## Interpretation Guardrail

Projection scores quantify Perturb-seq reference-like transcriptional states. They should be framed as spatially localized candidate regulatory states, not as proof that the tissue cells underwent the corresponding knockout, CRISPRi perturbation, pathway intervention, or drug treatment.
