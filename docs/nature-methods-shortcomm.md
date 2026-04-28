# Nature Methods Brief Communication Workflow

This workflow turns the breast Xenium WTA reference projection into a submission-grade analysis bundle.

## What It Adds

- Held-out Perturb-seq reference validation with AUROC/AUPRC.
- Expression-matched random programs and group-label shuffles for empirical calibration.
- Claim-level filtering with `min_claim_cells=50`; smaller ROI/cell-type groups remain supplementary.
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

The A100 runner first tries to create `/data/taobo.hu/SpatialPerturb/envs/r-seurat` with micromamba, mamba, or conda. If R/Seurat is available, `gse281048_pathway_atlas` is prepared and filtered to MCF7; otherwise the primary GSE241115 analysis still runs and the optional reference is recorded as blocked in the manifest.

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

## Interpretation Guardrail

Projection scores quantify transcriptional similarity to Perturb-seq-derived programs. They should be framed as spatially localized candidate regulatory states, not as proof that the tissue cells underwent the corresponding knockout, CRISPRi perturbation, or drug treatment.
