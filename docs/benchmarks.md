# Benchmarks

SpatialPerturb 当前把 benchmark 固定成两条主轨道：

- `shen_2026_core`
  目标是复现空间扰动数据上的 intrinsic / neighbor / ligand-receptor / power / figure 主链。
- `cross_platform_concordance`
  目标是比较 spatial 和 dissociated reference 中的 perturbation signatures 与 programs。
- `breast_reference_projection`
  目标是把 breast Perturb-seq reference programs 投影到未扰动 Xenium WTA tissue，并输出 cell-level、cell-type/ROI-level 和 neighborhood-level program scores。

## 查看 catalog

```python
import spatialperturb as sp

sp.available_datasets()
sp.available_benchmarks()
```

## Public benchmark backbone

### `gse241115_breast_cropseq`

- accession: `GSE241115`
- role: primary breast cancer CROP-seq reference for reference projection
- raw format: flat GEO `RAW.tar` with 10x `mtx/tsv` files and `protospacer_calls_per_cell.csv.gz`
- status: automatic `fetch -> prepare -> load` supported
- note: sgRNA / intergenic guide features are tracked as `barcode_columns` and excluded from expression DE/program genes

### `gse281048_pathway_atlas`

- accession: `GSE281048`
- role: optional pathway Perturb-seq atlas; default downstream filter is `cell_line == "MCF7"`
- raw format: Seurat `.rds.gz`
- status: automatic fetch supported; prepare requires `Rscript` and Seurat

### `shen_2026_scrnaseq`

- accession: `GSE274058`
- role: reference / cross-platform track
- raw format: nested `10x tar.gz`
- status: automatic `fetch -> prepare -> load` supported

### `shen_2026_stereoseq`

- accession: `GSE274447`
- role: spatial core track
- raw format: `tar of GEF`
- status: automatic fetch and extraction supported; final prepare still expects a preconverted `.h5ad` or tabular cell-level export

## 运行 core benchmark

```python
import spatialperturb as sp

results = sp.run_core_benchmark(
    "demo_spatialperturb",
    config={
        "cache_dir": ".spatialperturb-cache",
        "method": "pseudobulk",
        "sample_col": "sample",
        "reference_dataset": "demo_spatialperturb",
        "concordance_level": "both",
    },
    output_dir="reports/demo_spatialperturb",
)
```

这个入口会自动：

- 载入 prepared dataset
- 补 spatial graph（如果还没建）
- 运行 `intrinsic_de`
- 运行 `neighbor_de`
- 运行 `differential_lr`
- 运行 `power_curve`
- 如果给了 reference，再运行 `platform_concordance`
- 输出 tables、figures、`manifest.json` 和 `input.h5ad`

## 运行 cross-platform benchmark

```python
spatial, reference = sp.load_demo_dataset(paired=True)

spatial_de = sp.intrinsic_de(
    spatial,
    perturbation="Lrrk2",
    control="control",
    method="pseudobulk",
    sample_col="sample",
)

reference_de = sp.intrinsic_de(
    reference,
    perturbation="Lrrk2",
    control="control",
    method="pseudobulk",
    sample_col="sample",
)

concordance = sp.run_cross_platform_benchmark(
    spatial_de,
    reference_de,
    config={"top_n": 50, "level": "both"},
)
```

## 运行 breast reference projection benchmark

```python
import spatialperturb as sp

results = sp.run_reference_projection_benchmark(
    "/data/taobo.hu/SpatialPerturb/prepared/xenium_wta_breast.h5ad",
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

输出包括：

- `tables/program_scores_cell_level.tsv.gz`
- `tables/program_scores_by_group.tsv`
- `tables/neighbor_program_scores_cell_level.tsv.gz`
- `tables/neighbor_program_scores_by_group.tsv`
- `tables/reference_de.tsv`
- `tables/reference_program_membership.tsv`
- `figures/program_scores_heatmap.png`
- `figures/neighbor_program_scores_heatmap.png`
- `manifest.json`
- `biological_interpretation.md`

在 full-scale runs 中，`reference_effect_size_only=True` 会用 log2 fold-change 排名构建 programs；这适合 program projection，但 `reference_de.tsv` 中的 p-value/FDR 不应当用于显著性声明。

## Benchmark 输出目录

`run_core_benchmark(..., output_dir=...)` 会生成固定目录结构：

- `tables/intrinsic_de.tsv`
- `tables/neighbor_de.tsv`
- `tables/differential_lr.tsv`
- `tables/power_curve.tsv`
- `tables/platform_concordance.tsv`（如果提供 reference）
- `figures/workflow_schema.png`
- `figures/assignment_qc.png`
- `figures/own_vs_neighbor.png`
- `figures/lr_differential.png`
- `figures/platform_concordance.png`
- `figures/power_curve.png`
- `manifest.json`
- `config.json`
- `input.h5ad`
