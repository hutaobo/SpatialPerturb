# Benchmarks

SpatialPerturb 当前把 benchmark 固定成两条主轨道：

- `shen_2026_core`
  目标是复现空间扰动数据上的 intrinsic / neighbor / ligand-receptor / power / figure 主链。
- `cross_platform_concordance`
  目标是比较 spatial 和 dissociated reference 中的 perturbation signatures 与 programs。

## 查看 catalog

```python
import spatialperturb as sp

sp.available_datasets()
sp.available_benchmarks()
```

## Public benchmark backbone

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
