# SpatialPerturb 文档

SpatialPerturb 是一个面向空间扰动转录组分析的 `AnnData` 原生 Python 包，目标是把空间 Perturb-seq 的主分析链固化成可复现的方法框架，而不是只停留在一次性分析脚本。

## 当前版本可以做什么

- 统一 `AnnData` schema，并在 `uns["spatialperturb"]` 里保留 provenance。
- 导入 Xenium / Stereo-seq 风格 cell-level 数据，或直接从 expression/meta tables 构建对象。
- 基于 barcode features 做 perturbation assignment，并自动标记 `single` / `multiple` / `unassigned`。
- 构建 `knn` / `radius` 空间图，收集邻居边。
- 运行 `intrinsic_de`、`neighbor_de`、`differential_lr`、`platform_concordance`、`power_curve`。
- 支持 `simple` 和 `pseudobulk` 两种分析模式。
- 提供 public dataset lifecycle：`fetch_dataset()`、`prepare_dataset()`、`load_public_dataset()`。
- 跑 benchmark 并固定导出六类论文图和 manifest。

## 最短示例

```python
import spatialperturb as sp

adata = sp.load_demo_dataset()

results = sp.run_core_benchmark(
    adata,
    perturbations=["Lrrk2", "Srf"],
    control="control",
    config={
        "method": "pseudobulk",
        "sample_col": "sample",
        "concordance_level": "both",
    },
    output_dir="reports/demo",
)
```

## 公开 benchmark 目录

- `shen_2026_stereoseq` -> `GSE274447`
- `shen_2026_scrnaseq` -> `GSE274058`
- `demo_spatialperturb` -> synthetic paired dataset

更多细节见 [Workflow](workflow.md)、[Benchmarks](benchmarks.md)、[Paper Repro](paper-repro.md) 和 [API 参考](api.md)。

```{toctree}
:maxdepth: 1
:caption: 文档目录

workflow
benchmarks
paper-repro
gse274058-reference-results
api
```
