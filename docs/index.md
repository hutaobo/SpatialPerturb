# SpatialPerturb 文档

SpatialPerturb 是一个面向 Spatial Transcriptomics + Perturb-seq 工作流的 `AnnData` 原生 Python 包。它既支持已有 perturbation label 的空间扰动分析，也支持“未扰动 Xenium 组织 + Perturb-seq reference atlas”的 program projection 分析。

## 当前版本可以做什么

- 统一 `AnnData` schema，并在 `uns["spatialperturb"]` 里保留 provenance。
- 读取真实 10x Xenium outs，包括 `cell_feature_matrix.h5`、`cells.csv.gz`、cell-group CSV 和 Xenium Explorer ROI GeoJSON。
- 导入 Stereo-seq 风格 cell-level 数据，或直接从 expression/meta tables 构建对象。
- 基于 barcode features 做 perturbation assignment，并自动标记 `single` / `multiple` / `unassigned`。
- 构建 `knn` / `radius` 空间图，收集邻居边。
- 运行 `intrinsic_de`、`neighbor_de`、`differential_lr`、`platform_concordance`、`power_curve`。
- 构建 Perturb-seq reference programs，并投影到 unperturbed spatial tissue。
- 支持 `simple`、`pseudobulk` 和 full-scale reference projection 中的 effect-size-only program ranking。
- 支持 Nature Methods Brief Communication 风格的 reference validation、null calibration、spatial autocorrelation、ablation 和主图输出。
- 提供 public dataset lifecycle：`fetch_dataset()`、`prepare_dataset()`、`load_public_dataset()`。
- 跑 benchmark 并固定导出 tables、figures、manifest 和 biological interpretation。

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

## 公开 benchmark / reference 目录

- `gse241115_breast_cropseq` -> breast cancer CROP-seq reference
- `gse281048_pathway_atlas` -> pathway Perturb-seq atlas, including MCF7 when R/Seurat is available
- `shen_2026_stereoseq` -> `GSE274447`
- `shen_2026_scrnaseq` -> `GSE274058`
- `demo_spatialperturb` -> synthetic paired dataset

更多细节见 [Workflow](workflow.md)、[Benchmarks](benchmarks.md)、[Breast Reference Projection](breast-reference-projection.md)、[Nature Methods Short Communication](nature-methods-shortcomm.md)、[Paper Repro](paper-repro.md) 和 [API 参考](api.md)。

```{toctree}
:maxdepth: 1
:caption: 文档目录

workflow
benchmarks
breast-reference-projection
nature-methods-shortcomm
paper-repro
gse274058-reference-results
api
```
