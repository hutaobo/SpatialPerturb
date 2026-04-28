# Workflow

SpatialPerturb 的标准工作流现在固定为四段：

1. `fetch -> prepare -> load`
2. perturbation assignment / schema validation / graph construction
3. intrinsic / neighbor / LR / concordance / power analysis，或 reference program projection
4. paper-style figure rendering、score tables、manifest 和 interpretation export

## 1. 从公开数据开始

```python
import spatialperturb as sp

sp.available_datasets()

sp.fetch_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
sp.prepare_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
adata = sp.load_public_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
```

如果你要跑 `shen_2026_stereoseq`，当前版本支持自动下载和解压 raw GEF，但最后一步仍需要你先把 raw GEF 转成 `.h5ad` 或 tabular cell-level export，再放回 dataset raw 目录重新执行 `prepare_dataset()`。

Breast reference projection 相关数据集：

```python
sp.fetch_dataset("gse241115_breast_cropseq", cache_dir=".spatialperturb-cache")
sp.prepare_dataset("gse241115_breast_cropseq", cache_dir=".spatialperturb-cache")
reference = sp.load_public_dataset("gse241115_breast_cropseq", cache_dir=".spatialperturb-cache")
```

`gse281048_pathway_atlas` 需要 `Rscript` 和 Seurat 来转换 `.rds.gz` Seurat object；如果运行环境没有 R/Seurat，A100 workflow 会把它标记为 optional blocked，而不是让主分析失败。

## 2. 从自己的 cell-level 数据开始

```python
import spatialperturb as sp

adata = sp.from_tables(
    expression_df,
    obs=cell_metadata,
    spatial=cell_metadata[["x", "y"]],
    metadata={"platform": "xenium"},
)

sp.assign_perturbations(
    adata,
    barcode_columns=["CTRL_BARCODE", "LRRK2_BARCODE", "SRF_BARCODE"],
    barcode_to_perturbation={
        "CTRL_BARCODE": "control",
        "LRRK2_BARCODE": "Lrrk2",
        "SRF_BARCODE": "Srf",
    },
)

sp.build_spatial_graph(adata, mode="knn", k=15)
```

### 2b. 从真实 Xenium outs 开始

```python
adata = sp.read_xenium(
    "/path/to/xenium_outs",
    cell_group_path="/path/to/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv",
    roi_geojson_path="/path/to/xenium_explorer_annotations.geojson",
    sample_name="xenium_wta_breast",
)
```

`read_xenium()` 会优先读取 10x `cell_feature_matrix.h5` 和 `cells.csv.gz`，并把 cell-group CSV 中的 `group` 合并到 `obs["cell_type"]`。ROI GeoJSON 使用 cell centroid 做 point-in-polygon，未命中多边形的细胞保留 `roi="global"`。

## 3. 统计分析

轻量模式适合 demo、小样本和快速探索：

```python
intrinsic = sp.intrinsic_de(
    adata,
    perturbation="Lrrk2",
    control="control",
    method="simple",
    cell_type="neuron",
    roi="hippocampus",
)
```

论文默认建议用样本级 `pseudobulk`：

```python
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
    drop_shared_neighbors=False,
    weight_by_distance=False,
    cell_type="neuron",
    roi="hippocampus",
)

lr = sp.differential_lr(
    adata,
    perturbation="Lrrk2",
    control="control",
    lr_network="fallback",
)

power = sp.power_curve(
    adata,
    perturbation="Lrrk2",
    control="control",
    method="pseudobulk",
    sample_col="sample",
)
```

## 4. Program 和 cross-platform concordance

```python
programs = sp.derive_perturbation_programs(intrinsic, top_n=50, direction="both")
scores = sp.score_programs(adata, programs)

concordance = sp.platform_concordance(
    spatial_results,
    reference_results,
    top_n=50,
    level="both",
)
```

`level="both"` 会同时返回：

- gene-level correlation
- top-gene overlap
- program-level Jaccard concordance

## 5. Reference projection 到未扰动组织

当空间样本没有真实 perturbation label 时，可以把 Perturb-seq reference programs 投影到 Xenium WTA tissue：

```python
results = sp.run_reference_projection_benchmark(
    adata,
    reference_datasets=["gse241115_breast_cropseq"],
    config={
        "cache_dir": ".spatialperturb-cache",
        "k": 15,
        "groupby": ["cell_type", "roi"],
        "reference_effect_size_only": True,
    },
    output_dir="reports/breast_reference_projection",
)
```

这个入口会自动：

- 构建 Xenium `knn` spatial graph。
- 从 Perturb-seq reference 里构建每个 perturbation 的 top gene program。
- 计算 cell-level program scores。
- 计算 neighborhood program scores。
- 按 `cell_type` 和 `roi` 聚合。
- 写出 heatmaps、score tables、reference DE table 和 `manifest.json`。

解释时需要注意：projection score 表示“空间细胞表达状态与 reference perturbation program 相似”，不能解释为组织里发生了真实 knockout 或药物扰动。

## 6. 导出论文图

```python
results = {
    "adata": adata,
    "intrinsic_de": intrinsic,
    "neighbor_de": neighbor,
    "differential_lr": lr,
    "platform_concordance": concordance,
    "power_curve": power,
}

sp.render_paper_figures(results, output_dir="reports/figures")
```

固定输出六类图：

- workflow/schema
- perturbation assignment QC
- own-vs-neighbor
- ligand-receptor differential
- cross-platform concordance
- power and sensitivity
