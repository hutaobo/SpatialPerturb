# Workflow

SpatialPerturb 的标准工作流现在固定为四段：

1. `fetch -> prepare -> load`
2. perturbation assignment / schema validation / graph construction
3. intrinsic / neighbor / LR / concordance / power analysis
4. paper-style figure rendering and manifest export

## 1. 从公开数据开始

```python
import spatialperturb as sp

sp.available_datasets()

sp.fetch_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
sp.prepare_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
adata = sp.load_public_dataset("shen_2026_scrnaseq", cache_dir=".spatialperturb-cache")
```

如果你要跑 `shen_2026_stereoseq`，当前版本支持自动下载和解压 raw GEF，但最后一步仍需要你先把 raw GEF 转成 `.h5ad` 或 tabular cell-level export，再放回 dataset raw 目录重新执行 `prepare_dataset()`。

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

## 5. 导出论文图

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
