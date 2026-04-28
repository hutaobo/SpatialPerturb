# Paper Repro

这一页对应论文级复现包的固定入口，目标是把“从公开数据到主图和主表”变成一个稳定工作流。

## CLI workflow

### 1. 查看数据集

```bash
spatialperturb datasets
```

### 2. 下载和准备数据

```bash
spatialperturb fetch-dataset shen_2026_scrnaseq --cache-dir .spatialperturb-cache
spatialperturb prepare-dataset shen_2026_scrnaseq --cache-dir .spatialperturb-cache
```

对于 `shen_2026_stereoseq`，当前版本支持自动下载和解压，但如果你还没有把 GEF 转换成 `.h5ad` 或 tabular cell-level export，`prepare-dataset` 会在 dataset cache 目录写出说明文件，而不是伪造一个准备完成的对象。

### 3. 运行 benchmark

```bash
spatialperturb run-benchmark demo_spatialperturb \
  --cache-dir .spatialperturb-cache \
  --output-dir reports/demo_spatialperturb
```

如果你需要同时生成 cross-platform concordance，可以给 `run_core_benchmark()` 传 `reference_dataset` 或 `reference_adata`。

### 4. 重渲染主图

```bash
spatialperturb render-paper-figures demo_spatialperturb \
  --cache-dir .spatialperturb-cache \
  --output-dir reports/demo_spatialperturb
```

## Python workflow

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

## 固定输出

`output_dir` 下会包含：

- `tables/`: 所有 tidy benchmark tables
- `figures/`: 六张主图
- `manifest.json`: 运行时间、版本、配置、文件索引
- `config.json`: 参数快照
- `input.h5ad`: 本次运行实际使用的数据对象

## 推荐论文主图

当前固定导出：

1. `workflow_schema.png`
2. `assignment_qc.png`
3. `own_vs_neighbor.png`
4. `lr_differential.png`
5. `platform_concordance.png`
6. `power_curve.png`

## 当前边界

- `GSE274058` 已支持自动准备。
- `GSE274447` 目前仍需要你先做 GEF 到 cell-level AnnData 的转换。
- LR 分析默认使用内置 deterministic fallback network；正式 benchmark 最好传入你自己冻结版本的 curated LR table。
