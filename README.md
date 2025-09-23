# SpatialPerturb

Toolkit for combining **Spatial Transcriptomics** with **Perturb-seq** workflows — signatures, label transfer, spatial scoring, and graph/structure analysis.

## Install

```bash
pip install SpatialPerturb
# or with GNN extras:
pip install 'SpatialPerturb[gnn]'
```

## Quick start

```python
import spatialperturb as sp

print(sp.__version__)
```

CLI:
```bash
SpatialPerturb version
```

## Development

```bash
python -m pip install --upgrade build twine
python -m build
twine upload --repository testpypi dist/*
```

## Citation

Please cite the package if you find it useful. See `CITATION.cff`.
