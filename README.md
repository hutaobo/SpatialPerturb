# SpatialPerturb

Toolkit for combining **Spatial Transcriptomics** with **Perturb-seq** workflows — signatures, label transfer, spatial scoring, and graph/structure analysis.

The repository currently provides a lightweight, documented package foundation with signature-matrix utilities and a CLI entrypoint. Additional analysis modules can grow on top of this base.

## Install

```bash
pip install SpatialPerturb
# or with GNN extras:
pip install 'SpatialPerturb[gnn]'
```

## Quick start

```python
import spatialperturb as sp
from spatialperturb import build_signature_matrix

print(sp.__version__)

gene_sets = {
    "IFN_response": ["STAT1", "IRF1", "CXCL10"],
    "Cell_cycle": ["MKI67", "TOP2A"],
}

signature_matrix = build_signature_matrix(gene_sets)
print(signature_matrix)
```

CLI:
```bash
spatialperturb version
```

## Development

```bash
python -m pip install --upgrade build twine pytest
python -m build
pytest -q
twine upload --repository testpypi dist/*
```

## Citation

Please cite the package if you find it useful. See `CITATION.cff`.
