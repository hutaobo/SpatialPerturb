"""Internal helpers shared across SpatialPerturb modules."""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse


def to_dense(matrix: object) -> np.ndarray:
    """Convert sparse or array-like matrices to a dense NumPy array."""
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def resolve_var_names(adata: AnnData, var_names: Sequence[str] | None = None) -> list[str]:
    """Return feature names present in the object."""
    if var_names is None:
        return list(map(str, adata.var_names))
    missing = [name for name in var_names if name not in adata.var_names]
    if missing:
        raise KeyError(f"Features not found in AnnData.var_names: {missing}")
    return list(var_names)


def extract_matrix(
    adata: AnnData,
    *,
    layer: str | None = None,
    obs_names: Sequence[str] | None = None,
    var_names: Sequence[str] | None = None,
) -> np.ndarray:
    """Extract a dense expression matrix from an AnnData object."""
    obs_indexer = slice(None) if obs_names is None else list(obs_names)
    features = resolve_var_names(adata, var_names)
    view = adata[obs_indexer, features]
    matrix = view.layers[layer] if layer is not None else view.X
    return to_dense(matrix)


def benjamini_hochberg(pvalues: Iterable[float]) -> np.ndarray:
    """Compute Benjamini-Hochberg adjusted q-values."""
    pvalues = np.asarray(list(pvalues), dtype=float)
    if pvalues.size == 0:
        return pvalues

    order = np.argsort(pvalues)
    ranked = pvalues[order]
    n = ranked.size
    adjusted = np.empty(n, dtype=float)

    cumulative_min = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        value = ranked[idx] * n / rank
        cumulative_min = min(cumulative_min, value)
        adjusted[idx] = cumulative_min

    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def safe_log2_fold_change(case_mean: np.ndarray, control_mean: np.ndarray) -> np.ndarray:
    """Compute a stable log2 fold change with a unit pseudocount."""
    return np.log2((np.asarray(case_mean, dtype=float) + 1.0) / (np.asarray(control_mean, dtype=float) + 1.0))


def merge_uns_dict(target: dict, source: dict | None) -> dict:
    """Recursively merge two dictionaries."""
    if source is None:
        return target
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            merge_uns_dict(target[key], value)
        else:
            target[key] = value
    return target
