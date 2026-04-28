"""Schema helpers for AnnData objects consumed by SpatialPerturb."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from ._utils import merge_uns_dict

REQUIRED_OBS_COLUMNS = ("perturbation", "perturbation_status", "cell_type", "roi")
DEFAULT_OBS_VALUES = {
    "perturbation": "unassigned",
    "perturbation_status": "unassigned",
    "cell_type": "unknown",
    "roi": "global",
}
SPATIAL_KEY = "spatial"
GRAPH_KEYS = ("sp_knn", "sp_radius")
SPATIALPERTURB_UNS_KEY = "spatialperturb"
SCHEMA_VERSION = "0.3.0"


class SpatialPerturbSchemaError(ValueError):
    """Raised when an AnnData object does not satisfy the SpatialPerturb schema."""


def _default_spatial_coords(n_obs: int) -> np.ndarray:
    x = np.arange(n_obs, dtype=float)
    y = np.zeros(n_obs, dtype=float)
    return np.column_stack([x, y])


def ensure_spatialperturb_schema(
    adata: AnnData,
    *,
    metadata: dict[str, Any] | None = None,
    copy: bool = False,
) -> AnnData:
    """Ensure required schema fields exist, filling defaults when needed."""
    target = adata.copy() if copy else adata

    for column, default_value in DEFAULT_OBS_VALUES.items():
        if column not in target.obs:
            target.obs[column] = default_value
        target.obs[column] = target.obs[column].astype(str)

    if SPATIAL_KEY not in target.obsm:
        target.obsm[SPATIAL_KEY] = _default_spatial_coords(target.n_obs)

    spatial = np.asarray(target.obsm[SPATIAL_KEY], dtype=float)
    if spatial.ndim != 2 or spatial.shape[1] != 2:
        raise SpatialPerturbSchemaError(
            f"obsm['{SPATIAL_KEY}'] must be an (n_obs, 2) array, found shape {spatial.shape}."
        )

    if SPATIALPERTURB_UNS_KEY not in target.uns or not isinstance(target.uns[SPATIALPERTURB_UNS_KEY], dict):
        target.uns[SPATIALPERTURB_UNS_KEY] = {}

    merge_uns_dict(
        target.uns[SPATIALPERTURB_UNS_KEY],
        {
            "schema_version": SCHEMA_VERSION,
            "graph_keys": [key for key in GRAPH_KEYS if key in target.obsp],
        },
    )
    merge_uns_dict(target.uns[SPATIALPERTURB_UNS_KEY], metadata)
    return target


def validate_spatialperturb_schema(
    adata: AnnData,
    *,
    require_graph: bool = False,
) -> None:
    """Validate that an AnnData object follows the SpatialPerturb schema."""
    missing_obs = [column for column in REQUIRED_OBS_COLUMNS if column not in adata.obs]
    if missing_obs:
        raise SpatialPerturbSchemaError(f"Missing required obs columns: {missing_obs}")

    if SPATIAL_KEY not in adata.obsm:
        raise SpatialPerturbSchemaError("Missing required obsm['spatial'] coordinates.")

    spatial = np.asarray(adata.obsm[SPATIAL_KEY])
    if spatial.ndim != 2 or spatial.shape[0] != adata.n_obs or spatial.shape[1] != 2:
        raise SpatialPerturbSchemaError(
            f"obsm['{SPATIAL_KEY}'] must have shape (n_obs, 2); found {spatial.shape}."
        )

    if SPATIALPERTURB_UNS_KEY not in adata.uns or not isinstance(adata.uns[SPATIALPERTURB_UNS_KEY], dict):
        raise SpatialPerturbSchemaError("Missing required uns['spatialperturb'] metadata dictionary.")

    if require_graph and not any(key in adata.obsp for key in GRAPH_KEYS):
        raise SpatialPerturbSchemaError(
            "A spatial graph is required but neither obsp['sp_knn'] nor obsp['sp_radius'] is available."
        )


def schema_summary(adata: AnnData) -> pd.Series:
    """Return a concise schema summary for validation reports."""
    validate_spatialperturb_schema(adata, require_graph=False)
    return pd.Series(
        {
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "has_sp_knn": "sp_knn" in adata.obsp,
            "has_sp_radius": "sp_radius" in adata.obsp,
            "perturbations": int(adata.obs["perturbation"].nunique()),
            "cell_types": int(adata.obs["cell_type"].nunique()),
            "rois": int(adata.obs["roi"].nunique()),
        }
    )
