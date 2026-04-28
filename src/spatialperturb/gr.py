"""Spatial graph construction and neighborhood collection."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from .schema import ensure_spatialperturb_schema


def build_spatial_graph(
    adata: AnnData,
    *,
    mode: str = "knn",
    k: int = 15,
    radius: float | None = None,
    spatial_key: str = "spatial",
    graph_key: str | None = None,
    include_self: bool = False,
    copy: bool = False,
) -> AnnData:
    """Construct a spatial neighbor graph and store it in ``adata.obsp``."""
    target = adata.copy() if copy else adata
    ensure_spatialperturb_schema(target)
    coords = np.asarray(target.obsm[spatial_key], dtype=float)

    if mode not in {"knn", "radius"}:
        raise ValueError("mode must be 'knn' or 'radius'.")

    if mode == "knn":
        if k < 1:
            raise ValueError("k must be at least 1 for knn graphs.")
        key = graph_key or "sp_knn"
        neighbors = min(target.n_obs, k + 1)
        model = NearestNeighbors(n_neighbors=neighbors, metric="euclidean")
        model.fit(coords)
        distances, indices = model.kneighbors(coords)
        if not include_self:
            distances = distances[:, 1:]
            indices = indices[:, 1:]

        rows = np.repeat(np.arange(target.n_obs), indices.shape[1])
        graph = sparse.csr_matrix((distances.ravel(), (rows, indices.ravel())), shape=(target.n_obs, target.n_obs))
    else:
        if radius is None or radius <= 0:
            raise ValueError("A positive radius is required when mode='radius'.")
        key = graph_key or "sp_radius"
        model = NearestNeighbors(radius=radius, metric="euclidean")
        model.fit(coords)
        graph = model.radius_neighbors_graph(coords, radius=radius, mode="distance")
        if not include_self:
            graph = graph.tolil()
            graph.setdiag(0)
            graph = graph.tocsr()
            graph.eliminate_zeros()

    graph = graph.maximum(graph.T)
    target.obsp[key] = graph
    target.uns["spatialperturb"]["graph_keys"] = sorted({*target.uns["spatialperturb"].get("graph_keys", []), key})
    target.uns["spatialperturb"]["graph_params"] = {
        "mode": mode,
        "k": int(k),
        "radius": None if radius is None else float(radius),
        "include_self": include_self,
    }
    return target


def collect_neighbors(
    adata: AnnData,
    *,
    cells: Sequence[str] | None = None,
    graph_key: str | None = None,
    exclude_perturbed: bool = True,
    status_col: str = "perturbation_status",
    include_self: bool = False,
    as_frame: bool = False,
) -> dict[str, list[str]] | pd.DataFrame:
    """Collect neighbors for a set of cells from a stored spatial graph."""
    ensure_spatialperturb_schema(adata)
    key = graph_key or ("sp_knn" if "sp_knn" in adata.obsp else "sp_radius")
    if key not in adata.obsp:
        raise KeyError(f"Spatial graph {key!r} not found in adata.obsp.")

    graph = adata.obsp[key].tocsr()
    selected_cells = list(adata.obs_names if cells is None else map(str, cells))
    allowed_mask = np.ones(adata.n_obs, dtype=bool)
    if exclude_perturbed:
        allowed_mask = adata.obs[status_col].astype(str).to_numpy() == "unassigned"

    mapping: dict[str, list[str]] = {}
    records: list[dict] = []
    for cell in selected_cells:
        idx = int(adata.obs_names.get_loc(cell))
        row = graph.getrow(idx)
        neighbors = []
        for neighbor_idx, distance in zip(row.indices, row.data, strict=False):
            if not include_self and neighbor_idx == idx:
                continue
            if exclude_perturbed and not allowed_mask[neighbor_idx]:
                continue
            neighbor_name = str(adata.obs_names[neighbor_idx])
            neighbors.append(neighbor_name)
            records.append(
                {
                    "source": str(cell),
                    "neighbor": neighbor_name,
                    "distance": float(distance),
                    "graph_key": key,
                }
            )
        mapping[str(cell)] = neighbors

    if as_frame:
        return pd.DataFrame.from_records(records)
    return mapping
