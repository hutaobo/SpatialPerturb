"""AnnData I/O helpers for SpatialPerturb."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
from anndata import AnnData
import numpy as np
import pandas as pd
from scipy import sparse

from .schema import ensure_spatialperturb_schema, validate_spatialperturb_schema


def _as_dataframe(
    table: pd.DataFrame | None,
    *,
    index: pd.Index,
    default_prefix: str,
) -> pd.DataFrame:
    if table is None:
        return pd.DataFrame(index=index)
    target = table.copy()
    target.index = target.index.astype(str)
    if not target.index.equals(index):
        target = target.reindex(index)
    if target.index.isnull().any():
        target.index = pd.Index([f"{default_prefix}_{idx}" for idx in range(len(target))], dtype="object")
    return target


def _infer_spatial(obs: pd.DataFrame, spatial: Any | None) -> np.ndarray | None:
    if spatial is not None:
        coords = np.asarray(spatial, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("spatial coordinates must have shape (n_obs, 2).")
        return coords

    if {"x", "y"}.issubset(obs.columns):
        return obs.loc[:, ["x", "y"]].to_numpy(dtype=float)
    if {"X", "Y"}.issubset(obs.columns):
        return obs.loc[:, ["X", "Y"]].to_numpy(dtype=float)
    if {"x_centroid", "y_centroid"}.issubset(obs.columns):
        return obs.loc[:, ["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    return None


def from_tables(
    expression: pd.DataFrame | np.ndarray | sparse.spmatrix,
    *,
    obs: pd.DataFrame | None = None,
    var: pd.DataFrame | None = None,
    spatial: pd.DataFrame | np.ndarray | None = None,
    layers: dict[str, pd.DataFrame | np.ndarray | sparse.spmatrix] | None = None,
    metadata: dict[str, Any] | None = None,
) -> AnnData:
    """Build a schema-compliant AnnData object from tabular inputs."""
    if isinstance(expression, pd.DataFrame):
        obs_index = pd.Index(expression.index.astype(str), dtype="object")
        var_index = pd.Index(expression.columns.astype(str), dtype="object")
        matrix = expression.to_numpy(dtype=float)
    else:
        matrix = expression
        if obs is None or var is None:
            raise ValueError("obs and var must be provided when expression is not a DataFrame.")
        obs_index = pd.Index(obs.index.astype(str), dtype="object")
        var_index = pd.Index(var.index.astype(str), dtype="object")

    obs_frame = _as_dataframe(obs, index=obs_index, default_prefix="cell")
    var_frame = _as_dataframe(var, index=var_index, default_prefix="gene")

    if sparse.issparse(matrix):
        adata = ad.AnnData(X=matrix.tocsr(), obs=obs_frame, var=var_frame)
    else:
        adata = ad.AnnData(X=np.asarray(matrix), obs=obs_frame, var=var_frame)

    if layers:
        for layer_name, layer_values in layers.items():
            adata.layers[str(layer_name)] = layer_values.to_numpy() if isinstance(layer_values, pd.DataFrame) else layer_values

    coords = _infer_spatial(obs_frame, spatial)
    if coords is not None:
        adata.obsm["spatial"] = coords

    ensure_spatialperturb_schema(adata, metadata=metadata or {})
    validate_spatialperturb_schema(adata)
    return adata


def _read_expression_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    sep = "\t" if ".tsv" in suffixes else ","
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_csv(path, sep=sep, index_col=0, compression=compression)


def _read_cells_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    sep = "\t" if ".tsv" in suffixes else ","
    compression = "gzip" if path.suffix == ".gz" else None
    table = pd.read_csv(path, sep=sep, index_col=0, compression=compression)
    table.index = table.index.astype(str)
    return table


def _read_platform_directory(path: str | Path, *, platform: str) -> AnnData:
    directory = Path(path).expanduser().resolve()
    if directory.is_file():
        if directory.suffix == ".h5ad":
            adata = ad.read_h5ad(directory)
            ensure_spatialperturb_schema(adata, metadata={"platform": platform, "source_path": str(directory)})
            validate_spatialperturb_schema(adata)
            return adata
        raise ValueError(f"Unsupported file input for {platform}: {directory}")

    if not directory.exists():
        raise FileNotFoundError(directory)

    counts_candidates = [
        directory / "counts.csv",
        directory / "counts.tsv",
        directory / "counts.csv.gz",
        directory / "counts.tsv.gz",
        directory / "matrix.csv",
        directory / "matrix.tsv",
        directory / "matrix.csv.gz",
        directory / "matrix.tsv.gz",
    ]
    cells_candidates = [
        directory / "cells.csv",
        directory / "cells.tsv",
        directory / "cells.csv.gz",
        directory / "cells.tsv.gz",
        directory / "metadata.csv",
        directory / "metadata.tsv",
        directory / "metadata.csv.gz",
        directory / "metadata.tsv.gz",
    ]

    counts_path = next((candidate for candidate in counts_candidates if candidate.exists()), None)
    cells_path = next((candidate for candidate in cells_candidates if candidate.exists()), None)
    if counts_path is None or cells_path is None:
        raise FileNotFoundError(
            f"Could not locate paired counts/cells tables in {directory}. Expected files like counts.csv and cells.csv."
        )

    expression = _read_expression_table(counts_path)
    cells = _read_cells_table(cells_path).reindex(expression.index.astype(str))
    adata = from_tables(expression, obs=cells, metadata={"platform": platform, "source_path": str(directory)})

    molecules_candidates = [
        directory / "molecules.csv",
        directory / "molecules.tsv",
        directory / "molecules.csv.gz",
        directory / "molecules.tsv.gz",
    ]
    molecules_path = next((candidate for candidate in molecules_candidates if candidate.exists()), None)
    if molecules_path is not None:
        adata.uns.setdefault("spatialperturb", {})
        adata.uns["spatialperturb"]["molecules"] = {"path": str(molecules_path)}
    return adata


def _annotate_cell_groups(adata: AnnData, path: Path) -> None:
    table = pd.read_csv(path)
    if "cell_id" not in table.columns or "group" not in table.columns:
        raise ValueError("cell_group_path must contain at least 'cell_id' and 'group' columns.")
    mapping = dict(zip(table["cell_id"].astype(str), table["group"].astype(str), strict=False))
    matched = 0
    cell_types = []
    for cell_id, current in zip(adata.obs_names.astype(str), adata.obs["cell_type"].astype(str), strict=False):
        group = mapping.get(cell_id)
        if group is None:
            cell_types.append(current)
        else:
            cell_types.append(group)
            matched += 1
    adata.obs["cell_type"] = cell_types
    adata.uns.setdefault("spatialperturb", {})
    adata.uns["spatialperturb"]["cell_group_annotation"] = {"path": str(path), "matched_cells": matched}


def _point_in_polygon(x: float, y: float, polygon: list[list[float]]) -> bool:
    inside = False
    j = len(polygon) - 1
    for i in range(len(polygon)):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        intersects = ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / max((yj - yi), 1e-12) + xi)
        if intersects:
            inside = not inside
        j = i
    return inside


def _annotate_roi_geojson(adata: AnnData, path: Path) -> None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    features = payload.get("features", [])
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    rois = adata.obs["roi"].astype(str).tolist()
    matched = 0
    for idx, (x_coord, y_coord) in enumerate(coords):
        for feature in features:
            geometry = feature.get("geometry", {})
            if geometry.get("type") != "Polygon":
                continue
            rings = geometry.get("coordinates", [])
            if not rings:
                continue
            if _point_in_polygon(float(x_coord), float(y_coord), rings[0]):
                rois[idx] = str(feature.get("properties", {}).get("assigned_structure", rois[idx]))
                matched += 1
                break
    adata.obs["roi"] = rois
    adata.uns.setdefault("spatialperturb", {})
    adata.uns["spatialperturb"]["roi_annotation"] = {"path": str(path), "matched_cells": matched}


def read_xenium(
    path: str | Path,
    *,
    cell_group_path: str | Path | None = None,
    roi_geojson_path: str | Path | None = None,
    sample_name: str | None = None,
    load_molecules: bool = False,
) -> AnnData:
    """Read a simple Xenium-style cell-level export into AnnData."""
    adata = _read_platform_directory(path, platform="xenium")
    if sample_name is not None:
        adata.obs["sample"] = str(sample_name)
        adata.uns.setdefault("spatialperturb", {})
        adata.uns["spatialperturb"]["sample_name"] = str(sample_name)
    if cell_group_path is not None:
        _annotate_cell_groups(adata, Path(cell_group_path))
    if roi_geojson_path is not None:
        _annotate_roi_geojson(adata, Path(roi_geojson_path))
    if not load_molecules:
        adata.uns.get("spatialperturb", {}).pop("molecules", None)
    return adata


def read_stereoseq(path: str | Path, **_: Any) -> AnnData:
    """Read a simple Stereo-seq-style cell-level export into AnnData."""
    return _read_platform_directory(path, platform="stereoseq")
