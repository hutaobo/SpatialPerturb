"""Preprocessing and perturbation assignment utilities."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from ._utils import extract_matrix, safe_log2_fold_change
from .schema import ensure_spatialperturb_schema


def assign_perturbations(
    adata: AnnData,
    *,
    barcode_columns: Sequence[str] | None = None,
    barcode_prefix: str | None = None,
    barcode_to_perturbation: Mapping[str, str] | None = None,
    min_counts: float = 1.0,
    negative_label: str = "unassigned",
    multiple_label: str = "multiple",
    perturbation_col: str = "perturbation",
    status_col: str = "perturbation_status",
    include_col: str = "include_for_inference",
    copy: bool = False,
) -> AnnData:
    """Assign a perturbation label to each cell from barcode feature counts."""
    target = adata.copy() if copy else adata
    ensure_spatialperturb_schema(target)

    if barcode_columns is None:
        if barcode_prefix is None:
            barcode_columns = target.uns.get("spatialperturb", {}).get("barcode_columns")
        else:
            barcode_columns = [name for name in target.var_names if str(name).startswith(barcode_prefix)]

    if not barcode_columns:
        raise ValueError("No barcode features were provided or discovered.")

    barcode_columns = [str(column) for column in barcode_columns]
    barcode_matrix = extract_matrix(target, var_names=barcode_columns)
    detected = barcode_matrix >= float(min_counts)
    positive_counts = detected.sum(axis=1)
    top_index = barcode_matrix.argmax(axis=1)
    top_counts = barcode_matrix.max(axis=1)

    perturbations = np.full(target.n_obs, negative_label, dtype=object)
    statuses = np.full(target.n_obs, "unassigned", dtype=object)

    single_mask = positive_counts == 1
    multiple_mask = positive_counts > 1

    called_barcodes = np.array(barcode_columns, dtype=object)[top_index]
    if barcode_to_perturbation is None:
        barcode_to_perturbation = {}

    perturbations[single_mask] = [
        barcode_to_perturbation.get(barcode, barcode) for barcode in called_barcodes[single_mask]
    ]
    statuses[single_mask] = "single"
    perturbations[multiple_mask] = multiple_label
    statuses[multiple_mask] = "multiple"

    target.obs[perturbation_col] = perturbations
    target.obs[status_col] = statuses
    target.obs[include_col] = statuses == "single"
    target.obs["barcode_positive_count"] = positive_counts.astype(int)
    target.obs["barcode_max_count"] = top_counts.astype(float)

    target.uns["spatialperturb"]["barcode_columns"] = barcode_columns
    target.uns["spatialperturb"]["barcode_min_counts"] = float(min_counts)
    target.uns["spatialperturb"]["assignment_labels"] = {
        "negative": negative_label,
        "multiple": multiple_label,
    }
    return target


def qc_perturbations(
    adata: AnnData,
    *,
    control: str,
    target_map: Mapping[str, str] | None = None,
    perturbation_col: str = "perturbation",
    status_col: str = "perturbation_status",
    layer: str | None = None,
    min_cells: int = 5,
) -> pd.DataFrame:
    """Summarize perturbation calling quality and target knockdown sanity checks."""
    ensure_spatialperturb_schema(adata)
    target_map = dict(target_map or {})

    status_series = adata.obs[status_col].astype(str)
    perturbation_series = adata.obs[perturbation_col].astype(str)
    single_mask = status_series == "single"
    control_mask = single_mask & (perturbation_series == str(control))

    records: list[dict] = []
    for perturbation in sorted(perturbation_series.unique()):
        perturbation_mask = single_mask & (perturbation_series == perturbation)
        n_cells = int(perturbation_mask.sum())
        target_gene = target_map.get(perturbation)
        mean_case = np.nan
        mean_control = np.nan
        log2fc = np.nan
        if target_gene is not None and target_gene in adata.var_names and n_cells > 0 and control_mask.sum() > 0:
            case_expr = extract_matrix(adata, layer=layer, obs_names=adata.obs_names[perturbation_mask], var_names=[target_gene]).ravel()
            control_expr = extract_matrix(
                adata,
                layer=layer,
                obs_names=adata.obs_names[control_mask],
                var_names=[target_gene],
            ).ravel()
            mean_case = float(case_expr.mean())
            mean_control = float(control_expr.mean())
            log2fc = float(safe_log2_fold_change(np.array([mean_case]), np.array([mean_control]))[0])

        records.append(
            {
                "perturbation": perturbation,
                "n_cells": n_cells,
                "fraction_cells": n_cells / max(int(single_mask.sum()), 1),
                "target_gene": target_gene,
                "target_mean_case": mean_case,
                "target_mean_control": mean_control,
                "target_log2fc": log2fc,
                "valid_for_inference": bool(
                    perturbation not in {"unassigned", "multiple"} and n_cells >= min_cells
                ),
            }
        )

    return pd.DataFrame.from_records(records).sort_values(["valid_for_inference", "n_cells"], ascending=[False, False])
