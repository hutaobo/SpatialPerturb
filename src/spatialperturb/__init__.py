"""AnnData-native tools for spatial perturbation analysis."""

from __future__ import annotations

import importlib

__version__ = "0.3.0"

from .benchmarks import (
    available_benchmarks,
    run_core_benchmark,
    run_cross_platform_benchmark,
    run_nature_methods_breast_analysis,
    run_reference_projection_benchmark,
)
from .datasets import (
    available_datasets,
    fetch_dataset,
    get_dataset_card,
    load_demo_dataset,
    load_public_dataset,
    prepare_dataset,
)
from .gr import build_spatial_graph, collect_neighbors
from .io import from_tables, read_stereoseq, read_xenium
from .pp import assign_perturbations, qc_perturbations
from .reports import render_paper_figures
from .schema import (
    DEFAULT_OBS_VALUES,
    REQUIRED_OBS_COLUMNS,
    SpatialPerturbSchemaError,
    ensure_spatialperturb_schema,
    validate_spatialperturb_schema,
)
from .signatures import (
    aggregate_program_scores,
    bootstrap_program_score_intervals,
    build_reference_programs,
    build_signature_matrix,
    calibrate_program_scores,
    compare_program_concordance,
    derive_perturbation_programs,
    neighbor_program_scores,
    program_redundancy_table,
    score_programs,
    spatial_autocorrelation_scores,
    validate_reference_programs,
)
from .tl import differential_lr, intrinsic_de, neighbor_de, platform_concordance, power_curve

benchmarks = importlib.import_module(".benchmarks", __name__)
datasets = importlib.import_module(".datasets", __name__)
gr = importlib.import_module(".gr", __name__)
io = importlib.import_module(".io", __name__)
pl = importlib.import_module(".pl", __name__)
pp = importlib.import_module(".pp", __name__)
reports = importlib.import_module(".reports", __name__)
schema = importlib.import_module(".schema", __name__)
signatures = importlib.import_module(".signatures", __name__)
tl = importlib.import_module(".tl", __name__)

__all__ = [
    "__version__",
    "DEFAULT_OBS_VALUES",
    "REQUIRED_OBS_COLUMNS",
    "SpatialPerturbSchemaError",
    "aggregate_program_scores",
    "assign_perturbations",
    "available_benchmarks",
    "available_datasets",
    "benchmarks",
    "bootstrap_program_score_intervals",
    "build_reference_programs",
    "build_signature_matrix",
    "build_spatial_graph",
    "calibrate_program_scores",
    "collect_neighbors",
    "compare_program_concordance",
    "datasets",
    "derive_perturbation_programs",
    "differential_lr",
    "ensure_spatialperturb_schema",
    "fetch_dataset",
    "from_tables",
    "get_dataset_card",
    "gr",
    "intrinsic_de",
    "io",
    "load_demo_dataset",
    "load_public_dataset",
    "neighbor_de",
    "neighbor_program_scores",
    "platform_concordance",
    "pl",
    "power_curve",
    "pp",
    "prepare_dataset",
    "program_redundancy_table",
    "qc_perturbations",
    "read_stereoseq",
    "read_xenium",
    "render_paper_figures",
    "reports",
    "run_core_benchmark",
    "run_cross_platform_benchmark",
    "run_nature_methods_breast_analysis",
    "run_reference_projection_benchmark",
    "schema",
    "score_programs",
    "signatures",
    "spatial_autocorrelation_scores",
    "tl",
    "validate_reference_programs",
    "validate_spatialperturb_schema",
]
