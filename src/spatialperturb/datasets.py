"""Dataset registry, lifecycle helpers, and demo data for SpatialPerturb."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import gzip
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
from typing import Any
from urllib.request import urlopen

import anndata as ad
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io as spio
from scipy import sparse

from .gr import build_spatial_graph
from .io import from_tables
from .pp import assign_perturbations
from .schema import ensure_spatialperturb_schema, validate_spatialperturb_schema

_GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
_GSE274047_TARGETS = (
    "mSafe",
    "Trem2",
    "Rraga",
    "Myrf",
    "Fasn",
    "Clu",
    "Dpp6",
    "Tbk1",
    "Flcn",
    "Gfap",
    "C9orf72",
    "Cfap410",
    "Stk39",
    "Lrrk2",
    "Ndufaf",
    "Sh3gl2",
    "Srf",
    "Rbfox",
    "Olig2",
)
_GRNA_ALIAS_MAP = {
    "mSafe": "control",
    "Ndufaf": "Ndufaf2",
    "Rbfox": "Rbfox3",
}
_CONTROL_GUIDE_RE = re.compile(r"(?:^|[_\-\s])(NT|non.?target|intergenic|control)(?:$|[_\-\s])", re.IGNORECASE)
_GUIDE_SUFFIX_RE = re.compile(r"(?:_sgRNA\d+|_sg\d+)$", re.IGNORECASE)
_FLAT_10X_FILE_RE = re.compile(
    r"^(?P<gsm>GSM\d+)_(?P<kind>barcodes|features|matrix|protospacer_calls_per_cell|protospacer_umi_thresholds)_(?P<sample>.+?)\.(?:tsv|csv|mtx)\.gz$",
    re.IGNORECASE,
)
_CONTROL_COUNT_RE = re.compile(r"^(?P<gsm>GSM\d+)_counts_(?P<sample>.+?)\.csv\.gz$", re.IGNORECASE)
_KNOWN_CELL_LINES = ("A549", "BxPC3", "HAP1", "HCC1143", "HCC2157", "HCC38", "HT29", "K562", "MCF7")


@dataclass(frozen=True)
class DatasetCard:
    """Metadata for a benchmarkable dataset."""

    name: str
    title: str
    platform: str
    accession: str
    benchmark_track: str
    source_url: str
    description: str
    status: str = "metadata_only"
    raw_format: str = "none"
    parser: str = "demo"
    prepared_kind: str = "adata"
    raw_archive_name: str | None = None
    prepared_filename: str | None = None
    matrix_filename: str | None = None
    supports_auto_prepare: bool = True
    validation_note: str | None = None
    supplementary_files: tuple[str, ...] = ()


def _geo_series_bucket(accession: str) -> str:
    prefix = accession[:6]
    return f"{prefix}nnn"


def _geo_raw_url(accession: str) -> str:
    bucket = _geo_series_bucket(accession)
    return f"{_GEO_BASE}/{bucket}/{accession}/suppl/{accession}_RAW.tar"


def _geo_matrix_url(accession: str) -> str:
    bucket = _geo_series_bucket(accession)
    return f"{_GEO_BASE}/{bucket}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def _geo_supplementary_url(accession: str, filename: str) -> str:
    bucket = _geo_series_bucket(accession)
    return f"{_GEO_BASE}/{bucket}/{accession}/suppl/{filename}"


_DATASET_CATALOG = [
    DatasetCard(
        name="shen_2026_stereoseq",
        title="Spatial Perturb-Seq Stereo-seq hippocampus",
        platform="stereoseq",
        accession="GSE274447",
        benchmark_track="shen_2026_core",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274447",
        description="Primary spatial perturbation benchmark from Shen et al. 2026.",
        status="raw_downloadable",
        raw_format="tar_of_gef",
        parser="geo_gef_passthrough",
        prepared_kind="adata",
        raw_archive_name="GSE274447_RAW.tar",
        prepared_filename="shen_2026_stereoseq.h5ad",
        matrix_filename="GSE274447_series_matrix.txt.gz",
        supports_auto_prepare=False,
        validation_note="Automatic fetch and extraction are supported. Final preparation requires either a preconverted h5ad/tables or an external GEF converter.",
    ),
    DatasetCard(
        name="shen_2026_scrnaseq",
        title="Perturb-seq scRNA-seq companion dataset",
        platform="scrnaseq",
        accession="GSE274058",
        benchmark_track="cross_platform_concordance",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274058",
        description="Reference dissociated Perturb-seq benchmark from the same study.",
        status="raw_downloadable",
        raw_format="tar_of_10x_tar_gz",
        parser="geo_nested_10x_tar",
        prepared_kind="adata",
        raw_archive_name="GSE274058_RAW.tar",
        prepared_filename="shen_2026_scrnaseq.h5ad",
        matrix_filename="GSE274058_series_matrix.txt.gz",
        supports_auto_prepare=True,
        validation_note="Automatic fetch, extraction, and AnnData preparation are supported from the GEO raw archive.",
    ),
    DatasetCard(
        name="gse241115_breast_cropseq",
        title="Breast cancer CROP-seq master regulator screen",
        platform="cropseq",
        accession="GSE241115",
        benchmark_track="breast_reference_projection",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241115",
        description="Breast cancer stem-like state CROP-seq reference dataset with HCC38 and HCC1143 perturbations plus control profiles.",
        status="raw_downloadable",
        raw_format="tar_of_csv_mtx_tsv",
        parser="geo_flat_10x_tar",
        prepared_kind="adata",
        raw_archive_name="GSE241115_RAW.tar",
        prepared_filename="gse241115_breast_cropseq.h5ad",
        supports_auto_prepare=True,
        validation_note="Automatic preparation parses control count CSVs, sparse expression matrices, and protospacer calls from the GEO tarball.",
    ),
    DatasetCard(
        name="gse281048_pathway_atlas",
        title="Perturb-seq signaling pathway atlas",
        platform="perturbseq_pathway_atlas",
        accession="GSE281048",
        benchmark_track="breast_reference_projection",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281048",
        description="Multi-cell-line pathway Perturb-seq atlas containing MCF7 responses across IFNB, IFNG, TGFB, TNFA, and INS.",
        status="raw_downloadable",
        raw_format="gzipped_seurat_rds",
        parser="geo_rds_collection",
        prepared_kind="adata",
        prepared_filename="gse281048_pathway_atlas.h5ad",
        supports_auto_prepare=True,
        validation_note="Automatic preparation uses an R-backed converter that exports Seurat objects to sparse tables before Python standardization.",
        supplementary_files=(
            "GSE281048_Seurat_object_IFNB_Perturb_seq.rds.gz",
            "GSE281048_Seurat_object_IFNG_Perturb_seq.rds.gz",
            "GSE281048_Seurat_object_INS_Perturb_seq.rds.gz",
            "GSE281048_Seurat_object_TGFB_Perturb_seq.rds.gz",
            "GSE281048_Seurat_object_TNFA_Perturb_seq.rds.gz",
        ),
    ),
    DatasetCard(
        name="demo_spatialperturb",
        title="Synthetic paired demo dataset",
        platform="paired_demo",
        accession="synthetic",
        benchmark_track="demo_end_to_end",
        source_url="https://github.com/hutaobo/SpatialPerturb",
        description="Deterministic demo data covering intrinsic, neighbor, and concordance workflows.",
        status="built_in",
        raw_format="synthetic",
        parser="demo",
        prepared_kind="adata",
        prepared_filename="demo_spatialperturb.h5ad",
        supports_auto_prepare=True,
    ),
]


def _dataset_layout(name: str, cache_dir: str | Path) -> dict[str, Path]:
    root = Path(cache_dir).expanduser().resolve() / name
    raw_dir = root / "raw"
    prepared_dir = root / "prepared"
    reports_dir = root / "reports"
    return {
        "root": root,
        "raw_dir": raw_dir,
        "prepared_dir": prepared_dir,
        "reports_dir": reports_dir,
        "metadata_path": root / "dataset.json",
    }


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=_json_default), encoding="utf-8")


def _download_file(url: str, destination: Path, *, force: bool = False) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        return destination

    with urlopen(url) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return destination


def _extract_tar_once(archive_path: Path, destination: Path) -> list[Path]:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path) as archive:
        extracted = [destination / member.name for member in archive.getmembers() if member.isfile()]
        if extracted and all(path.exists() for path in extracted):
            return extracted
        try:
            archive.extractall(destination, filter="data")
        except TypeError:
            archive.extractall(destination)
    return extracted


def _guess_table_pair(raw_dir: Path) -> tuple[Path, Path] | None:
    count_candidates = [
        *raw_dir.glob("*.h5ad"),
        *raw_dir.glob("*.csv"),
        *raw_dir.glob("*.csv.gz"),
        *raw_dir.glob("*.tsv"),
        *raw_dir.glob("*.tsv.gz"),
    ]
    h5ad_candidates = [path for path in count_candidates if path.suffix == ".h5ad"]
    if len(h5ad_candidates) == 1:
        return h5ad_candidates[0], h5ad_candidates[0]

    counts = [path for path in count_candidates if "count" in path.name.lower() or "expr" in path.name.lower()]
    cells = [path for path in count_candidates if "cell" in path.name.lower() or "meta" in path.name.lower()]
    if counts and cells:
        return counts[0], cells[0]
    return None


def _infer_gRNA_map(feature_names: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for feature in feature_names:
        if feature.endswith("_gRNA"):
            perturbation = feature.removesuffix("_gRNA")
            perturbation = _GRNA_ALIAS_MAP.get(perturbation, perturbation)
            mapping[feature] = perturbation
    return mapping


def _extract_nested_10x_archives(archive_path: Path, extracted_root: Path) -> list[Path]:
    extracted_root.mkdir(parents=True, exist_ok=True)
    sample_dirs = [path for path in extracted_root.iterdir() if path.is_dir()]
    if sample_dirs:
        return sorted(sample_dirs)

    with tarfile.open(archive_path) as outer:
        for member in outer.getmembers():
            if not member.isfile() or not member.name.endswith(".tar.gz"):
                continue
            sample_name = member.name.replace("_filtered_feature_bc_matrix.tar.gz", "")
            sample_root = extracted_root / sample_name
            sample_root.mkdir(parents=True, exist_ok=True)
            inner_stream = outer.extractfile(member)
            if inner_stream is None:
                continue
            with tarfile.open(fileobj=inner_stream, mode="r:gz") as inner:
                try:
                    inner.extractall(sample_root, filter="data")
                except TypeError:
                    inner.extractall(sample_root)
    return sorted([path for path in extracted_root.iterdir() if path.is_dir()])


def _sample_matrix_dir(sample_root: Path) -> Path:
    if (sample_root / "filtered_feature_bc_matrix").exists():
        return sample_root / "filtered_feature_bc_matrix"
    return sample_root


def _normalize_obs_index(values: pd.Index | list[str] | np.ndarray) -> pd.Index:
    index = pd.Index(values).astype(str)
    series = pd.Series(index.to_list(), index=index, dtype="string")
    return pd.Index(series.str.replace(r"\.1$", "-1", regex=True).tolist(), dtype="object")


def _is_control_perturbation(name: str) -> bool:
    return bool(_CONTROL_GUIDE_RE.search(str(name)))


def _normalize_perturbation_name(guide_id: str | None, *, control_label: str = "control") -> str:
    if guide_id is None:
        return "unassigned"
    guide = str(guide_id).strip()
    if not guide:
        return "unassigned"
    base = _GUIDE_SUFFIX_RE.sub("", guide)
    if _is_control_perturbation(base):
        return control_label
    return base


def _finalize_perturbation_obs(
    obs: pd.DataFrame,
    *,
    guide_column: str = "guide_id",
    feature_count_column: str | None = "num_features",
    control_label: str = "control",
) -> pd.DataFrame:
    target = obs.copy()
    guide_values = target.get(guide_column, pd.Series(index=target.index, dtype="object")).astype("string")
    feature_counts = (
        pd.to_numeric(target[feature_count_column], errors="coerce")
        if feature_count_column is not None and feature_count_column in target.columns
        else pd.Series(1.0, index=target.index, dtype=float)
    )

    perturbations: list[str] = []
    statuses: list[str] = []
    for cell_id in target.index.astype(str):
        guide_id = guide_values.get(cell_id)
        if pd.isna(guide_id) or str(guide_id).strip() == "":
            perturbations.append("unassigned")
            statuses.append("unassigned")
            continue
        guide_id = str(guide_id)
        multiple = bool(feature_counts.get(cell_id, 1.0) > 1)
        if not multiple and any(token in guide_id for token in ("|", ";", ",")):
            multiple = True
        if multiple:
            perturbations.append("multiple")
            statuses.append("multiple")
            continue
        perturbations.append(_normalize_perturbation_name(guide_id, control_label=control_label))
        statuses.append("single")

    target["guide_id"] = guide_values.astype("object")
    target["perturbation"] = perturbations
    target["perturbation_status"] = statuses
    target["include_for_inference"] = target["perturbation_status"].astype(str) == "single"
    target.loc[target["perturbation"].astype(str) == control_label, "include_for_inference"] = True
    return target


def _read_control_counts_table(path: Path, *, sample_name: str, cell_line: str, dataset_name: str) -> AnnData:
    expression = pd.read_csv(path, index_col=0, compression="gzip").transpose()
    expression.index = _normalize_obs_index(expression.index)
    expression.columns = expression.columns.astype(str)

    obs = pd.DataFrame(
        {
            "sample": sample_name,
            "cell_line": cell_line,
            "guide_id": "control",
            "perturbation": "control",
            "perturbation_status": "single",
            "include_for_inference": True,
            "reference_source": dataset_name,
            "cell_type": "unknown",
            "roi": "global",
        },
        index=expression.index.astype(str),
    )
    var = pd.DataFrame(index=expression.columns.astype(str))
    var.index.name = None
    matrix = sparse.csr_matrix(expression.to_numpy(dtype=float))
    adata = from_tables(
        matrix,
        obs=obs,
        var=var,
        metadata={
            "dataset_name": dataset_name,
            "platform": "cropseq",
            "sample": sample_name,
            "cell_line": cell_line,
            "reference_source": dataset_name,
            "data_origin": "geo_control_counts",
        },
    )
    adata.var_names_make_unique()
    prefixed = [f"{sample_name}:{barcode}" for barcode in adata.obs_names.astype(str)]
    adata.obs_names = prefixed
    return adata


def _read_sparse_10x_sample(
    *,
    matrix_path: Path,
    barcodes_path: Path,
    features_path: Path,
    protospacer_path: Path | None,
    sample_name: str,
    cell_line: str,
    dataset_name: str,
) -> AnnData:
    with gzip.open(matrix_path, "rb") as handle:
        matrix = spio.mmread(handle).transpose().tocsr()
    features = pd.read_csv(features_path, sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="gzip")

    gene_names = features.iloc[:, 1].astype(str) if features.shape[1] > 1 else features.iloc[:, 0].astype(str)
    feature_ids = features.iloc[:, 0].astype(str)
    feature_types = features.iloc[:, 2].astype(str) if features.shape[1] > 2 else pd.Series(["Gene Expression"] * len(features))
    obs_index = _normalize_obs_index(barcodes.iloc[:, 0].astype(str))
    var = pd.DataFrame(
        {
            "feature_id": feature_ids.to_numpy(),
            "feature_type": feature_types.to_numpy(),
        },
        index=gene_names.astype(str),
    )
    var.index.name = None
    obs = pd.DataFrame(
        {
            "sample": sample_name,
            "cell_line": cell_line,
            "reference_source": dataset_name,
            "cell_type": "unknown",
            "roi": "global",
        },
        index=obs_index.astype(str),
    )

    if protospacer_path is not None and protospacer_path.exists():
        calls = pd.read_csv(protospacer_path, compression="gzip")
        if "cell_barcode" in calls.columns:
            calls["cell_barcode"] = _normalize_obs_index(calls["cell_barcode"].astype(str))
            calls = calls.set_index("cell_barcode")
        else:
            calls.index = _normalize_obs_index(calls.index.astype(str))
        rename_map = {}
        if "feature_call" in calls.columns:
            rename_map["feature_call"] = "guide_id"
        calls = calls.rename(columns=rename_map)
        obs = obs.join(calls, how="left")
    obs = _finalize_perturbation_obs(obs, guide_column="guide_id", feature_count_column="num_features")

    adata = from_tables(
        matrix,
        obs=obs,
        var=var,
        metadata={
            "dataset_name": dataset_name,
            "platform": "cropseq",
            "sample": sample_name,
            "cell_line": cell_line,
            "reference_source": dataset_name,
            "data_origin": "geo_sparse_10x",
        },
    )
    adata.var_names_make_unique()
    prefixed = [f"{sample_name}:{barcode}" for barcode in adata.obs_names.astype(str)]
    adata.obs_names = prefixed
    return adata


def _prepare_flat_geo_10x_dataset(card: DatasetCard, layout: dict[str, Path]) -> dict[str, Any]:
    raw_archive = layout["raw_dir"] / str(card.raw_archive_name)
    extracted_root = layout["raw_dir"] / "extracted"
    extracted_files = _extract_tar_once(raw_archive, extracted_root)
    file_lookup = {path.name: path for path in extracted_files}

    control_count_paths: list[tuple[str, str, Path]] = []
    sparse_samples: dict[str, dict[str, Path | str]] = {}
    for file_name, path in file_lookup.items():
        control_match = _CONTROL_COUNT_RE.match(file_name)
        if control_match:
            control_sample = str(control_match.group("sample"))
            control_count_paths.append(
                (
                    f"{control_match.group('gsm')}_{control_sample}_control",
                    control_sample.split("_", 1)[0],
                    path,
                )
            )
            continue

        sparse_match = _FLAT_10X_FILE_RE.match(file_name)
        if sparse_match:
            sample_name = f"{sparse_match.group('gsm')}_{sparse_match.group('sample')}"
            sparse_samples.setdefault(sample_name, {"sample_label": sparse_match.group("sample")})
            sparse_samples[sample_name][str(sparse_match.group("kind")).lower()] = path

    adatas: list[AnnData] = []
    for sample_name, sample_label, path in sorted(control_count_paths):
        adatas.append(_read_control_counts_table(path, sample_name=sample_name, cell_line=sample_label, dataset_name=card.name))

    for sample_name, payload in sorted(sparse_samples.items()):
        required = {"matrix", "barcodes", "features"}
        if missing := required.difference(payload):
            raise FileNotFoundError(f"Sample {sample_name} is missing required sparse inputs: {sorted(missing)}")
        sample_label = str(payload["sample_label"])
        adatas.append(
            _read_sparse_10x_sample(
                matrix_path=Path(payload["matrix"]),
                barcodes_path=Path(payload["barcodes"]),
                features_path=Path(payload["features"]),
                protospacer_path=Path(payload["protospacer_calls_per_cell"]) if "protospacer_calls_per_cell" in payload else None,
                sample_name=sample_name,
                cell_line=sample_label.split("_", 1)[0],
                dataset_name=card.name,
            )
        )

    if not adatas:
        raise FileNotFoundError(f"Could not find any parseable control counts or sparse 10x samples in {raw_archive}")

    adata = ad.concat(adatas, join="outer", merge="same")
    adata.var_names_make_unique()
    adata.uns["spatialperturb"] = {
        "dataset_name": card.name,
        "platform": card.platform,
        "source_accession": card.accession,
        "parser": card.parser,
        "sample_count": len(adatas),
        "guide_normalization": {
            "suffix_regex": _GUIDE_SUFFIX_RE.pattern,
            "control_regex": _CONTROL_GUIDE_RE.pattern,
        },
    }
    ensure_spatialperturb_schema(adata)
    validate_spatialperturb_schema(adata)
    prepared_path = layout["prepared_dir"] / str(card.prepared_filename)
    prepared_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(prepared_path)
    return {
        "status": "ready",
        "prepared_path": prepared_path,
        "prepared_kind": card.prepared_kind,
        "raw_dir": layout["raw_dir"],
        "sample_count": len(adatas),
    }


def _r_converter_script_path() -> Path:
    return Path(__file__).resolve().parent / "resources" / "convert_seurat_to_tables.R"


def _run_r_table_export(input_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    done_marker = output_dir / "export.done"
    if done_marker.exists():
        return

    script_path = _r_converter_script_path()
    command = ["Rscript", str(script_path), str(input_path), str(output_dir)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "R-backed Seurat conversion failed.\n"
            f"command={' '.join(command)}\n"
            f"stdout={result.stdout}\n"
            f"stderr={result.stderr}"
        )
    done_marker.write_text("ok\n", encoding="utf-8")


def _find_obs_column(obs: pd.DataFrame, candidates: tuple[str, ...], *, contains: tuple[str, ...] = ()) -> str | None:
    lower_to_original = {str(column).lower(): str(column) for column in obs.columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    for column in obs.columns.astype(str):
        lower = column.lower()
        if contains and any(token in lower for token in contains):
            return str(column)
    return None


def _infer_cell_line_series(obs: pd.DataFrame) -> pd.Series:
    explicit = _find_obs_column(
        obs,
        ("cell_line", "cellline", "line", "cell.line", "CellLine", "cellLine"),
        contains=("cell_line", "cellline"),
    )
    if explicit is not None:
        return obs[explicit].astype(str)

    object_columns = [column for column in obs.columns if pd.api.types.is_object_dtype(obs[column]) or pd.api.types.is_string_dtype(obs[column])]
    text_frame = obs.loc[:, object_columns].fillna("").astype(str)
    values = []
    for _, row in text_frame.iterrows():
        joined = " ".join(row.tolist())
        matched = next((cell_line for cell_line in _KNOWN_CELL_LINES if cell_line.lower() in joined.lower()), "unknown")
        values.append(matched)
    return pd.Series(values, index=obs.index, dtype="object")


def _infer_stimulus_from_name(file_name: str) -> str:
    for stimulus in ("IFNB", "IFNG", "INS", "TGFB", "TNFA"):
        if stimulus in file_name.upper():
            return stimulus
    return "unknown"


def _standardize_pathway_atlas_obs(
    obs: pd.DataFrame,
    *,
    file_name: str,
    dataset_name: str,
) -> pd.DataFrame:
    target = obs.copy()
    target.index = target.index.astype(str)

    sample_col = _find_obs_column(target, ("sample", "orig.ident", "orig_ident", "replicate", "batch"), contains=("sample", "orig", "replicate", "batch"))
    if sample_col is not None:
        target["sample"] = target[sample_col].astype(str)
    else:
        target["sample"] = file_name

    target["cell_line"] = _infer_cell_line_series(target).astype(str)
    target["stimulus"] = _infer_stimulus_from_name(file_name)

    guide_col = _find_obs_column(
        target,
        ("guide_id", "guide", "feature_call", "sgRNA", "sgrna", "gRNA", "grna", "target_sgRNA", "target_sgrna"),
        contains=("guide", "grna", "sgrna", "feature_call"),
    )
    if guide_col is None:
        raise ValueError(
            "Could not identify a guide assignment column in the exported Seurat metadata. "
            "Expected one of guide_id/guide/feature_call/sgRNA/gRNA-like columns."
        )
    target = target.rename(columns={guide_col: "guide_id"})

    num_features_col = _find_obs_column(target, ("num_features", "num_guides"), contains=("num_features", "num_guides"))
    target = _finalize_perturbation_obs(target, guide_column="guide_id", feature_count_column=num_features_col)
    target["reference_source"] = dataset_name
    target["cell_type"] = target.get("cell_type", "unknown").astype(str) if "cell_type" in target.columns else "unknown"
    target["roi"] = target.get("roi", "global").astype(str) if "roi" in target.columns else "global"
    return target


def _load_exported_seurat_tables(path: Path, *, dataset_name: str) -> AnnData:
    obs = pd.read_csv(path / "obs.csv", index_col=0)
    var = pd.read_csv(path / "var.csv", index_col=0)
    matrix = spio.mmread(path / "matrix.mtx").transpose().tocsr()
    obs = _standardize_pathway_atlas_obs(obs, file_name=path.name, dataset_name=dataset_name)
    var.index = var.index.astype(str)
    adata = from_tables(
        matrix,
        obs=obs,
        var=var,
        metadata={
            "dataset_name": dataset_name,
            "platform": "perturbseq_pathway_atlas",
            "reference_source": dataset_name,
            "data_origin": "seurat_rds_export",
            "stimulus": obs["stimulus"].iloc[0] if not obs.empty else "unknown",
        },
    )
    stimulus = obs["stimulus"].iloc[0] if not obs.empty else path.name
    prefixed = [f"{stimulus}:{barcode}" for barcode in adata.obs_names.astype(str)]
    adata.obs_names = prefixed
    return adata


def _prepare_rds_collection_dataset(card: DatasetCard, layout: dict[str, Path]) -> dict[str, Any]:
    raw_dir = layout["raw_dir"]
    extracted_root = raw_dir / "extracted"
    adatas: list[AnnData] = []
    exported_dirs: list[Path] = []
    for filename in card.supplementary_files:
        input_path = raw_dir / filename
        if not input_path.exists():
            raise FileNotFoundError(f"Missing expected RDS file for {card.name}: {input_path}")
        export_name = filename.removesuffix(".gz").removesuffix(".rds")
        export_dir = extracted_root / export_name
        _run_r_table_export(input_path, export_dir)
        exported_dirs.append(export_dir)
        adatas.append(_load_exported_seurat_tables(export_dir, dataset_name=card.name))

    if not adatas:
        raise FileNotFoundError(f"No Seurat exports were produced for {card.name}")

    adata = ad.concat(adatas, join="outer", merge="same")
    adata.var_names_make_unique()
    adata.uns["spatialperturb"] = {
        "dataset_name": card.name,
        "platform": card.platform,
        "source_accession": card.accession,
        "parser": card.parser,
        "sample_count": len(adatas),
        "supplementary_files": list(card.supplementary_files),
    }
    ensure_spatialperturb_schema(adata)
    validate_spatialperturb_schema(adata)
    prepared_path = layout["prepared_dir"] / str(card.prepared_filename)
    prepared_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(prepared_path)
    return {
        "status": "ready",
        "prepared_path": prepared_path,
        "prepared_kind": card.prepared_kind,
        "raw_dir": layout["raw_dir"],
        "exported_dirs": exported_dirs,
    }


def _prepare_nested_10x_dataset(card: DatasetCard, layout: dict[str, Path]) -> dict[str, Any]:
    raw_archive = layout["raw_dir"] / str(card.raw_archive_name)
    extracted_root = layout["raw_dir"] / "extracted"
    sample_dirs = _extract_nested_10x_archives(raw_archive, extracted_root)
    if not sample_dirs:
        raise FileNotFoundError(f"No sample archives were extracted from {raw_archive}")

    adatas: list[AnnData] = []
    barcode_columns: list[str] = []
    barcode_map: dict[str, str] = {}
    for sample_root in sample_dirs:
        matrix_dir = _sample_matrix_dir(sample_root)
        sample_name = sample_root.name
        sample_adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True, gex_only=False)
        sample_adata.obs_names = [f"{sample_name}:{obs_name}" for obs_name in sample_adata.obs_names]
        sample_adata.obs["sample"] = sample_name
        sample_adata.obs["cell_type"] = "unknown"
        sample_adata.obs["roi"] = "global"
        feature_names = list(map(str, sample_adata.var_names))
        sample_barcode_columns = [feature for feature in feature_names if feature.endswith("_gRNA")]
        barcode_columns.extend(sample_barcode_columns)
        barcode_map.update(_infer_gRNA_map(sample_barcode_columns))
        adatas.append(sample_adata)

    adata = ad.concat(adatas, join="outer", merge="same")
    adata.var_names_make_unique()
    adata.uns["spatialperturb"] = {
        "dataset_name": card.name,
        "platform": card.platform,
        "source_accession": card.accession,
        "sample_count": len(sample_dirs),
        "parser": card.parser,
    }
    ensure_spatialperturb_schema(adata)
    if barcode_map:
        assign_perturbations(
            adata,
            barcode_columns=sorted(set(barcode_columns)),
            barcode_to_perturbation=barcode_map,
            negative_label="unassigned",
            multiple_label="multiple",
        )
    adata.uns["spatialperturb"]["prepared_kind"] = card.prepared_kind
    prepared_path = layout["prepared_dir"] / str(card.prepared_filename)
    prepared_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(prepared_path)
    return {
        "status": "ready",
        "prepared_path": prepared_path,
        "prepared_kind": card.prepared_kind,
        "raw_dir": layout["raw_dir"],
    }


def _prepare_gef_passthrough_dataset(card: DatasetCard, layout: dict[str, Path]) -> dict[str, Any]:
    raw_archive = layout["raw_dir"] / str(card.raw_archive_name)
    extracted_root = layout["raw_dir"] / "extracted"
    _extract_tar_once(raw_archive, extracted_root)

    guessed = _guess_table_pair(layout["raw_dir"]) or _guess_table_pair(extracted_root)
    if guessed is not None:
        counts_path, metadata_path = guessed
        if counts_path.suffix == ".h5ad":
            adata = ad.read_h5ad(counts_path)
        else:
            sep = "\t" if ".tsv" in "".join(counts_path.suffixes) else ","
            expression = pd.read_csv(counts_path, sep=sep, index_col=0)
            metadata = pd.read_csv(metadata_path, sep=sep, index_col=0)
            adata = from_tables(expression, obs=metadata)
        ensure_spatialperturb_schema(
            adata,
            metadata={
                "dataset_name": card.name,
                "platform": card.platform,
                "source_accession": card.accession,
                "parser": "preconverted_tables",
            },
        )
        validate_spatialperturb_schema(adata)
        prepared_path = layout["prepared_dir"] / str(card.prepared_filename)
        prepared_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(prepared_path)
        return {
            "status": "ready",
            "prepared_path": prepared_path,
            "prepared_kind": card.prepared_kind,
            "raw_dir": layout["raw_dir"],
        }

    note = (
        "Raw GEF files were downloaded and extracted, but automatic conversion is not available in the core "
        "package. Place a preconverted .h5ad file or counts/cells tables in the dataset raw directory and rerun prepare_dataset."
    )
    placeholder = layout["prepared_dir"] / "PREPARE_NOTES.txt"
    placeholder.parent.mkdir(parents=True, exist_ok=True)
    placeholder.write_text(note, encoding="utf-8")
    return {
        "status": "requires_external_conversion",
        "prepared_path": None,
        "prepared_kind": card.prepared_kind,
        "raw_dir": layout["raw_dir"],
        "note": note,
    }


def available_datasets() -> pd.DataFrame:
    """Return the bundled dataset catalog."""
    rows = []
    for card in _DATASET_CATALOG:
        payload = asdict(card)
        if card.accession == "synthetic":
            payload["download_url"] = None
            payload["matrix_url"] = None
            payload["supplementary_urls"] = []
        else:
            payload["download_url"] = _geo_raw_url(card.accession) if card.raw_archive_name is not None else None
            payload["matrix_url"] = _geo_matrix_url(card.accession) if card.matrix_filename is not None else None
            payload["supplementary_urls"] = [_geo_supplementary_url(card.accession, name) for name in card.supplementary_files]
        rows.append(payload)
    return pd.DataFrame(rows)


def get_dataset_card(name: str) -> DatasetCard:
    """Fetch a dataset card by name."""
    for card in _DATASET_CATALOG:
        if card.name == name:
            return card
    raise KeyError(f"Dataset {name!r} is not registered.")


def fetch_dataset(name: str, *, cache_dir: str | Path, force: bool = False) -> dict[str, Any]:
    """Download and cache the raw public files for a registered dataset."""
    card = get_dataset_card(name)
    layout = _dataset_layout(name, cache_dir)
    layout["root"].mkdir(parents=True, exist_ok=True)
    layout["raw_dir"].mkdir(parents=True, exist_ok=True)

    if card.accession == "synthetic":
        payload = {
            "dataset": card.name,
            "status": "built_in",
            "raw_dir": layout["raw_dir"],
            "downloaded_files": [],
        }
        _write_json(layout["metadata_path"], {"card": asdict(card), "fetch": payload})
        return payload

    downloaded_files: list[Path] = []
    download_urls: list[str] = []
    if card.raw_archive_name is not None:
        raw_archive = layout["raw_dir"] / card.raw_archive_name
        downloaded_files.append(_download_file(_geo_raw_url(card.accession), raw_archive, force=force))
        download_urls.append(_geo_raw_url(card.accession))
    if card.matrix_filename is not None:
        matrix_file = layout["raw_dir"] / card.matrix_filename
        downloaded_files.append(_download_file(_geo_matrix_url(card.accession), matrix_file, force=force))
        download_urls.append(_geo_matrix_url(card.accession))
    for filename in card.supplementary_files:
        file_path = layout["raw_dir"] / filename
        file_url = _geo_supplementary_url(card.accession, filename)
        downloaded_files.append(_download_file(file_url, file_path, force=force))
        download_urls.append(file_url)

    payload = {
        "dataset": card.name,
        "status": "downloaded",
        "raw_dir": layout["raw_dir"],
        "downloaded_files": downloaded_files,
        "download_urls": download_urls,
    }
    _write_json(layout["metadata_path"], {"card": asdict(card), "fetch": payload})
    return payload


def prepare_dataset(
    name: str,
    *,
    cache_dir: str | Path,
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Standardize a public dataset into SpatialPerturb-compatible prepared outputs."""
    card = get_dataset_card(name)
    layout = _dataset_layout(name, cache_dir)
    layout["root"].mkdir(parents=True, exist_ok=True)
    layout["prepared_dir"].mkdir(parents=True, exist_ok=True)
    if output_dir is not None:
        layout["prepared_dir"] = Path(output_dir).expanduser().resolve()
        layout["prepared_dir"].mkdir(parents=True, exist_ok=True)

    if card.parser == "demo":
        adata = load_demo_dataset()
        prepared_path = layout["prepared_dir"] / str(card.prepared_filename)
        adata.write_h5ad(prepared_path)
        payload = {
            "dataset": card.name,
            "status": "ready",
            "prepared_path": prepared_path,
            "prepared_kind": card.prepared_kind,
        }
    else:
        expected_files = []
        if card.raw_archive_name is not None:
            expected_files.append(layout["raw_dir"] / card.raw_archive_name)
        expected_files.extend(layout["raw_dir"] / filename for filename in card.supplementary_files)
        if any(not path.exists() for path in expected_files):
            fetch_dataset(name, cache_dir=cache_dir, force=False)

        if card.parser == "geo_nested_10x_tar":
            payload = _prepare_nested_10x_dataset(card, layout)
        elif card.parser == "geo_gef_passthrough":
            payload = _prepare_gef_passthrough_dataset(card, layout)
        elif card.parser == "geo_flat_10x_tar":
            payload = _prepare_flat_geo_10x_dataset(card, layout)
        elif card.parser == "geo_rds_collection":
            payload = _prepare_rds_collection_dataset(card, layout)
        else:
            raise NotImplementedError(f"Unsupported dataset parser: {card.parser}")
        payload["dataset"] = card.name

    _write_json(layout["metadata_path"], {"card": asdict(card), "prepare": payload})
    return payload


def load_public_dataset(
    name: str,
    *,
    cache_dir: str | Path,
    prepared: bool = True,
) -> AnnData:
    """Load a prepared public dataset as an AnnData object."""
    card = get_dataset_card(name)
    if card.parser == "demo":
        return load_demo_dataset()

    if not prepared:
        raise ValueError("load_public_dataset requires prepared=True for public datasets.")

    layout = _dataset_layout(name, cache_dir)
    prepared_path = layout["prepared_dir"] / str(card.prepared_filename)
    if not prepared_path.exists():
        result = prepare_dataset(name, cache_dir=cache_dir)
        if result.get("prepared_path") is None:
            raise FileNotFoundError(
                f"Prepared dataset for {name!r} is not available yet. {result.get('note', '')}".strip()
            )
        prepared_path = Path(result["prepared_path"])

    adata = ad.read_h5ad(prepared_path)
    ensure_spatialperturb_schema(
        adata,
        metadata={
            "dataset_name": card.name,
            "platform": card.platform,
            "prepared_path": str(prepared_path),
        },
    )
    validate_spatialperturb_schema(adata)
    return adata


def _generate_demo_expression(platform: str, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed if platform == "spatial" else seed + 13)
    genes = [
        "STAT1",
        "IRF1",
        "CXCL10",
        "TOP2A",
        "MKI67",
        "LRRK2",
        "SRF",
        "LIG1",
        "REC1",
        "CTRL_BARCODE",
        "LRRK2_BARCODE",
        "SRF_BARCODE",
    ]
    cells = [f"{platform}_cell_{idx:02d}" for idx in range(20)]
    base = rng.poisson(lam=2.0, size=(len(cells), len(genes))).astype(float)
    expression = pd.DataFrame(base, index=cells, columns=genes)
    expression.loc[:, ["CTRL_BARCODE", "LRRK2_BARCODE", "SRF_BARCODE"]] = 0.0

    control_cells = cells[:4]
    lrrk2_cells = cells[4:8]
    srf_cells = cells[8:12]
    multiple_cells = cells[12:14]
    unassigned_cells = cells[14:]

    expression.loc[control_cells, "CTRL_BARCODE"] = 10
    expression.loc[lrrk2_cells, "LRRK2_BARCODE"] = 12
    expression.loc[srf_cells, "SRF_BARCODE"] = 12
    expression.loc[multiple_cells, ["LRRK2_BARCODE", "SRF_BARCODE"]] = 9
    expression.loc[unassigned_cells, ["CTRL_BARCODE", "LRRK2_BARCODE", "SRF_BARCODE"]] = 0

    expression.loc[lrrk2_cells, ["STAT1", "IRF1", "CXCL10", "LIG1"]] += 5
    expression.loc[lrrk2_cells, "LRRK2"] = np.clip(expression.loc[lrrk2_cells, "LRRK2"] - 1.5, 0, None)
    expression.loc[srf_cells, ["TOP2A", "MKI67"]] += 5
    expression.loc[srf_cells, "SRF"] = np.clip(expression.loc[srf_cells, "SRF"] - 1.5, 0, None)

    if platform == "spatial":
        expression.loc[cells[14:18], "REC1"] += 3
        expression.loc[cells[16:20], "CXCL10"] += 2
    else:
        expression.loc[lrrk2_cells, "REC1"] += 2
        expression.loc[srf_cells, "CXCL10"] += 2
    return expression


def _generate_demo_obs(platform: str) -> pd.DataFrame:
    cells = [f"{platform}_cell_{idx:02d}" for idx in range(20)]
    cell_type = ["neuron"] * 12 + ["astrocyte"] * 8
    roi = ["hippocampus"] * 18 + ["cortex"] * 2
    x_coords = [0, 0, 1, 1, 4, 4, 5, 5, 8, 8, 9, 9, 12, 12, 2, 2, 6, 6, 10, 10]
    y_coords = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    sample = ["S1", "S2", "S3", "S4"] * 5
    cell_line = ["MCF7"] * 20 if platform == "reference" else ["xenium_breast"] * 20
    stimulus = ["TNFA"] * 20 if platform == "reference" else ["baseline"] * 20
    return pd.DataFrame(
        {
            "cell_type": cell_type,
            "roi": roi,
            "sample": sample,
            "cell_line": cell_line,
            "stimulus": stimulus,
            "reference_source": f"demo_{platform}",
            "x": x_coords,
            "y": y_coords,
        },
        index=cells,
    )


def load_demo_dataset(
    *,
    platform: str = "spatial",
    paired: bool = False,
    annotate: bool = True,
    seed: int = 0,
) -> AnnData | tuple[AnnData, AnnData]:
    """Return a deterministic demo dataset for tests, docs, and end-to-end examples."""
    if paired:
        return (
            load_demo_dataset(platform="spatial", annotate=annotate, seed=seed),
            load_demo_dataset(platform="reference", annotate=annotate, seed=seed),
        )

    if platform not in {"spatial", "reference"}:
        raise ValueError("platform must be 'spatial', 'reference', or use paired=True.")

    expression = _generate_demo_expression(platform=platform, seed=seed)
    obs = _generate_demo_obs(platform)
    adata = from_tables(
        expression,
        obs=obs,
        spatial=obs[["x", "y"]],
        metadata={"platform": platform, "dataset_name": "demo_spatialperturb"},
    )
    if annotate:
        assign_perturbations(
            adata,
            barcode_columns=["CTRL_BARCODE", "LRRK2_BARCODE", "SRF_BARCODE"],
            barcode_to_perturbation={
                "CTRL_BARCODE": "control",
                "LRRK2_BARCODE": "Lrrk2",
                "SRF_BARCODE": "Srf",
            },
        )
        adata.obs["guide_id"] = adata.obs["perturbation"].astype(str)
        adata.obs["reference_source"] = f"demo_{platform}"
    build_spatial_graph(adata, mode="knn", k=5)
    return adata
