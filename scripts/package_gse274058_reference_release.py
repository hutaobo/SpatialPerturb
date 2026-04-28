"""Package GSE274058 reference outputs into release assets and docs-friendly summaries."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import shutil
import tarfile
from pathlib import Path
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT_DIR = REPO_ROOT / "reports" / "gse274058_reference_run"
DEFAULT_ARTIFACTS_DIR = REPO_ROOT / "artifacts" / "gse274058_reference_release"
DEFAULT_DOCS_DATA_DIR = REPO_ROOT / "docs" / "results" / "gse274058_reference"
DEFAULT_DOCS_STATIC_DIR = REPO_ROOT / "docs" / "_static" / "results" / "gse274058_reference"
REPO_RELEASE_PREFIX = "https://github.com/hutaobo/SpatialPerturb/releases/latest/download"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--artifacts-dir", type=Path, default=DEFAULT_ARTIFACTS_DIR)
    parser.add_argument("--docs-data-dir", type=Path, default=DEFAULT_DOCS_DATA_DIR)
    parser.add_argument("--docs-static-dir", type=Path, default=DEFAULT_DOCS_STATIC_DIR)
    parser.add_argument("--comparison-json", type=Path, default=None)
    parser.add_argument(
        "--authoritative-source",
        default=None,
        help="Human-readable source path for the authoritative run, e.g. an A100 report directory.",
    )
    return parser.parse_args()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _write_text(path: Path, content: str) -> None:
    _ensure_dir(path.parent)
    path.write_text(content, encoding="utf-8")


def _write_json(path: Path, payload: Any) -> None:
    _ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_table(path: Path, table: pd.DataFrame) -> None:
    _ensure_dir(path.parent)
    table.to_csv(path, sep="\t", index=False)


def _copy_file(src: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    shutil.copy2(src, dst)


def _gzip_copy(src: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    with src.open("rb") as source, gzip.open(dst, "wb", compresslevel=9) as target:
        shutil.copyfileobj(source, target)


def _tar_directory(src_dir: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    with tarfile.open(dst, "w:gz") as archive:
        archive.add(src_dir, arcname=src_dir.name)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _markdown_table(table: pd.DataFrame) -> str:
    if table.empty:
        return "_No rows available._\n"
    columns = [str(column) for column in table.columns]
    values = [[str(value) for value in row] for row in table.itertuples(index=False, name=None)]
    widths = [len(column) for column in columns]
    for row in values:
        widths = [max(width, len(cell)) for width, cell in zip(widths, row)]

    def format_row(row: list[str]) -> str:
        return "| " + " | ".join(cell.ljust(width) for cell, width in zip(row, widths)) + " |"

    header = format_row(columns)
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [format_row(row) for row in values]
    return "\n".join([header, separator, *body]) + "\n"


def _focus_targets(top_hits: pd.DataFrame) -> list[str]:
    preferred = [name for name in ["Lrrk2", "Srf"] if name in set(top_hits["perturbation"].astype(str))]
    if preferred:
        return preferred
    ranking = (
        top_hits.groupby("perturbation")
        .agg(best_fdr=("fdr", "min"), best_abs_log2fc=("log2fc", lambda values: values.abs().max()))
        .sort_values(["best_fdr", "best_abs_log2fc"], ascending=[True, False])
    )
    return ranking.head(2).index.astype(str).tolist()


def _format_float(value: Any, digits: int = 3) -> str:
    if pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}"


def _format_sci(value: Any) -> str:
    if pd.isna(value):
        return ""
    numeric = float(value)
    if numeric == 0.0:
        return "0"
    if abs(numeric) < 0.001 or abs(numeric) >= 1000:
        return f"{numeric:.3e}"
    return f"{numeric:.3g}"


def _load_comparison_status(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.exists():
        return None
    return _load_json(path)


def main() -> None:
    args = _parse_args()
    report_dir = args.report_dir.resolve()
    artifacts_dir = args.artifacts_dir.resolve()
    docs_data_dir = args.docs_data_dir.resolve()
    docs_static_dir = args.docs_static_dir.resolve()

    dataset_summary = _load_json(report_dir / "dataset_summary.json")
    qc = pd.read_csv(report_dir / "qc" / "perturbation_qc.tsv", sep="\t")
    top_hits = pd.read_csv(report_dir / "de" / "top_hits.tsv", sep="\t")
    program_matrix = pd.read_csv(report_dir / "programs" / "program_matrix.tsv", sep="\t")
    intrinsic_de = pd.read_csv(report_dir / "de" / "intrinsic_de.tsv", sep="\t")

    valid = (
        qc.loc[(qc["valid_for_inference"]) & (qc["perturbation"].astype(str) != "control")]
        .copy()
        .assign(knockdown_adequate=lambda df: df["target_log2fc"].fillna(0.0) < 0.0)
        .sort_values(["n_cells", "perturbation"], ascending=[False, True])
        .reset_index(drop=True)
    )
    valid["target_log2fc"] = valid["target_log2fc"].map(lambda value: _format_float(value, digits=3))
    valid["fraction_cells"] = valid["fraction_cells"].map(lambda value: _format_float(value, digits=3))
    valid["knockdown_adequate"] = valid["knockdown_adequate"].map(lambda value: "yes" if bool(value) else "no")
    valid_table = valid[["perturbation", "n_cells", "fraction_cells", "target_gene", "target_log2fc", "knockdown_adequate"]]

    focus_targets = _focus_targets(top_hits)
    focus_table = (
        top_hits.loc[top_hits["perturbation"].astype(str).isin(focus_targets), ["perturbation", "gene", "log2fc", "fdr", "case_n", "control_n"]]
        .copy()
    )
    focus_table["log2fc"] = focus_table["log2fc"].map(_format_float)
    focus_table["fdr"] = focus_table["fdr"].map(_format_sci)

    target_rows = intrinsic_de.loc[
        intrinsic_de.apply(lambda row: str(row["perturbation"]) == str(row["gene"]), axis=1),
        ["perturbation", "gene", "log2fc", "fdr", "mean_case", "mean_control"],
    ].copy()
    target_rows = target_rows[target_rows["perturbation"].astype(str).isin(["Lrrk2", "Srf"])]
    target_rows["log2fc"] = target_rows["log2fc"].map(_format_float)
    target_rows["fdr"] = target_rows["fdr"].map(_format_sci)
    target_rows["mean_case"] = target_rows["mean_case"].map(_format_float)
    target_rows["mean_control"] = target_rows["mean_control"].map(_format_float)

    index_column = "perturbation" if "perturbation" in program_matrix.columns else "program"
    program_gene_counts = (
        program_matrix.set_index(index_column)
        .fillna(0)
        .astype(int)
        .sum(axis=1)
        .rename("program_gene_count")
        .reset_index()
        .rename(columns={index_column: "perturbation"})
        .sort_values(["program_gene_count", "perturbation"], ascending=[False, True])
        .reset_index(drop=True)
    )

    comparison_status = _load_comparison_status(args.comparison_json)
    a100_status_lines = [
        "- Status: pending authoritative A100 rerun.",
        "- Expected rerun path: `/data/taobo.hu/SpatialPerturb`.",
        "- Current public summary reflects the local draft run committed with this branch.",
    ]
    if comparison_status is not None:
        a100_status_lines = [
            f"- Status: {comparison_status.get('status', 'unknown')}.",
            f"- Baseline report: `{comparison_status.get('baseline_report_dir', 'unknown')}`.",
            f"- Candidate report: `{comparison_status.get('candidate_report_dir', 'unknown')}`.",
            f"- Compared at (UTC): `{comparison_status.get('compared_at_utc', 'unknown')}`.",
        ]
        if args.authoritative_source:
            a100_status_lines.append(f"- Authoritative source: `{args.authoritative_source}`.")
        if comparison_status.get("status") == "replaced":
            a100_status_lines.append("- Outcome: A100 rerun replaced the local draft package for release.")

    runtime = dataset_summary.get("runtime", {})
    source_report_dir = args.authoritative_source or str(report_dir)
    overview_lines = [
        f"- Dataset card: `{dataset_summary['dataset']}` (`{dataset_summary['accession']}`).",
        "- Result type: dissociated reference-side pseudobulk intrinsic DE package.",
        f"- Generated at (UTC): `{runtime.get('generated_at_utc', 'unknown')}`.",
        f"- Duration: `{runtime.get('duration_seconds', 'unknown')}` seconds.",
        f"- Command: `{runtime.get('command', 'unknown')}`.",
        f"- Python: `{runtime.get('python_version', 'unknown')}` on `{runtime.get('platform', 'unknown')}`.",
        f"- SpatialPerturb version: `{runtime.get('spatialperturb_version', 'unknown')}`.",
        f"- Source report directory: `{source_report_dir}`.",
    ]
    qc_summary_lines = [
        f"- Cells: `{dataset_summary['n_obs']}`",
        f"- Genes: `{dataset_summary['n_vars']}`",
        f"- Samples: `{dataset_summary['n_samples']}`",
        f"- Barcode status counts: `{dataset_summary['barcode_status_counts']}`",
        f"- Single-cell control count: `{dataset_summary['single_cell_perturbation_counts'].get('control', 0)}`",
        f"- Valid perturbations: `{len(dataset_summary['valid_perturbations'])}`",
        f"- Successful perturbations: `{len(dataset_summary['successful_perturbations'])}`",
    ]

    knockdown_pass_n = int((valid["knockdown_adequate"] == "yes").sum())
    improvement_lines = [
        (
            f"- Target knockdown quality is unstable: only `{knockdown_pass_n}` of "
            f"`{len(valid_table)}` valid perturbations show `target_log2fc < 0`, and neither `Lrrk2` nor `Srf` "
            "shows a convincing target-gene decrease in this run. Add a `knockdown_adequate` QC flag and separate "
            "weak perturbations from main claims."
        ),
        (
            "- Sample and cell balance is tight: `control` has only `24` single cells across `3` samples, and several "
            "perturbations sit between `14` and `23` cells. Add per-sample minimum cell thresholds and stricter "
            "sample-level exclusion rules before treating reference signatures as stable."
        ),
        (
            "- The raw DE table is large for routine browsing (`intrinsic_de.tsv` is about `52.7 MB`). Keep publishing "
            "compressed `.tsv.gz` assets by default, and consider adding `.parquet` export for downstream reuse."
        ),
        (
            "- Documentation is still dual-stack today: RTD uses Sphinx while local site configuration still carries "
            "MkDocs nav. This publish keeps Sphinx for stability, but the long-term maintenance path should converge "
            "to one docs stack."
        ),
    ]

    _write_text(docs_data_dir / "overview.md", "\n".join(overview_lines) + "\n")
    _write_text(docs_data_dir / "qc_summary.md", "\n".join(qc_summary_lines) + "\n")
    _write_text(docs_data_dir / "a100_status.md", "\n".join(a100_status_lines) + "\n")
    _write_text(docs_data_dir / "valid_perturbations.md", _markdown_table(valid_table))
    _write_text(docs_data_dir / "top_hits.md", _markdown_table(focus_table))
    _write_text(docs_data_dir / "target_gene_sanity.md", _markdown_table(target_rows))
    _write_text(docs_data_dir / "program_summary.md", _markdown_table(program_gene_counts))
    _write_text(docs_data_dir / "improvement.md", "\n".join(improvement_lines) + "\n")

    _write_json(docs_data_dir / "dataset_summary.json", dataset_summary)
    _write_table(docs_data_dir / "valid_perturbations.tsv", valid_table)
    _write_table(docs_data_dir / "top_hits.tsv", focus_table)
    _write_table(docs_data_dir / "program_summary.tsv", program_gene_counts)
    _copy_file(report_dir / "qc" / "barcode_spread.png", docs_static_dir / "barcode_spread.png")

    intrinsic_gz = artifacts_dir / "gse274058_reference_run_intrinsic_de.tsv.gz"
    bundle_path = artifacts_dir / "gse274058_reference_run_bundle.tar.gz"
    sha_path = artifacts_dir / "SHA256SUMS.txt"
    manifest_path = artifacts_dir / "release_manifest.json"
    _gzip_copy(report_dir / "de" / "intrinsic_de.tsv", intrinsic_gz)
    _tar_directory(report_dir, bundle_path)

    sha_lines = []
    for path in [bundle_path, intrinsic_gz]:
        sha_lines.append(f"{_sha256(path)}  {path.name}")
    _write_text(sha_path, "\n".join(sha_lines) + "\n")

    release_manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "report_dir": str(report_dir),
        "artifacts_dir": str(artifacts_dir),
        "docs_data_dir": str(docs_data_dir),
        "docs_static_dir": str(docs_static_dir),
        "release_assets": {
            "bundle": {
                "path": str(bundle_path),
                "url": f"{REPO_RELEASE_PREFIX}/{bundle_path.name}",
                "sha256": _sha256(bundle_path),
            },
            "intrinsic_de_gz": {
                "path": str(intrinsic_gz),
                "url": f"{REPO_RELEASE_PREFIX}/{intrinsic_gz.name}",
                "sha256": _sha256(intrinsic_gz),
            },
            "sha256sums": {
                "path": str(sha_path),
                "url": f"{REPO_RELEASE_PREFIX}/{sha_path.name}",
                "sha256": _sha256(sha_path),
            },
        },
        "focus_perturbations": focus_targets,
        "knockdown_adequate_count": knockdown_pass_n,
    }
    _write_json(manifest_path, release_manifest)

    print(f"Docs summaries written to: {docs_data_dir}")
    print(f"Docs static assets written to: {docs_static_dir}")
    print(f"Release artifacts written to: {artifacts_dir}")


if __name__ == "__main__":
    main()
