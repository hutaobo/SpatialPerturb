"""Run the GSE274058 reference-side analysis package."""

from __future__ import annotations

from datetime import datetime, timezone
import importlib.metadata
import json
import platform
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import spatialperturb as sp


CACHE_DIR = REPO_ROOT / ".spatialperturb-cache"
OUTPUT_DIR = REPO_ROOT / "reports" / "gse274058_reference_run"
CONTROL = "control"
MIN_CELLS = 10
MIN_SAMPLES = 2
METHOD = "pseudobulk"
SAMPLE_COL = "sample"
PROGRAM_TOP_N = 50
TOP_HITS_PER_PERTURBATION = 10

TARGET_MAP = {
    "Trem2": "Trem2",
    "Rraga": "Rraga",
    "Myrf": "Myrf",
    "Fasn": "Fasn",
    "Clu": "Clu",
    "Dpp6": "Dpp6",
    "Tbk1": "Tbk1",
    "Flcn": "Flcn",
    "Gfap": "Gfap",
    "C9orf72": "C9orf72",
    "Cfap410": "Cfap410",
    "Stk39": "Stk39",
    "Lrrk2": "Lrrk2",
    "Ndufaf2": "Ndufaf2",
    "Sh3gl2": "Sh3gl2",
    "Srf": "Srf",
    "Rbfox3": "Rbfox3",
    "Olig2": "Olig2",
}


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=_json_default), encoding="utf-8")


def _write_table(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def _status_counts(adata) -> dict[str, int]:
    return adata.obs["perturbation_status"].astype(str).value_counts().sort_index().astype(int).to_dict()


def _perturbation_counts(adata) -> dict[str, int]:
    single = adata.obs.loc[adata.obs["perturbation_status"].astype(str) == "single", "perturbation"].astype(str)
    return single.value_counts().sort_index().astype(int).to_dict()


def _focus_perturbations(de_results: pd.DataFrame) -> list[str]:
    preferred = [name for name in ["Lrrk2", "Srf"] if name in set(de_results["perturbation"].astype(str))]
    if preferred:
        return preferred
    ranking = (
        de_results.assign(abs_log2fc=lambda df: df["log2fc"].abs())
        .groupby("perturbation")
        .agg(best_fdr=("fdr", "min"), best_abs_log2fc=("abs_log2fc", "max"))
        .sort_values(["best_fdr", "best_abs_log2fc"], ascending=[True, False])
    )
    return ranking.head(2).index.astype(str).tolist()


def _render_barcode_plot(adata, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    sp.pl.barcode_spread(adata, ax=ax)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _cell_and_sample_counts(adata, perturbation: str) -> tuple[int, int]:
    mask = (
        (adata.obs["perturbation"].astype(str) == str(perturbation))
        & (adata.obs["perturbation_status"].astype(str) == "single")
    )
    cells = int(mask.sum())
    samples = int(adata.obs.loc[mask, SAMPLE_COL].astype(str).nunique()) if SAMPLE_COL in adata.obs.columns else 0
    return cells, samples


def _dependency_versions() -> dict[str, str]:
    packages = [
        "anndata",
        "matplotlib",
        "numpy",
        "pandas",
        "scanpy",
        "scipy",
        "seaborn",
        "sklearn",
        "spatialperturb",
        "statsmodels",
        "typer",
    ]
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            if package == "spatialperturb":
                versions[package] = getattr(sp, "__version__", "unknown")
            else:
                versions[package] = "not-installed"
    return versions


def _summary_markdown(
    *,
    adata,
    qc: pd.DataFrame,
    de_results: pd.DataFrame,
    top_hits: pd.DataFrame,
    skipped: pd.DataFrame,
    focus: list[str],
) -> str:
    lines: list[str] = []
    lines.append("# GSE274058 reference run")
    lines.append("")
    lines.append("## Dataset summary")
    lines.append(f"- Cells: {adata.n_obs}")
    lines.append(f"- Genes: {adata.n_vars}")
    lines.append(f"- Samples: {adata.obs[SAMPLE_COL].astype(str).nunique() if SAMPLE_COL in adata.obs.columns else 0}")
    lines.append(f"- Barcode status: {_status_counts(adata)}")
    valid = qc.loc[
        (qc["valid_for_inference"]) & (qc["perturbation"].astype(str) != CONTROL),
        ["perturbation", "n_cells"],
    ].copy()
    lines.append("")
    lines.append("## Valid perturbations")
    if valid.empty:
        lines.append("No perturbations passed inference QC.")
    else:
        for row in valid.sort_values("n_cells", ascending=False).itertuples(index=False):
            lines.append(f"- {row.perturbation}: {row.n_cells} cells")

    lines.append("")
    lines.append("## Focus perturbations")
    if top_hits.empty:
        lines.append("No successful intrinsic DE results were produced.")
    else:
        for perturbation in focus:
            subset = top_hits[top_hits["perturbation"].astype(str) == perturbation].head(5)
            if subset.empty:
                continue
            lines.append(f"### {perturbation}")
            for row in subset.itertuples(index=False):
                lines.append(f"- {row.gene}: log2fc={row.log2fc:.3f}, fdr={row.fdr:.3g}")
            lines.append("")

    if not skipped.empty:
        lines.append("## Skipped perturbations")
        for row in skipped.itertuples(index=False):
            lines.append(
                f"- {row.perturbation}: {row.reason} "
                f"(case_n={row.case_n}, control_n={row.control_n}, case_sample_n={row.case_sample_n}, control_sample_n={row.control_sample_n})"
            )
        lines.append("")
    return "\n".join(lines).strip() + "\n"


def main() -> None:
    started_at = datetime.now(timezone.utc)

    print("[1/9] Fetching GSE274058 raw files...", flush=True)
    fetch_result = sp.fetch_dataset("shen_2026_scrnaseq", cache_dir=CACHE_DIR)

    print("[2/9] Preparing GSE274058 into AnnData...", flush=True)
    prepare_result = sp.prepare_dataset("shen_2026_scrnaseq", cache_dir=CACHE_DIR)
    prepared_path = prepare_result.get("prepared_path")
    if prepared_path is None:
        raise RuntimeError(f"prepare_dataset did not produce a prepared path: {prepare_result}")

    print("[3/9] Loading prepared dataset...", flush=True)
    adata = sp.load_public_dataset("shen_2026_scrnaseq", cache_dir=CACHE_DIR)

    print("[4/9] Validating schema and writing QC outputs...", flush=True)
    sp.validate_spatialperturb_schema(adata)
    qc_dir = OUTPUT_DIR / "qc"
    de_dir = OUTPUT_DIR / "de"
    programs_dir = OUTPUT_DIR / "programs"
    for path in [qc_dir, de_dir, programs_dir]:
        path.mkdir(parents=True, exist_ok=True)

    qc = sp.qc_perturbations(adata, control=CONTROL, target_map=TARGET_MAP, min_cells=MIN_CELLS)
    _write_table(qc, qc_dir / "perturbation_qc.tsv")
    _render_barcode_plot(adata, qc_dir / "barcode_spread.png")

    print("[5/9] Selecting valid perturbations...", flush=True)
    valid_perturbations = (
        qc.loc[(qc["valid_for_inference"]) & (qc["perturbation"].astype(str) != CONTROL), "perturbation"]
        .astype(str)
        .tolist()
    )

    print(f"[6/9] Running intrinsic DE for {len(valid_perturbations)} perturbations...", flush=True)
    de_tables: list[pd.DataFrame] = []
    skipped_records: list[dict[str, Any]] = []
    for perturbation in valid_perturbations:
        case_n, case_sample_n = _cell_and_sample_counts(adata, perturbation)
        control_n, control_sample_n = _cell_and_sample_counts(adata, CONTROL)
        try:
            table = sp.intrinsic_de(
                adata,
                perturbation=perturbation,
                control=CONTROL,
                method=METHOD,
                sample_col=SAMPLE_COL,
                min_cells_per_group=MIN_CELLS,
                min_samples_per_group=MIN_SAMPLES,
            )
            de_tables.append(table)
            print(f"  - {perturbation}: ok ({len(table)} genes)", flush=True)
        except Exception as exc:
            skipped_records.append(
                {
                    "perturbation": perturbation,
                    "reason": str(exc),
                    "case_n": case_n,
                    "control_n": control_n,
                    "case_sample_n": case_sample_n,
                    "control_sample_n": control_sample_n,
                }
            )
            print(f"  - {perturbation}: skipped ({exc})", flush=True)

    de_results = pd.concat(de_tables, ignore_index=True) if de_tables else pd.DataFrame()
    skipped = pd.DataFrame.from_records(skipped_records)
    _write_table(de_results, de_dir / "intrinsic_de.tsv")
    _write_table(skipped, de_dir / "skipped_perturbations.tsv")

    print("[7/9] Deriving perturbation programs and top hits...", flush=True)
    if not de_results.empty:
        successful_perturbations = sorted(de_results["perturbation"].astype(str).unique().tolist())
        derived_programs = sp.derive_perturbation_programs(de_results, top_n=PROGRAM_TOP_N, direction="both")
        programs = {perturbation: list(derived_programs.get(perturbation, [])) for perturbation in successful_perturbations}
        program_matrix = sp.build_signature_matrix(programs)
        focus = _focus_perturbations(de_results)
        top_hits = (
            de_results[de_results["perturbation"].astype(str).isin(focus)]
            .assign(abs_log2fc=lambda df: df["log2fc"].abs())
            .sort_values(["perturbation", "fdr", "abs_log2fc"], ascending=[True, True, False])
            .groupby("perturbation", as_index=False, group_keys=False)
            .head(TOP_HITS_PER_PERTURBATION)
            .drop(columns=["abs_log2fc"])
            .reset_index(drop=True)
        )
    else:
        programs = {}
        program_matrix = pd.DataFrame()
        focus = []
        top_hits = pd.DataFrame()

    (programs_dir / "programs.json").write_text(json.dumps(programs, indent=2), encoding="utf-8")
    _write_table(program_matrix.reset_index(names="perturbation"), programs_dir / "program_matrix.tsv")
    _write_table(top_hits, de_dir / "top_hits.tsv")

    print("[8/9] Writing dataset summary and markdown report...", flush=True)
    finished_at = datetime.now(timezone.utc)
    dataset_summary = {
        "dataset": "shen_2026_scrnaseq",
        "accession": "GSE274058",
        "method": METHOD,
        "sample_col": SAMPLE_COL,
        "control": CONTROL,
        "min_cells_per_group": MIN_CELLS,
        "min_samples_per_group": MIN_SAMPLES,
        "fetch_result": fetch_result,
        "prepare_result": prepare_result,
        "prepared_path": prepared_path,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_samples": int(adata.obs[SAMPLE_COL].astype(str).nunique()) if SAMPLE_COL in adata.obs.columns else 0,
        "barcode_status_counts": _status_counts(adata),
        "single_cell_perturbation_counts": _perturbation_counts(adata),
        "valid_perturbations": valid_perturbations,
        "successful_perturbations": sorted(de_results["perturbation"].astype(str).unique().tolist()) if not de_results.empty else [],
        "skipped_perturbations": skipped_records,
        "runtime": {
            "generated_at_utc": finished_at.isoformat(),
            "started_at_utc": started_at.isoformat(),
            "duration_seconds": round((finished_at - started_at).total_seconds(), 3),
            "command": " ".join(sys.argv),
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "spatialperturb_version": getattr(sp, "__version__", "unknown"),
            "dependency_versions": _dependency_versions(),
        },
    }
    _write_json(OUTPUT_DIR / "dataset_summary.json", dataset_summary)
    summary = _summary_markdown(
        adata=adata,
        qc=qc,
        de_results=de_results,
        top_hits=top_hits,
        skipped=skipped,
        focus=focus,
    )
    (OUTPUT_DIR / "summary.md").write_text(summary, encoding="utf-8")

    print("[9/9] Done.", flush=True)
    print(f"Output written to: {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
