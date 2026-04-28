"""Compare two GSE274058 reference run directories and optionally replace the baseline."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, default=None)
    parser.add_argument("--replace-baseline", action="store_true")
    return parser.parse_args()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _choose_focus(top_hits: pd.DataFrame) -> list[str]:
    preferred = [name for name in ["Lrrk2", "Srf"] if name in set(top_hits["perturbation"].astype(str))]
    if preferred:
        return preferred
    ranking = (
        top_hits.groupby("perturbation")
        .agg(best_fdr=("fdr", "min"), best_abs_log2fc=("log2fc", lambda values: values.abs().max()))
        .sort_values(["best_fdr", "best_abs_log2fc"], ascending=[True, False])
    )
    return ranking.head(2).index.astype(str).tolist()


def _load_program_coverage(path: Path) -> dict[str, int]:
    table = pd.read_csv(path, sep="\t")
    index_column = "perturbation" if "perturbation" in table.columns else "program"
    coverage = (
        table.set_index(index_column)
        .fillna(0)
        .astype(int)
        .sum(axis=1)
        .astype(int)
        .to_dict()
    )
    return {str(key): int(value) for key, value in coverage.items()}


def _top_hits_signature(path: Path, focus: list[str]) -> list[dict[str, Any]]:
    table = pd.read_csv(path, sep="\t")
    subset = table[table["perturbation"].astype(str).isin(focus)].copy()
    subset = subset.sort_values(["perturbation", "fdr", "log2fc", "gene"], ascending=[True, True, False, True])
    signature: list[dict[str, Any]] = []
    for row in subset.itertuples(index=False):
        signature.append(
            {
                "perturbation": str(row.perturbation),
                "gene": str(row.gene),
                "log2fc": round(float(row.log2fc), 6),
                "fdr": round(float(row.fdr), 12),
            }
        )
    return signature


def _safe_replace_tree(src: Path, dst: Path) -> None:
    dst = dst.resolve()
    src = src.resolve()
    if not src.exists():
        raise FileNotFoundError(f"Candidate directory does not exist: {src}")
    if dst == dst.anchor:
        raise ValueError(f"Refusing to replace filesystem root: {dst}")
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def main() -> None:
    args = _parse_args()
    baseline_dir = args.baseline_dir.resolve()
    candidate_dir = args.candidate_dir.resolve()

    baseline_summary = _load_json(baseline_dir / "dataset_summary.json")
    candidate_summary = _load_json(candidate_dir / "dataset_summary.json")
    baseline_top_hits = pd.read_csv(baseline_dir / "de" / "top_hits.tsv", sep="\t")
    candidate_top_hits = pd.read_csv(candidate_dir / "de" / "top_hits.tsv", sep="\t")
    focus = sorted(set(_choose_focus(baseline_top_hits)) | set(_choose_focus(candidate_top_hits)))

    baseline_core = {
        "n_obs": baseline_summary.get("n_obs"),
        "n_vars": baseline_summary.get("n_vars"),
        "n_samples": baseline_summary.get("n_samples"),
        "barcode_status_counts": baseline_summary.get("barcode_status_counts"),
        "valid_perturbations": baseline_summary.get("valid_perturbations"),
    }
    candidate_core = {
        "n_obs": candidate_summary.get("n_obs"),
        "n_vars": candidate_summary.get("n_vars"),
        "n_samples": candidate_summary.get("n_samples"),
        "barcode_status_counts": candidate_summary.get("barcode_status_counts"),
        "valid_perturbations": candidate_summary.get("valid_perturbations"),
    }

    baseline_programs = _load_program_coverage(baseline_dir / "programs" / "program_matrix.tsv")
    candidate_programs = _load_program_coverage(candidate_dir / "programs" / "program_matrix.tsv")
    baseline_signature = _top_hits_signature(baseline_dir / "de" / "top_hits.tsv", focus)
    candidate_signature = _top_hits_signature(candidate_dir / "de" / "top_hits.tsv", focus)

    match = (
        baseline_core == candidate_core
        and baseline_programs == candidate_programs
        and baseline_signature == candidate_signature
    )

    status = "match" if match else "failed"
    if not match and args.replace_baseline:
        _safe_replace_tree(candidate_dir, baseline_dir)
        status = "replaced"

    payload = {
        "status": status,
        "compared_at_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_report_dir": str(baseline_dir),
        "candidate_report_dir": str(candidate_dir),
        "focus_perturbations": focus,
        "baseline_core": baseline_core,
        "candidate_core": candidate_core,
        "baseline_program_coverage": baseline_programs,
        "candidate_program_coverage": candidate_programs,
        "baseline_top_hits_signature": baseline_signature,
        "candidate_top_hits_signature": candidate_signature,
    }

    output = json.dumps(payload, indent=2)
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(output, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
