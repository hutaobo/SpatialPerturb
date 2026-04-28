"""Summarize biological meaning of breast reference projection outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_csv(path, sep="\t", compression=compression)


def _program_theme(program: str, genes: list[str] | None = None) -> str:
    text = " ".join([program, *(genes or [])]).upper()
    if any(token in text for token in ("IFN", "STAT1", "IRF", "CXCL", "TNFA", "TNF")):
        return "immune/cytokine signaling, often compatible with interferon or inflammatory tumor microenvironment states"
    if any(token in text for token in ("TGFB", "TGF", "COL", "VIM", "FN1", "MMP")):
        return "TGF-beta, extracellular matrix, stromal remodeling, or EMT-like biology"
    if any(token in text for token in ("ERBB", "EGFR", "IGF1R", "PIK3", "AKT", "MAPK")):
        return "growth-factor and oncogenic signaling programs"
    if any(token in text for token in ("B2M", "HLA", "NLRC5", "TAP1", "TAP2")):
        return "antigen presentation and tumor-immune visibility"
    if any(token in text for token in ("ARID", "CREBBP", "EP300", "SMAR", "KMT", "HDAC")):
        return "chromatin regulation and transcriptional state control"
    if any(token in text for token in ("MSH", "MLH", "BRCA", "RAD", "ATM", "ATR", "CHEK")):
        return "DNA repair, genome stability, or replication stress"
    if any(token in text for token in ("MKI67", "TOP2A", "CDK", "E2F", "AURK")):
        return "cell cycle and proliferative activity"
    return "reference perturbation-associated transcriptional state"


def _top_reference_genes(reference_de: pd.DataFrame) -> dict[str, list[str]]:
    if reference_de.empty or "program" not in reference_de.columns or "gene" not in reference_de.columns:
        return {}
    score_col = "log2fc" if "log2fc" in reference_de.columns else None
    genes: dict[str, list[str]] = {}
    for program, frame in reference_de.groupby("program"):
        subset = frame.copy()
        if score_col is not None:
            subset = subset.assign(_abs=subset[score_col].abs()).sort_values("_abs", ascending=False)
        genes[str(program)] = subset["gene"].astype(str).drop_duplicates().head(8).tolist()
    return genes


def _rank_scores(table: pd.DataFrame, *, top_n: int) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame(columns=["group", "program", "mean_score", "n_cells", "theme"])
    required = {"group", "program", "mean_score"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Score table is missing required columns: {sorted(missing)}")
    ranked = table.copy()
    ranked["mean_score"] = pd.to_numeric(ranked["mean_score"], errors="coerce")
    ranked = ranked.dropna(subset=["mean_score"]).sort_values(["group", "mean_score"], ascending=[True, False])
    ranked = ranked.groupby("group", group_keys=False).head(top_n).reset_index(drop=True)
    ranked["theme"] = ranked["program"].astype(str).map(_program_theme)
    return ranked


def _markdown_list(table: pd.DataFrame, *, max_rows: int = 12) -> list[str]:
    lines: list[str] = []
    for row in table.head(max_rows).itertuples(index=False):
        n_cells = getattr(row, "n_cells", "NA")
        lines.append(
            f"- `{row.group}` is highest for `{row.program}` "
            f"(mean score {float(row.mean_score):.4g}, n={n_cells}); interpretation: {row.theme}."
        )
    return lines


def interpret_report(report_dir: str | Path, *, top_n: int = 5) -> dict[str, Any]:
    report_root = Path(report_dir).expanduser().resolve()
    tables_dir = report_root / "tables"
    manifest_path = report_root / "manifest.json"
    reference_status_path = report_root / "reference_status.json"

    grouped_scores = _read_table(tables_dir / "program_scores_by_group.tsv")
    neighbor_scores = _read_table(tables_dir / "neighbor_program_scores_by_group.tsv")
    reference_de = _read_table(tables_dir / "reference_de.tsv")
    if reference_de.empty:
        reference_de = _read_table(tables_dir / "reference_de_results.tsv.gz")

    manifest: dict[str, Any] = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    reference_status: dict[str, Any] = {}
    if reference_status_path.exists():
        reference_status = json.loads(reference_status_path.read_text(encoding="utf-8"))

    top_programs = _rank_scores(grouped_scores, top_n=top_n)
    top_neighbors = _rank_scores(neighbor_scores, top_n=top_n)
    genes_by_program = _top_reference_genes(reference_de)
    if not top_programs.empty:
        top_programs["top_reference_genes"] = top_programs["program"].astype(str).map(lambda program: ", ".join(genes_by_program.get(program, [])))
    if not top_neighbors.empty:
        top_neighbors["top_reference_genes"] = top_neighbors["program"].astype(str).map(lambda program: ", ".join(genes_by_program.get(program, [])))

    top_programs_path = tables_dir / "top_programs_by_roi_cell_type.tsv"
    top_neighbors_path = tables_dir / "top_neighbor_programs.tsv"
    top_programs.to_csv(top_programs_path, sep="\t", index=False)
    top_neighbors.to_csv(top_neighbors_path, sep="\t", index=False)

    summary = manifest.get("summary", {})
    config = manifest.get("config", {})
    lines = [
        "# Breast Xenium Reference Projection: Biological Interpretation",
        "",
        "## Working Conclusion",
        (
            "The reported scores measure transcriptional similarity between Xenium WTA cells and Perturb-seq-derived "
            "reference programs. A high score means a spatial cell state resembles that perturbation program; it does "
            "not prove that the tissue cell has the corresponding genetic perturbation or drug response."
        ),
        "",
        "## Run Context",
        f"- Report directory: `{report_root}`",
        f"- Spatial cells: `{summary.get('n_obs', 'NA')}`",
        f"- Spatial genes/features: `{summary.get('n_vars', 'NA')}`",
        f"- Reference program count: `{summary.get('program_count', 'NA')}`",
        f"- Reference datasets: `{', '.join(map(str, manifest.get('reference_datasets', []))) or 'NA'}`",
        "",
        "## Top Cell-Type/ROI Programs",
    ]
    lines.extend(_markdown_list(top_programs) or ["- No group-level program score table was available."])
    lines.extend(["", "## Top Neighbor Programs"])
    lines.extend(_markdown_list(top_neighbors) or ["- No neighbor program score table was available."])
    lines.extend(
        [
            "",
            "## Interpretation Guide",
            "- `GSE241115` programs are breast cancer CROP-seq references and are best used to nominate candidate regulators of tumor cell state, stem-like behavior, chromatin control, growth signaling, antigen presentation, and inflammatory-like programs.",
            "- `GSE281048` programs, when available, provide MCF7 pathway context for IFNB, IFNG, TGFB, TNFA, and INS responses.",
            "- ROI/cell-type enrichment should be interpreted as spatial localization of a reference-like state.",
            "- Neighbor enrichment should be interpreted as local microenvironment context around cells with similar program activity.",
            "",
            "## Caveats",
            "- The references are cell-line Perturb-seq datasets, while the query is FFPE breast tissue; cell-line biology may not fully match primary tissue.",
            "- Projection scores are association-style readouts, not causal perturbation evidence.",
            "- ROI and cell-group annotation quality directly affects biological interpretation.",
            "- Programs with broad stress, proliferation, interferon, or extracellular-matrix genes can reflect shared state biology rather than a single upstream regulator.",
        ]
    )
    if config.get("reference_effect_size_only"):
        lines.append(
            "- This run used effect-size-only reference DE for speed on full-scale data; program ranking is based on log2 fold-change, while p-values/FDR in `reference_de.tsv` should not be used for statistical claims."
        )
    if reference_status:
        blocked = {
            name: payload
            for name, payload in reference_status.items()
            if isinstance(payload, dict) and str(payload.get("status")) == "blocked"
        }
        if blocked:
            lines.extend(["", "## Blocked Optional References"])
            for name, payload in blocked.items():
                lines.append(f"- `{name}`: `{payload.get('reason', 'blocked')}`; {payload.get('message', '')}")

    markdown_path = report_root / "biological_interpretation.md"
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "markdown_path": markdown_path,
        "top_programs_path": top_programs_path,
        "top_neighbors_path": top_neighbors_path,
        "top_programs": top_programs,
        "top_neighbors": top_neighbors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", default="/data/taobo.hu/SpatialPerturb/reports/breast_reference_projection")
    parser.add_argument("--top-n", type=int, default=5)
    args = parser.parse_args()
    result = interpret_report(args.report_dir, top_n=args.top_n)
    print(json.dumps({key: str(value) for key, value in result.items() if key.endswith("_path")}, indent=2))


if __name__ == "__main__":
    main()
