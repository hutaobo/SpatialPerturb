"""Package Nature Methods Brief Communication-facing report text."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

PATHWAY_REFERENCE = "gse281048_pathway_atlas"
BREAST_REFERENCE = "gse241115_breast_cropseq"
SHORTCOMM_CANDIDATE_PROGRAMS = (
    ("Mast Cells", f"{PATHWAY_REFERENCE}:FOS"),
    ("Basal-like Structured DCIS Cells", f"{PATHWAY_REFERENCE}:CEBPB"),
    ("Dendritic Cells", f"{PATHWAY_REFERENCE}:SP1"),
    ("Dendritic Cells", f"{PATHWAY_REFERENCE}:MTOR"),
    ("Dendritic Cells", f"{PATHWAY_REFERENCE}:RPS6KB1"),
    ("Dendritic Cells", f"{PATHWAY_REFERENCE}:MAPK3"),
    ("Luminal-like Amorphous DCIS Cells", f"{PATHWAY_REFERENCE}:PTGS2"),
    ("CAFs, Invasive Associated", f"{PATHWAY_REFERENCE}:MAPK8"),
    ("11q13 Invasive Tumor Cells (Mitotic)", f"{PATHWAY_REFERENCE}:IFNAR1"),
    ("11q13 Invasive Tumor Cells (Mitotic)", f"{PATHWAY_REFERENCE}:TYK2"),
)


def _read_rows(path: Path, *, limit: int = 12) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            rows.append(row)
            if len(rows) >= limit:
                break
    return rows


def _is_true(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _float_value(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _short_reference_label(value: str) -> str:
    if value == BREAST_REFERENCE:
        return "GSE241115"
    if value == PATHWAY_REFERENCE:
        return "GSE281048"
    return value


def _short_program_label(program: str) -> str:
    if ":" not in program:
        return program
    reference, gene = program.split(":", 1)
    return f"{_short_reference_label(reference)}:{gene}"


def _candidate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    ready = [
        row
        for row in rows
        if _is_true(row.get("is_claim_level", "true")) and row.get("claim_status", "claim_ready") == "claim_ready"
    ]
    if not ready:
        ready = rows
    selected = []
    for rank, (cell_type, program) in enumerate(SHORTCOMM_CANDIDATE_PROGRAMS):
        matches = [
            row
            for row in ready
            if row.get("program") == program
            and (row.get("cell_type") == cell_type or f"cell_type={cell_type}" in row.get("group", ""))
        ]
        if not matches:
            continue
        best = sorted(matches, key=lambda row: (_float_value(row.get("z_score", "0")), _float_value(row.get("mean_score", "0"))), reverse=True)[0]
        best = dict(best)
        best["candidate_rank"] = str(rank)
        best["program_label"] = _short_program_label(program)
        selected.append(best)
    if selected:
        return selected
    return sorted(ready, key=lambda row: _float_value(row.get("z_score", "0")), reverse=True)[:8]


def package_report(report_dir: str | Path, output_dir: str | Path | None = None) -> Path:
    root = Path(report_dir).expanduser().resolve()
    output_root = Path(output_dir).expanduser().resolve() if output_dir is not None else root
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    tables = root / "tables"
    claims = _read_rows(tables / "calibrated_program_scores_by_group.tsv", limit=2000)
    validation = _read_rows(tables / "reference_validation.tsv", limit=20)
    spatial = _read_rows(tables / "spatial_autocorrelation.tsv", limit=8)

    claim_rows = _candidate_rows(claims)
    min_fdr = min((_float_value(row.get("fdr"), 1.0) for row in claims), default=1.0)
    program_count = manifest.get("summary", {}).get("program_count", 268)

    lines = [
        "# Nature Methods Brief Communication Draft Scaffold",
        "",
        "## Working Title",
        "SpatialPerturb maps Perturb-seq programs onto unperturbed spatial transcriptomes",
        "",
        "## 70-word Abstract Draft",
        (
            "SpatialPerturb links perturbational single-cell references to unperturbed spatial transcriptomes. "
            "The toolkit derives Perturb-seq gene programs, calibrates spatial projection scores against matched nulls, "
            "and tests spatial organization on tissue graphs. Applied to Xenium WTA breast cancer, it nominates localized "
            "luminal/secretory invasive tumor-state programs while explicitly separating transcriptional similarity from "
            "causal perturbation. The workflow produces reproducible tables, figures and biological interpretation."
        ),
        "",
        "## Main Figure 1",
        "- Method schematic showing GSE241115 and GSE281048 as ready references, program expansion from 50 to 268, and GSE281048 contributing 218 MCF7 pathway programs.",
        "- Held-out reference AUROC, null calibration and ablation robustness are shown as credibility checks; candidate panels explicitly display the global FDR caveat.",
        "",
        "## Main Figure 2",
        "- Xenium breast application uses selected ranked candidate spatial programs: Mast-cell FOS, basal-like structured DCIS CEBPB, dendritic SP1/MTOR/RPS6KB1/MAPK3, luminal-like amorphous DCIS PTGS2, invasive-associated CAF MAPK8, and mitotic invasive tumor IFNAR1/TYK2.",
        "- These are presented as candidate Perturb-seq reference-like transcriptional states, not causal perturbation calls.",
        "",
        "## Figure Language Guardrail",
        f"- Use `candidate spatial programs`, `reference-like states` and `ranked calibrated projections`; avoid discovery language because the best calibrated global FDR is approximately `{min_fdr:.3g}`.",
        f"- Full source data retain all `{program_count}` programs, bootstrap confidence intervals, ablations, redundancy and FDR values for supplementary or Extended Data.",
        "",
        "## Figure Legend Drafts",
        (
            "Figure 1. SpatialPerturb projects perturbation programs from GSE241115 and GSE281048 onto Xenium WTA breast tissue. "
            "The validation panels show held-out recovery, null calibration and ablation robustness for 268 programs, including 218 MCF7 pathway programs from GSE281048. "
            f"The calibrated projections are ranked candidate programs because the best global FDR is approximately {min_fdr:.3g}."
        ),
        (
            "Figure 2. Selected ranked candidate spatial programs localize to interpretable breast tissue contexts, including Mast-cell FOS, basal-like structured DCIS CEBPB, "
            "dendritic SP1/MTOR/RPS6KB1/MAPK3, luminal-like amorphous DCIS PTGS2, invasive-associated CAF MAPK8 and mitotic invasive tumor IFNAR1/TYK2. "
            "These labels indicate Perturb-seq reference-like transcriptional states, not observed genetic perturbations or drug actions in the tissue."
        ),
        "",
        "## Key Results",
    ]
    if claim_rows:
        for row in claim_rows:
            lines.append(
                f"- `{row.get('group')}` resembles `{row.get('program_label', row.get('program'))}` "
                f"(z={_float_value(row.get('z_score')):.3g}, FDR={_float_value(row.get('fdr'), 1.0):.3g}, n={row.get('n_cells')})."
            )
    else:
        lines.append("- No ranked candidate rows were found in the calibrated table.")

    lines.extend(["", "## Reference Validation Snapshot"])
    if validation:
        for row in validation:
            lines.append(
                f"- `{row.get('reference_dataset', 'reference')}` `{row.get('program')}`: "
                f"AUROC={row.get('auroc', 'NA')}, AUPRC={row.get('auprc', 'NA')}, query coverage={row.get('query_gene_coverage', 'NA')}."
            )
    else:
        lines.append("- Reference validation table is unavailable.")

    lines.extend(["", "## Spatial Statistics Snapshot"])
    if spatial:
        for row in spatial:
            lines.append(f"- `{row.get('program')}`: Moran-style I={row.get('moran_i')}, FDR={row.get('fdr')}.")
    else:
        lines.append("- Spatial autocorrelation table is unavailable.")

    lines.extend(
        [
            "",
            "## Required Caveat Sentence",
            "Projection scores quantify Perturb-seq reference-like transcriptional states and do not prove that the spatial tissue contains the corresponding genetic perturbation, pathway intervention, or drug response.",
            "",
            "## Manifest",
            f"- Benchmark: `{manifest.get('benchmark', 'NA')}`",
            f"- References: `{', '.join(map(str, manifest.get('reference_datasets', [])))}`",
            f"- Figures: `{', '.join(manifest.get('figures', {}).keys())}`",
        ]
    )
    output = output_root / "nature_methods_shortcomm_scaffold.md"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", default="/data/taobo.hu/SpatialPerturb/reports/nature_methods_breast_shortcomm")
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()
    path = package_report(args.report_dir, output_dir=args.output_dir)
    print(path)


if __name__ == "__main__":
    main()
