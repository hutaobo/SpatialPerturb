"""Package Nature Methods Brief Communication-facing report text."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


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


def package_report(report_dir: str | Path) -> Path:
    root = Path(report_dir).expanduser().resolve()
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    tables = root / "tables"
    claims = _read_rows(tables / "calibrated_program_scores_by_group.tsv", limit=20)
    validation = _read_rows(tables / "reference_validation.tsv", limit=8)
    spatial = _read_rows(tables / "spatial_autocorrelation.tsv", limit=8)

    claim_rows = [
        row
        for row in claims
        if row.get("is_claim_level", "").lower() in {"true", "1"} and row.get("claim_status") == "claim_ready"
    ]
    claim_rows = sorted(claim_rows, key=lambda row: float(row.get("z_score") or 0), reverse=True)[:8]

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
        "- Method schematic, held-out reference recovery, null-calibrated score distribution and spatial autocorrelation summary.",
        "",
        "## Main Figure 2",
        "- Xenium breast map, calibrated ROI/cell-type heatmap, spatial organization and top robust biological programs.",
        "",
        "## Key Results",
    ]
    if claim_rows:
        for row in claim_rows:
            lines.append(
                f"- `{row.get('group')}` resembles `{row.get('program')}` "
                f"(z={float(row.get('z_score') or 0):.3g}, FDR={float(row.get('fdr') or 1):.3g}, n={row.get('n_cells')})."
            )
    else:
        lines.append("- No claim-level rows were found in the calibrated table.")

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
            "Projection scores quantify transcriptional similarity to Perturb-seq-derived programs and do not prove that the spatial tissue contains the corresponding genetic perturbation or drug response.",
            "",
            "## Manifest",
            f"- Benchmark: `{manifest.get('benchmark', 'NA')}`",
            f"- References: `{', '.join(map(str, manifest.get('reference_datasets', [])))}`",
            f"- Figures: `{', '.join(manifest.get('figures', {}).keys())}`",
        ]
    )
    output = root / "nature_methods_shortcomm_scaffold.md"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", default="/data/taobo.hu/SpatialPerturb/reports/nature_methods_breast_shortcomm")
    args = parser.parse_args()
    path = package_report(args.report_dir)
    print(path)


if __name__ == "__main__":
    main()
