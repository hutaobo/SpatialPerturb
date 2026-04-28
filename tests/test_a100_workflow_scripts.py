from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


def _load_script_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_interpret_breast_reference_projection_writes_biology_outputs(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    module = _load_script_module(
        "interpret_breast_reference_projection",
        repo_root / "scripts" / "interpret_breast_reference_projection.py",
    )
    report_dir = tmp_path / "report"
    tables_dir = report_dir / "tables"
    tables_dir.mkdir(parents=True)
    (report_dir / "manifest.json").write_text(
        json.dumps(
            {
                "reference_datasets": ["gse241115_breast_cropseq"],
                "summary": {"n_obs": 10, "n_vars": 4, "program_count": 2},
            }
        ),
        encoding="utf-8",
    )
    (report_dir / "reference_status.json").write_text(
        json.dumps(
            {
                "gse281048_pathway_atlas": {
                    "status": "blocked",
                    "reason": "GSE281048_BLOCKED_RSCRIPT_MISSING",
                    "message": "Rscript is not available on PATH.",
                }
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {"group": "cell_type=Tumor | roi=Core", "program": "gse241115_breast_cropseq:STAT1", "mean_score": 3.2, "n_cells": 6},
            {"group": "cell_type=Stroma | roi=Edge", "program": "gse241115_breast_cropseq:TGFB1", "mean_score": 2.1, "n_cells": 4},
        ]
    ).to_csv(tables_dir / "program_scores_by_group.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"group": "cell_type=Tumor | roi=Core", "program": "gse241115_breast_cropseq:TGFB1", "mean_score": 1.5, "n_cells": 6},
        ]
    ).to_csv(tables_dir / "neighbor_program_scores_by_group.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"program": "gse241115_breast_cropseq:STAT1", "gene": "CXCL10", "log2fc": 2.0},
            {"program": "gse241115_breast_cropseq:TGFB1", "gene": "COL1A1", "log2fc": 1.7},
        ]
    ).to_csv(tables_dir / "reference_de.tsv", sep="\t", index=False)

    result = module.interpret_report(report_dir)

    assert result["markdown_path"].exists()
    assert result["top_programs_path"].exists()
    assert result["top_neighbors_path"].exists()
    text = result["markdown_path"].read_text(encoding="utf-8")
    assert "transcriptional similarity" in text
    assert "GSE281048_BLOCKED_RSCRIPT_MISSING" in text
    top_programs = pd.read_csv(result["top_programs_path"], sep="\t")
    assert "theme" in top_programs.columns


def test_a100_monitor_status_detects_complete_with_blocked_optional_reference(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    module = _load_script_module("a100_monitor_status", repo_root / "scripts" / "a100_monitor_status.py")
    report_dir = tmp_path / "report"
    tables_dir = report_dir / "tables"
    figures_dir = report_dir / "figures"
    tables_dir.mkdir(parents=True)
    figures_dir.mkdir()
    (report_dir / "manifest.json").write_text(json.dumps({"summary": {"program_count": 1}}), encoding="utf-8")
    (report_dir / "biological_interpretation.md").write_text("# Interpretation\n", encoding="utf-8")
    (report_dir / "reference_status.json").write_text(
        json.dumps({"gse281048_pathway_atlas": {"status": "blocked", "reason": "GSE281048_BLOCKED_RSCRIPT_MISSING"}}),
        encoding="utf-8",
    )
    (tables_dir / "program_scores_by_group.tsv").write_text("group\tprogram\tmean_score\n", encoding="utf-8")
    (tables_dir / "neighbor_program_scores_by_group.tsv").write_text("group\tprogram\tmean_score\n", encoding="utf-8")
    (figures_dir / "program_scores_heatmap.png").write_bytes(b"fake")

    status = module.build_status(report_dir, check_runtime=False)
    paths = module.write_status_files(status, report_dir)

    assert status["state"] == "COMPLETE_WITH_BLOCKED_OPTIONAL_REFERENCES"
    assert paths["json"].exists()
    assert paths["markdown"].exists()
