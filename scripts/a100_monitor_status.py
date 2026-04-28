"""Build status files for the A100 breast reference projection run."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any


def _run(command: list[str], *, timeout: int = 10) -> tuple[int, str]:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout)
    except FileNotFoundError:
        return 127, ""
    except subprocess.TimeoutExpired:
        return 124, ""
    return result.returncode, (result.stdout + result.stderr).strip()


def _tail(path: Path, *, lines: int = 40) -> list[str]:
    if not path.exists():
        return []
    return path.read_text(encoding="utf-8", errors="replace").splitlines()[-lines:]


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {"_error": f"Could not parse {path}"}


def _path_status(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() and path.is_file() else None,
    }


def build_status(
    report_dir: str | Path,
    *,
    session: str = "sp_breast_ref",
    check_runtime: bool = True,
) -> dict[str, Any]:
    report_root = Path(report_dir).expanduser().resolve()
    tables_dir = report_root / "tables"
    figures_dir = report_root / "figures"
    manifest_path = report_root / "manifest.json"
    interpretation_path = report_root / "biological_interpretation.md"
    run_log_path = report_root / "run.log"
    reference_status_path = report_root / "reference_status.json"

    tmux_running = False
    process_text = ""
    gpu_text = ""
    disk_text = ""
    if check_runtime:
        tmux_code, _ = _run(["tmux", "has-session", "-t", session])
        tmux_running = tmux_code == 0
        _, process_text = _run(["pgrep", "-af", "run_breast_reference_projection|spatialperturb|python"], timeout=5)
        _, gpu_text = _run(["nvidia-smi", "--query-gpu=name,memory.total,memory.used,utilization.gpu", "--format=csv,noheader"], timeout=10)
        _, disk_text = _run(["df", "-h", str(report_root)], timeout=5)

    expected_outputs = {
        "manifest": _path_status(manifest_path),
        "input_spatial": _path_status(report_root / "input_spatial.h5ad"),
        "program_scores_by_group": _path_status(tables_dir / "program_scores_by_group.tsv"),
        "neighbor_program_scores_by_group": _path_status(tables_dir / "neighbor_program_scores_by_group.tsv"),
        "program_scores_cell_level": _path_status(tables_dir / "program_scores_cell_level.tsv.gz"),
        "program_scores_heatmap": _path_status(figures_dir / "program_scores_heatmap.png"),
        "biological_interpretation": _path_status(interpretation_path),
    }
    if "nature_methods" in report_root.name:
        expected_outputs.update(
            {
                "nature_methods_summary": _path_status(report_root / "nature_methods_summary.md"),
                "main_figure_1": _path_status(figures_dir / "main_figure_1.png"),
                "main_figure_2": _path_status(figures_dir / "main_figure_2.png"),
                "reference_validation": _path_status(tables_dir / "reference_validation.tsv"),
                "calibrated_program_scores_by_group": _path_status(tables_dir / "calibrated_program_scores_by_group.tsv"),
                "spatial_autocorrelation": _path_status(tables_dir / "spatial_autocorrelation.tsv"),
                "ablation_summary": _path_status(tables_dir / "ablation_summary.tsv"),
            }
        )

    log_tail = _tail(run_log_path)
    log_text = "\n".join(log_tail)
    reference_status = _load_json(reference_status_path)
    if not reference_status and (report_root / "reference_prepare_status.json").exists():
        reference_status = _load_json(report_root / "reference_prepare_status.json")
    manifest = _load_json(manifest_path)
    if not reference_status and isinstance(manifest, dict):
        reference_status = manifest.get("reference_status", {}) or {}

    if manifest_path.exists() and interpretation_path.exists():
        state = "COMPLETE"
    elif "Traceback" in log_text or "Error" in log_text or "ERROR" in log_text:
        state = "FAILED" if not tmux_running else "RUNNING_WITH_ERRORS"
    elif tmux_running:
        state = "RUNNING"
    else:
        state = "PENDING"

    blocked_references = {
        name: payload
        for name, payload in reference_status.items()
        if isinstance(payload, dict) and str(payload.get("status")) == "blocked"
    }
    if state == "COMPLETE" and blocked_references:
        state = "COMPLETE_WITH_BLOCKED_OPTIONAL_REFERENCES"

    status = {
        "state": state,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "report_dir": str(report_root),
        "tmux_session": {"name": session, "running": tmux_running},
        "expected_outputs": expected_outputs,
        "reference_status": reference_status,
        "blocked_optional_references": blocked_references,
        "manifest_summary": manifest.get("summary", {}) if isinstance(manifest, dict) else {},
        "processes": process_text.splitlines() if process_text else [],
        "gpu": gpu_text.splitlines() if gpu_text else [],
        "disk": disk_text.splitlines() if disk_text else [],
        "log_tail": log_tail,
    }
    return status


def write_status_files(status: dict[str, Any], report_dir: str | Path) -> dict[str, Path]:
    report_root = Path(report_dir).expanduser().resolve()
    report_root.mkdir(parents=True, exist_ok=True)
    json_path = report_root / "status.json"
    md_path = report_root / "status.md"
    json_path.write_text(json.dumps(status, indent=2, default=str), encoding="utf-8")

    lines = [
        "# A100 Breast Reference Projection Status",
        "",
        f"- State: `{status['state']}`",
        f"- Checked at: `{status['checked_at']}`",
        f"- Report dir: `{status['report_dir']}`",
        f"- tmux `{status['tmux_session']['name']}` running: `{status['tmux_session']['running']}`",
        "",
        "## Expected Outputs",
    ]
    for name, payload in status["expected_outputs"].items():
        lines.append(f"- `{name}`: `{payload['exists']}` {payload['path']}")
    if status.get("blocked_optional_references"):
        lines.extend(["", "## Blocked Optional References"])
        for name, payload in status["blocked_optional_references"].items():
            lines.append(f"- `{name}`: `{payload.get('reason', 'blocked')}` {payload.get('message', '')}")
    if status.get("gpu"):
        lines.extend(["", "## GPU", "```text", *status["gpu"], "```"])
    if status.get("log_tail"):
        lines.extend(["", "## Log Tail", "```text", *status["log_tail"][-20:], "```"])
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"json": json_path, "markdown": md_path}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", default="/data/taobo.hu/SpatialPerturb/reports/breast_reference_projection")
    parser.add_argument("--session", default="sp_breast_ref")
    parser.add_argument("--no-runtime", action="store_true")
    args = parser.parse_args()
    status = build_status(args.report_dir, session=args.session, check_runtime=not args.no_runtime)
    paths = write_status_files(status, args.report_dir)
    print(json.dumps({"state": status["state"], "status_json": str(paths["json"]), "status_md": str(paths["markdown"])}, indent=2))


if __name__ == "__main__":
    main()
