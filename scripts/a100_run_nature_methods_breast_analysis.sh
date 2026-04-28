#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/data/taobo.hu/SpatialPerturb}"
REPO_DIR="${REPO_DIR:-$DATA_ROOT/code/SpatialPerturb}"
BRANCH="${BRANCH:-master}"
GITHUB_URL="${GITHUB_URL:-https://github.com/hutaobo/SpatialPerturb.git}"
ENV_DIR="${ENV_DIR:-$DATA_ROOT/envs/spatialperturb-py310}"
R_ENV_DIR="${R_ENV_DIR:-$DATA_ROOT/envs/r-seurat}"
REPORT_DIR="${REPORT_DIR:-$DATA_ROOT/reports/nature_methods_breast_shortcomm}"

mkdir -p "$DATA_ROOT/code" "$REPORT_DIR"

if [[ -d "$REPO_DIR/.git" ]]; then
  git -C "$REPO_DIR" fetch origin "$BRANCH"
  git -C "$REPO_DIR" checkout "$BRANCH"
  git -C "$REPO_DIR" pull --ff-only origin "$BRANCH"
else
  git clone --branch "$BRANCH" "$GITHUB_URL" "$REPO_DIR"
fi

REPORT_DIR="$REPORT_DIR" bash "$REPO_DIR/scripts/a100_setup_env.sh"
REPORT_DIR="$REPORT_DIR" R_ENV_DIR="$R_ENV_DIR" bash "$REPO_DIR/scripts/a100_setup_r_seurat_env.sh" || true

source "$ENV_DIR/bin/activate"
if [[ -x "$R_ENV_DIR/bin/Rscript" ]]; then
  export PATH="$R_ENV_DIR/bin:$PATH"
fi

python - <<PY
import json
import shutil
from pathlib import Path

checks = {
    "input_dir": Path("$DATA_ROOT/inputs/xenium_wta_breast").exists(),
    "cell_feature_matrix_h5": Path("$DATA_ROOT/inputs/xenium_wta_breast/cell_feature_matrix.h5").exists(),
    "cells_csv_gz": Path("$DATA_ROOT/inputs/xenium_wta_breast/cells.csv.gz").exists(),
    "rscript": shutil.which("Rscript"),
}
try:
    import spatialperturb
    import anndata
    import scanpy
    checks["python_imports"] = "ok"
except Exception as exc:
    checks["python_imports"] = f"failed: {exc}"
path = Path("$REPORT_DIR") / "a100_smoke.json"
path.write_text(json.dumps(checks, indent=2, default=str), encoding="utf-8")
print(json.dumps(checks, indent=2, default=str))
if checks["python_imports"] != "ok" or not checks["input_dir"] or not checks["cell_feature_matrix_h5"] or not checks["cells_csv_gz"]:
    raise SystemExit("A100 smoke checks failed.")
PY

python "$REPO_DIR/scripts/run_nature_methods_breast_analysis.py"
python "$REPO_DIR/scripts/a100_monitor_status.py" --report-dir "$REPORT_DIR" --session sp_nm_breast
