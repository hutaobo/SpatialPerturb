#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/data/taobo.hu/SpatialPerturb}"
REPO_DIR="${REPO_DIR:-$DATA_ROOT/code/SpatialPerturb}"
ENV_DIR="${ENV_DIR:-$DATA_ROOT/envs/spatialperturb-py310}"
REPORT_DIR="${REPORT_DIR:-$DATA_ROOT/reports/breast_reference_projection}"

mkdir -p "$DATA_ROOT/envs" "$REPORT_DIR"

if [[ ! -d "$ENV_DIR" ]]; then
  python3 -m venv "$ENV_DIR"
fi

source "$ENV_DIR/bin/activate"
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e "$REPO_DIR"

R_STATUS="missing"
R_MESSAGE="Rscript is not available on PATH."
if command -v Rscript >/dev/null 2>&1; then
  if Rscript -e "suppressPackageStartupMessages(library(Seurat)); cat('Seurat OK\n')" >/tmp/spatialperturb_r_check.log 2>&1; then
    R_STATUS="ready"
    R_MESSAGE="$(cat /tmp/spatialperturb_r_check.log)"
  else
    R_STATUS="blocked"
    R_MESSAGE="$(cat /tmp/spatialperturb_r_check.log)"
  fi
fi

python - <<PY
import json
from pathlib import Path
payload = {
    "env_dir": "$ENV_DIR",
    "repo_dir": "$REPO_DIR",
    "python": "$(command -v python)",
    "r_status": "$R_STATUS",
    "r_message": """$R_MESSAGE""",
}
path = Path("$REPORT_DIR") / "env_status.json"
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
PY
