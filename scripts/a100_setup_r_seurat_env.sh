#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/data/taobo.hu/SpatialPerturb}"
R_ENV_DIR="${R_ENV_DIR:-$DATA_ROOT/envs/r-seurat}"
REPORT_DIR="${REPORT_DIR:-$DATA_ROOT/reports/nature_methods_breast_shortcomm}"
STATUS_PATH="$REPORT_DIR/r_env_status.json"

mkdir -p "$DATA_ROOT/envs" "$REPORT_DIR"

manager=""
if command -v micromamba >/dev/null 2>&1; then
  manager="micromamba"
elif command -v mamba >/dev/null 2>&1; then
  manager="mamba"
elif command -v conda >/dev/null 2>&1; then
  manager="conda"
fi

if [[ -z "$manager" ]]; then
  python3 - <<PY
import json
from pathlib import Path
payload = {
    "status": "blocked",
    "reason": "CONDA_OR_MICROMAMBA_MISSING",
    "message": "Neither micromamba, mamba, nor conda is available on PATH.",
    "r_env_dir": "$R_ENV_DIR",
}
Path("$STATUS_PATH").write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
PY
  exit 0
fi

if [[ ! -x "$R_ENV_DIR/bin/Rscript" ]]; then
  if [[ "$manager" == "micromamba" ]]; then
    micromamba create -y -p "$R_ENV_DIR" -c conda-forge -c bioconda \
      r-base=4.3 r-seurat r-seuratobject r-matrix r-jsonlite
  elif [[ "$manager" == "mamba" ]]; then
    mamba create -y -p "$R_ENV_DIR" -c conda-forge -c bioconda \
      r-base=4.3 r-seurat r-seuratobject r-matrix r-jsonlite
  else
    conda create -y -p "$R_ENV_DIR" -c conda-forge -c bioconda \
      r-base=4.3 r-seurat r-seuratobject r-matrix r-jsonlite
  fi
fi

export PATH="$R_ENV_DIR/bin:$PATH"
if Rscript -e "suppressPackageStartupMessages({library(Seurat); library(Matrix)}); cat('R/Seurat OK\n')" >/tmp/spatialperturb_r_seurat_check.log 2>&1; then
  status="ready"
  reason=""
else
  status="blocked"
  reason="SEURAT_LOAD_FAILED"
fi
message="$(cat /tmp/spatialperturb_r_seurat_check.log)"

python3 - <<PY
import json
from pathlib import Path
payload = {
    "status": "$status",
    "reason": "$reason",
    "manager": "$manager",
    "r_env_dir": "$R_ENV_DIR",
    "rscript": "$R_ENV_DIR/bin/Rscript",
    "message": """$message""",
}
Path("$STATUS_PATH").write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
PY
