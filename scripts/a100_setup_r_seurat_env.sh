#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/data/taobo.hu/SpatialPerturb}"
R_ENV_DIR="${R_ENV_DIR:-$DATA_ROOT/envs/r-seurat}"
MICROMAMBA_DIR="${MICROMAMBA_DIR:-$DATA_ROOT/envs/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$DATA_ROOT/envs/micromamba-root}"
REPORT_DIR="${REPORT_DIR:-$DATA_ROOT/reports/nature_methods_breast_shortcomm}"
STATUS_PATH="$REPORT_DIR/r_env_status.json"

mkdir -p "$DATA_ROOT/envs" "$REPORT_DIR"

manager=""
if command -v micromamba >/dev/null 2>&1; then
  manager="$(command -v micromamba)"
elif command -v mamba >/dev/null 2>&1; then
  manager="$(command -v mamba)"
elif command -v conda >/dev/null 2>&1; then
  manager="$(command -v conda)"
fi

if [[ -z "$manager" ]]; then
  mkdir -p "$MICROMAMBA_DIR" "$MAMBA_ROOT_PREFIX"
  if [[ ! -x "$MICROMAMBA_DIR/bin/micromamba" ]]; then
    tmp_archive="$(mktemp --suffix=.tar.bz2)"
    curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest -o "$tmp_archive"
    tar -xjf "$tmp_archive" -C "$MICROMAMBA_DIR" bin/micromamba
    rm -f "$tmp_archive"
  fi
  manager="$MICROMAMBA_DIR/bin/micromamba"
fi

if [[ ! -x "$R_ENV_DIR/bin/Rscript" ]]; then
  manager_name="$(basename "$manager")"
  if [[ "$manager_name" == "micromamba" ]]; then
    "$manager" create -y -p "$R_ENV_DIR" -c conda-forge -c bioconda \
      r-base=4.3 r-seurat r-seuratobject r-matrix r-jsonlite
  elif [[ "$manager_name" == "mamba" ]]; then
    "$manager" create -y -p "$R_ENV_DIR" -c conda-forge -c bioconda \
      r-base=4.3 r-seurat r-seuratobject r-matrix r-jsonlite
  else
    "$manager" create -y -p "$R_ENV_DIR" -c conda-forge -c bioconda \
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
    "micromamba_dir": "$MICROMAMBA_DIR",
    "mamba_root_prefix": "$MAMBA_ROOT_PREFIX",
    "r_env_dir": "$R_ENV_DIR",
    "rscript": "$R_ENV_DIR/bin/Rscript",
    "message": """$message""",
}
Path("$STATUS_PATH").write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
PY
