#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/data/taobo.hu/SpatialPerturb}"
REPO_DIR="${REPO_DIR:-$DATA_ROOT/code/SpatialPerturb}"
ENV_DIR="${ENV_DIR:-$DATA_ROOT/envs/spatialperturb-py310}"
REPORT_DIR="${REPORT_DIR:-$DATA_ROOT/reports/breast_reference_projection}"
SESSION="${SESSION:-sp_breast_ref}"
INTERVAL="${INTERVAL:-60}"
WATCH=0

for arg in "$@"; do
  case "$arg" in
    --watch) WATCH=1 ;;
    --once) WATCH=0 ;;
    *) echo "Unknown argument: $arg" >&2; exit 2 ;;
  esac
done

PYTHON="python3"
if [[ -x "$ENV_DIR/bin/python" ]]; then
  PYTHON="$ENV_DIR/bin/python"
fi

run_once() {
  "$PYTHON" "$REPO_DIR/scripts/a100_monitor_status.py" --report-dir "$REPORT_DIR" --session "$SESSION"
}

if [[ "$WATCH" == "1" ]]; then
  while true; do
    run_once
    sleep "$INTERVAL"
  done
else
  run_once
fi
