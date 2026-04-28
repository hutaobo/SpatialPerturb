#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/data/taobo.hu/SpatialPerturb}"
REPO_DIR="${REPO_DIR:-$DATA_ROOT/code/SpatialPerturb}"
REPORT_DIR="${REPORT_DIR:-$DATA_ROOT/reports/nature_methods_breast_shortcomm}"
SESSION="${SESSION:-sp_nm_breast}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-300}"

watch_mode="false"
if [[ "${1:-}" == "--watch" ]]; then
  watch_mode="true"
fi

while true; do
  python "$REPO_DIR/scripts/a100_monitor_status.py" --report-dir "$REPORT_DIR" --session "$SESSION"
  if [[ "$watch_mode" != "true" ]]; then
    break
  fi
  sleep "$INTERVAL_SECONDS"
done
