#!/usr/bin/env bash
# Carlos Path-CBRP MTZ (no SEC) over Alto Santo + Limoeiro.
# Usage (from external/CBRP_original):
#   bash run_carlos_path_mtz.sh                         # default groups
#   bash run_carlos_path_mtz.sh AS
#   bash run_carlos_path_mtz.sh AS --instance 2016
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/path_mtz_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus carlos \
  --batch-name path_cbrp_mtz.batch \
  --log-tag path_mtz \
  --solutions-dir solutions/path_cbrp_mtz \
  "$@"
