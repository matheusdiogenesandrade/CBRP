#!/usr/bin/env bash
# Campinas Path-CBRP MTZ+SIMP (no SEC) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_path_mtz.sh                         # default groups
#   bash run_campinas_path_mtz.sh campinas-sparse
#   bash run_campinas_path_mtz.sh campinas-sparse --instance 1,3-5
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/path_mtz_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus campinas \
  --batch-name path_cbrp_mtz_simp.batch \
  --log-tag path_mtz \
  --solutions-dir solutions/path_cbrp_mtz_simp \
  "$@"
