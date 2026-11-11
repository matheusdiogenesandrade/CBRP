#!/usr/bin/env bash
# Campinas Path-CBRP Cut+SIMP (no MTZ, user+lazy SECs) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_path_cut.sh                         # default groups
#   bash run_campinas_path_cut.sh campinas-sparse
#   bash run_campinas_path_cut.sh campinas-sparse --instance 1,3-5
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/path_cut_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus campinas \
  --batch-name path_cbrp_cut_simp.batch \
  --log-tag path_cut \
  --solutions-dir solutions/path_cbrp_cut_simp \
  "$@"
