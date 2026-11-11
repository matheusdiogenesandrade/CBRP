#!/usr/bin/env bash
# Campinas Path-CBRP Lazy+SIMP (no MTZ, lazy SECs only) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_path_lazy.sh                         # default groups
#   bash run_campinas_path_lazy.sh campinas-sparse
#   bash run_campinas_path_lazy.sh campinas-sparse --instance 1,3-5
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/path_lazy_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus campinas \
  --batch-name path_cbrp_lazy_simp.batch \
  --log-tag path_lazy \
  --solutions-dir solutions/path_cbrp_lazy_simp \
  "$@"
