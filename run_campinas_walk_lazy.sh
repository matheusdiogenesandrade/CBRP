#!/usr/bin/env bash
# Campinas Walk-CBRP Lazy (no MTZ, lazy SECs) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_walk_lazy.sh                         # default groups
#   bash run_campinas_walk_lazy.sh campinas-sparse
#   bash run_campinas_walk_lazy.sh campinas-sparse --instance 1,3-5
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/walk_lazy_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus campinas \
  --batch-name walk_cbrp_lazy.batch \
  --log-tag walk_lazy \
  --solutions-dir solutions/walk_cbrp_lazy \
  "$@"
