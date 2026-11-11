#!/usr/bin/env bash
# Campinas Walk-CBRP VI (MTZ + intersection + user SECs) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_walk_vi.sh                         # default groups
#   bash run_campinas_walk_vi.sh campinas-sparse
#   bash run_campinas_walk_vi.sh campinas-sparse --instance 1,3-5
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/walk_vi_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus campinas \
  --batch-name walk_cbrp_vi.batch \
  --log-tag walk_vi \
  --solutions-dir solutions/walk_cbrp_vi \
  "$@"
