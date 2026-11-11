#!/usr/bin/env bash
# Carlos Walk-CBRP Cut over Alto Santo + Limoeiro.
# Usage (from external/CBRP_original):
#   bash run_carlos_walk_cut.sh                         # default groups
#   bash run_carlos_walk_cut.sh AS
#   bash run_carlos_walk_cut.sh AS --instance 2016
# --instance: Campinas = N of data/<set>/N.sbrp; Carlos = stem / year suffix
# Logs: logs/walk_cut_<group>.log
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
exec bash "$ROOT/scripts/run_paper_batches.sh" \
  --corpus carlos \
  --batch-name walk_cbrp_cut.batch \
  --log-tag walk_cut \
  --solutions-dir solutions/walk_cbrp_cut \
  "$@"
