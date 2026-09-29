#!/usr/bin/env bash
# Run Campinas Complete-IP (metric-closure complete digraph + intersection cuts + SEC callback)
# over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_complete_callback.sh              # all 8 groups sequentially
#   bash run_campinas_complete_callback.sh campinas-random campinas-sparse
# SEC threshold: --sec-min-violation default (1e-4). Callback stdout: COMPLETE_SEC_CALLBACK_LOG=1.
# Logs: logs/complete_callback_<group>.log
# Solutions: solutions/complete_callback/<group><instance>
# Summaries: grep -A1 '^instance,' logs/complete_callback_*.log
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
THREADS="${JULIA_THREADS:-1}"
mkdir -p logs solutions/complete_callback

# Do not name this GROUPS: bash reserves GROUPS for the user's Unix GIDs;
# assignments to GROUPS are ignored, so loops would iterate GIDs (e.g. 1009).
ALL_BATCH_GROUPS=(
  campinas-random
  campinas-sparse
  campinas-random-more-components
  campinas-sparse-more-components
  campinas-random-unitary-profits
  campinas-sparse-unitary-profits
  campinas-random-more-components-unitary-profits
  campinas-sparse-more-components-unitary-profits
)

if [[ $# -gt 0 ]]; then
  BATCH_GROUPS=("$@")
else
  BATCH_GROUPS=("${ALL_BATCH_GROUPS[@]}")
fi

for g in "${BATCH_GROUPS[@]}"; do
  batch="batchs/${g}/complete_callback.batch"
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi
  log="logs/complete_callback_${g}.log"
  echo "=== ${g} ($(wc -l < "$batch") instances) → ${log} ==="
  julia --threads="${THREADS}" --project=. src/run.jl --batch "$batch" \
    >"$log" 2>&1
  echo "=== done ${g} ==="
done

echo "All requested Campinas Complete-IP callback batches finished."
