#!/usr/bin/env bash
# Run Campinas Path-CBRP (simplified street digraph + SEC callback) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_path_simp_callback.sh              # all 8 groups sequentially
#   bash run_campinas_path_simp_callback.sh campinas-random campinas-sparse
# Logs: logs/path_simp_callback_<group>.log
# Summaries: grep -A1 '^instance,' logs/path_simp_callback_*.log
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
THREADS="${JULIA_THREADS:-1}"
mkdir -p logs solutions/path_cbrp_simp_callback

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
  batch="batchs/${g}/path_cbrp_simp_callback.batch"
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi
  log="logs/path_simp_callback_${g}.log"
  echo "=== ${g} ($(wc -l < "$batch") instances) → ${log} ==="
  julia --threads="${THREADS}" --project=. src/run.jl --batch "$batch" \
    >"$log" 2>&1
  echo "=== done ${g} ==="
done

echo "All requested Campinas Path simp+callback batches finished."
