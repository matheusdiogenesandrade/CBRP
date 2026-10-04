#!/usr/bin/env bash
# Run Campinas Complete-IP (metric-closure complete digraph + intersection cuts + SEC callback)
# over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_complete_callback.sh              # all 8 groups sequentially
#   bash run_campinas_complete_callback.sh campinas-random campinas-sparse
#   bash run_campinas_complete_callback.sh campinas-random --instance 3,7,40-45
# --instance: ids N of data/<set>/N.sbrp (comma list, ranges a-b), run in the given order,
#   applied to every listed group; logs go to logs/complete_callback_<group>.rerun-<timestamp>.log
#   so the full-batch log is not overwritten.
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

BATCH_GROUPS=()
INSTANCE_SPEC=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --instance | --instances)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 1; }
      INSTANCE_SPEC="$2"
      shift 2
      ;;
    --instance=* | --instances=*)
      INSTANCE_SPEC="${1#*=}"
      shift
      ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
    *)
      BATCH_GROUPS+=("$1")
      shift
      ;;
  esac
done
[[ ${#BATCH_GROUPS[@]} -gt 0 ]] || BATCH_GROUPS=("${ALL_BATCH_GROUPS[@]}")

WANTED=()
if [[ -n "$INSTANCE_SPEC" ]]; then
  IFS=',' read -ra parts <<<"${INSTANCE_SPEC//[[:space:]]/}"
  for p in "${parts[@]}"; do
    [[ -z "$p" ]] && continue
    if [[ "$p" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      lo=$((10#${BASH_REMATCH[1]}))
      hi=$((10#${BASH_REMATCH[2]}))
      ((lo <= hi)) || { echo "Bad range: $p" >&2; exit 1; }
      for ((k = lo; k <= hi; k++)); do WANTED+=("$k"); done
    else
      WANTED+=("$p")
    fi
  done
  [[ ${#WANTED[@]} -gt 0 ]] || { echo "--instance selected nothing" >&2; exit 1; }
fi

for g in "${BATCH_GROUPS[@]}"; do
  batch="batchs/${g}/complete_callback.batch"
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi
  run_batch="$(mktemp "logs/.complete_callback_${g}.XXXXXX.batch")"
  if [[ ${#WANTED[@]} -gt 0 ]]; then
    declare -A LINE_OF=()
    while IFS= read -r line || [[ -n "$line" ]]; do
      [[ "$line" =~ ^[[:space:]]*(#|$) ]] && continue
      LINE_OF["$(basename "${line%%[[:space:]]*}" .sbrp)"]="$line"
    done <"$batch"
    for id in "${WANTED[@]}"; do
      if [[ -z "${LINE_OF[$id]:-}" ]]; then
        echo "Instance ${id} not found in ${batch}" >&2
        rm -f "$run_batch"
        exit 1
      fi
      printf '%s\n' "${LINE_OF[$id]}" >>"$run_batch"
    done
    unset LINE_OF
    log="logs/complete_callback_${g}.rerun-$(date +%Y%m%d-%H%M%S).log"
  else
    cp "$batch" "$run_batch"
    log="logs/complete_callback_${g}.log"
  fi
  n_runs="$(grep -cvE '^[[:space:]]*(#|$)' "$run_batch" || true)"
  echo "=== ${g} (${n_runs} instances) → ${log} ==="
  julia --threads="${THREADS}" --project=. src/run.jl --batch "$run_batch" \
    >"$log" 2>&1
  rm -f "$run_batch"
  echo "=== done ${g} ==="
done

echo "All requested Campinas Complete-IP callback batches finished."
