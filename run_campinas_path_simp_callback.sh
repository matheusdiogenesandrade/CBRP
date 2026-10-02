#!/usr/bin/env bash
# Run Campinas Path-CBRP (simplified street digraph + SEC callback) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_path_simp_callback.sh              # all 8 groups sequentially
#   bash run_campinas_path_simp_callback.sh campinas-random campinas-sparse
#   bash run_campinas_path_simp_callback.sh campinas-random --instance 3,7,40-45
# --instance: ids N of data/<set>/N.sbrp (comma list, ranges a-b), run in the given order,
#   applied to every listed group; logs go to logs/path_simp_callback_<group>.rerun-<timestamp>.log
#   so the full-batch log is not overwritten.
# Logs: logs/path_simp_callback_<group>.log
# Summaries: grep -A1 '^instance,' logs/path_simp_callback_*.log
# Fractional SEC depth: SEC_USER_CUT_MAX_DEPTH (default 1000000 = every node; 0 = root only).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
THREADS="${JULIA_THREADS:-1}"
SEC_USER_CUT_MAX_DEPTH="${SEC_USER_CUT_MAX_DEPTH:-1000000}"
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
  batch="batchs/${g}/path_cbrp_simp_callback.batch"
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi
  selected="$(mktemp "logs/.path_simp_callback_${g}.XXXXXX.sel")"
  if [[ ${#WANTED[@]} -gt 0 ]]; then
    declare -A LINE_OF=()
    while IFS= read -r line || [[ -n "$line" ]]; do
      [[ "$line" =~ ^[[:space:]]*(#|$) ]] && continue
      LINE_OF["$(basename "${line%%[[:space:]]*}" .sbrp)"]="$line"
    done <"$batch"
    for id in "${WANTED[@]}"; do
      if [[ -z "${LINE_OF[$id]:-}" ]]; then
        echo "Instance ${id} not found in ${batch}" >&2
        rm -f "$selected"
        exit 1
      fi
      printf '%s\n' "${LINE_OF[$id]}" >>"$selected"
    done
    unset LINE_OF
    log="logs/path_simp_callback_${g}.rerun-$(date +%Y%m%d-%H%M%S).log"
  else
    cp "$batch" "$selected"
    log="logs/path_simp_callback_${g}.log"
  fi
  # run.jl parses each batch line on its own, so the flag must be on every line.
  run_batch="$(mktemp "logs/.path_simp_callback_${g}.XXXXXX.batch")"
  sed -E "/^[[:space:]]*(#|$)/! s/[[:space:]]*$/ --sec-user-cut-max-depth ${SEC_USER_CUT_MAX_DEPTH}/" \
    "$selected" >"$run_batch"
  rm -f "$selected"
  n_runs="$(grep -cvE '^[[:space:]]*(#|$)' "$run_batch" || true)"
  echo "=== ${g} (${n_runs} instances, sec-user-cut-max-depth=${SEC_USER_CUT_MAX_DEPTH}) → ${log} ==="
  julia --threads="${THREADS}" --project=. src/run.jl --batch "$run_batch" \
    >"$log" 2>&1
  rm -f "$run_batch"
  echo "=== done ${g} ==="
done

echo "All requested Campinas Path simp+callback batches finished."
