#!/usr/bin/env bash
# Run Campinas Path-CBRP (simplified street digraph + MTZ only (no SEC)) over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_path_mtz_simp.sh              # all 8 groups sequentially
#   bash run_campinas_path_mtz_simp.sh campinas-random campinas-sparse
#   bash run_campinas_path_mtz_simp.sh campinas-random --instance 3,7,40-45
# --instance: ids N of data/<set>/N.sbrp (comma list, ranges a-b), run in the given order,
#   applied to every listed group; logs go to logs/path_mtz_simp_<group>.rerun-<timestamp>.log
#   so the full-batch log is not overwritten.
# Logs: logs/path_mtz_simp_<group>.log
# Summaries: grep -A1 '^instance,' logs/path_mtz_simp_*.log
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

# Path-CBRP-MTZ uses no SEC callback; keep env harmless if set
export PATH_CBRP_SEC_CALLBACK_LOG="${PATH_CBRP_SEC_CALLBACK_LOG:-0}"
THREADS="${JULIA_THREADS:-1}"
mkdir -p logs solutions/path_cbrp_mtz_simp

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
  batch="batchs/${g}/path_cbrp_mtz_simp.batch"
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi
  selected="$(mktemp "logs/.path_mtz_simp_${g}.XXXXXX.sel")"
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
    log="logs/path_mtz_simp_${g}.rerun-$(date +%Y%m%d-%H%M%S).log"
  else
    cp "$batch" "$selected"
    log="logs/path_mtz_simp_${g}.log"
  fi
  n_runs="$(grep -cvE '^[[:space:]]*(#|$)' "$selected" || true)"
  echo "=== ${g} (${n_runs} instances) → ${log} ==="
  julia --threads="${THREADS}" --project=. src/run.jl --batch "$selected" \
    >"$log" 2>&1
  rm -f "$selected"
  echo "=== done ${g} ==="
done

echo "All requested Campinas Path MTZ+simp batches finished."
