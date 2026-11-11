#!/usr/bin/env bash
# Run Campinas Walk-CBRP (complete digraph) paper variants over all 8 groups.
# Usage (from external/CBRP_original):
#   bash run_campinas_walk_cbrp.sh --variant vi
#   bash run_campinas_walk_cbrp.sh --variant lazy campinas-sparse
#   bash run_campinas_walk_cbrp.sh --variant cut campinas-random --instance 3,7,40-45
# --variant: lazy | cut | vi (required)
#   lazy = no MTZ, no intersection, complete-sec-mode lazy
#   cut  = no MTZ, no intersection, user+lazy SECs (user cuts root-only)
#   vi   = MTZ + intersection + user SECs root-only (no lazy callback context)
# --instance: ids N of data/<set>/N.sbrp (comma list, ranges a-b)
# Logs: logs/walk_cbrp_<variant>_<group>.log
# Summaries: grep -A1 '^instance,' logs/walk_cbrp_<variant>_*.log
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export COMPLETE_SEC_CALLBACK_LOG="${COMPLETE_SEC_CALLBACK_LOG:-0}"
THREADS="${JULIA_THREADS:-1}"

# Do not name this GROUPS: bash reserves GROUPS for the user's Unix GIDs.
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
VARIANT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --variant)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 1; }
      VARIANT="$2"
      shift 2
      ;;
    --variant=*)
      VARIANT="${1#*=}"
      shift
      ;;
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

[[ -n "$VARIANT" ]] || { echo "--variant is required (lazy|cut|vi)" >&2; exit 1; }
case "$VARIANT" in
  lazy | cut | vi) ;;
  *)
    echo "Unknown --variant: $VARIANT (expected lazy, cut, or vi)" >&2
    exit 1
    ;;
esac

[[ ${#BATCH_GROUPS[@]} -gt 0 ]] || BATCH_GROUPS=("${ALL_BATCH_GROUPS[@]}")

TAG="walk_cbrp_${VARIANT}"
BATCH_NAME="walk_cbrp_${VARIANT}.batch"
mkdir -p logs "solutions/${TAG}"

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
  batch="batchs/${g}/${BATCH_NAME}"
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi
  selected="$(mktemp "logs/.${TAG}_${g}.XXXXXX.sel")"
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
    log="logs/${TAG}_${g}.rerun-$(date +%Y%m%d-%H%M%S).log"
  else
    cp "$batch" "$selected"
    log="logs/${TAG}_${g}.log"
  fi
  n_runs="$(grep -cvE '^[[:space:]]*(#|$)' "$selected" || true)"
  echo "=== ${g} [${VARIANT}] (${n_runs} instances) → ${log} ==="
  julia --threads="${THREADS}" --project=. src/run.jl --batch "$selected" \
    >"$log" 2>&1
  rm -f "$selected"
  echo "=== done ${g} ==="
done

echo "All requested Campinas Walk-CBRP (${VARIANT}) batches finished."
