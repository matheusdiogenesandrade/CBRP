#!/usr/bin/env bash
# Shared engine for paper Path/Walk batch runners (Campinas + Carlos).
# Invoked by thin wrappers at the package root; not usually called directly.
#
# Required:
#   --corpus campinas|carlos
#   --batch-name FILE.batch   (under batchs/<group>/ or batchs/carlos/<group>/)
#   --log-tag TAG            (logs/${TAG}_<group>.log)
# Optional:
#   --solutions-dir DIR      (mkdir -p; default solutions/<batch-stem>)
#   --instance SPEC          (comma list / ranges; matched against instance stem)
#   positional group names   (defaults: all 8 Campinas or AS+LN)
#
# Env:
#   JULIA_BIN      Julia executable (default: julia on PATH)
#   JULIA_THREADS  --threads value (default: 1)
#
# Example:
#   JULIA_BIN=/opt/julia-1.10.5/bin/julia \
#     bash scripts/run_paper_batches.sh --corpus campinas \
#     --batch-name path_cbrp_mtz_simp.batch --log-tag path_mtz \
#     campinas-sparse --instance 1,3-5
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

JULIA_BIN="${JULIA_BIN:-julia}"
THREADS="${JULIA_THREADS:-1}"

if ! command -v "$JULIA_BIN" >/dev/null 2>&1 && [[ ! -x "$JULIA_BIN" ]]; then
  echo "Julia executable not found: ${JULIA_BIN} (set JULIA_BIN or PATH)" >&2
  exit 1
fi

CAMPINAS_GROUPS=(
  campinas-random
  campinas-sparse
  campinas-random-more-components
  campinas-sparse-more-components
  campinas-random-unitary-profits
  campinas-sparse-unitary-profits
  campinas-random-more-components-unitary-profits
  campinas-sparse-more-components-unitary-profits
)

CARLOS_GROUPS=(
  notified-alto-santo
  notified-limoeiro
)

resolve_carlos_group() {
  case "$1" in
    AS | as | alto-santo | notified-alto-santo) echo notified-alto-santo ;;
    LN | ln | limoeiro | notified-limoeiro) echo notified-limoeiro ;;
    *) echo "$1" ;;
  esac
}

instance_stem() {
  local path="$1"
  local base
  base="$(basename "$path")"
  base="${base%.sbrp}"
  base="${base%.txt}"
  printf '%s' "$base"
}

instance_matches() {
  local stem="$1"
  local id="$2"
  [[ "$stem" == "$id" ]] && return 0
  [[ "$stem" == *"-${id}" ]] && return 0
  return 1
}

CORPUS=""
BATCH_NAME=""
LOG_TAG=""
SOLUTIONS_DIR=""
INSTANCE_SPEC=""
BATCH_GROUPS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --corpus)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 1; }
      CORPUS="$2"
      shift 2
      ;;
    --corpus=*)
      CORPUS="${1#*=}"
      shift
      ;;
    --batch-name)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 1; }
      BATCH_NAME="$2"
      shift 2
      ;;
    --batch-name=*)
      BATCH_NAME="${1#*=}"
      shift
      ;;
    --log-tag)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 1; }
      LOG_TAG="$2"
      shift 2
      ;;
    --log-tag=*)
      LOG_TAG="${1#*=}"
      shift
      ;;
    --solutions-dir)
      [[ $# -ge 2 ]] || { echo "$1 requires a value" >&2; exit 1; }
      SOLUTIONS_DIR="$2"
      shift 2
      ;;
    --solutions-dir=*)
      SOLUTIONS_DIR="${1#*=}"
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

[[ -n "$CORPUS" ]] || { echo "--corpus is required (campinas|carlos)" >&2; exit 1; }
[[ -n "$BATCH_NAME" ]] || { echo "--batch-name is required" >&2; exit 1; }
[[ -n "$LOG_TAG" ]] || { echo "--log-tag is required" >&2; exit 1; }

case "$CORPUS" in
  campinas | carlos) ;;
  *)
    echo "Unknown --corpus: $CORPUS (expected campinas or carlos)" >&2
    exit 1
    ;;
esac

if [[ -z "$SOLUTIONS_DIR" ]]; then
  SOLUTIONS_DIR="solutions/${BATCH_NAME%.batch}"
fi

mkdir -p logs "$SOLUTIONS_DIR"

if [[ ${#BATCH_GROUPS[@]} -eq 0 ]]; then
  if [[ "$CORPUS" == campinas ]]; then
    BATCH_GROUPS=("${CAMPINAS_GROUPS[@]}")
  else
    BATCH_GROUPS=("${CARLOS_GROUPS[@]}")
  fi
fi

# Normalize Carlos aliases
if [[ "$CORPUS" == carlos ]]; then
  resolved=()
  for g in "${BATCH_GROUPS[@]}"; do
    resolved+=("$(resolve_carlos_group "$g")")
  done
  BATCH_GROUPS=("${resolved[@]}")
fi

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
  if [[ "$CORPUS" == campinas ]]; then
    batch="batchs/${g}/${BATCH_NAME}"
  else
    batch="batchs/carlos/${g}/${BATCH_NAME}"
  fi
  if [[ ! -f "$batch" ]]; then
    echo "Missing batch: $batch" >&2
    exit 1
  fi

  selected="$(mktemp "logs/.${LOG_TAG}_${g}.XXXXXX.sel")"
  if [[ ${#WANTED[@]} -gt 0 ]]; then
    declare -A LINE_OF=()
    declare -a STEMS=()
    while IFS= read -r line || [[ -n "$line" ]]; do
      [[ "$line" =~ ^[[:space:]]*(#|$) ]] && continue
      stem="$(instance_stem "${line%%[[:space:]]*}")"
      LINE_OF["$stem"]="$line"
      STEMS+=("$stem")
    done <"$batch"

    for id in "${WANTED[@]}"; do
      matched=""
      if [[ -n "${LINE_OF[$id]:-}" ]]; then
        matched="${LINE_OF[$id]}"
      else
        for stem in "${STEMS[@]}"; do
          if instance_matches "$stem" "$id"; then
            matched="${LINE_OF[$stem]}"
            break
          fi
        done
      fi
      if [[ -z "$matched" ]]; then
        echo "Instance ${id} not found in ${batch}" >&2
        rm -f "$selected"
        exit 1
      fi
      printf '%s\n' "$matched" >>"$selected"
    done
    unset LINE_OF
    log="logs/${LOG_TAG}_${g}.rerun-$(date +%Y%m%d-%H%M%S).log"
  else
    cp "$batch" "$selected"
    log="logs/${LOG_TAG}_${g}.log"
  fi

  n_runs="$(grep -cvE '^[[:space:]]*(#|$)' "$selected" || true)"
  echo "=== ${g} (${n_runs} instances) → ${log} ==="
  echo "    julia: ${JULIA_BIN} (threads=${THREADS})"
  "$JULIA_BIN" --threads="${THREADS}" --project=. src/run.jl --batch "$selected" \
    >"$log" 2>&1
  rm -f "$selected"
  echo "=== done ${g} ==="
done

echo "All requested ${CORPUS} batches (${BATCH_NAME}) finished."
