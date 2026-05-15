#!/bin/bash
set -euo pipefail

# Stages this repository into ARTIFACT_ROOT/CBRP_original and instantiates Julia + CPLEX (IP).
# Intended to run inside the ip-builder container with /workspace = project root (read-only).

echo "=== Staging CBRP_original IP artifacts ==="

ARTIFACT_ROOT="${ARTIFACT_ROOT:-/artifacts/cbrp-original-ip}"
WORKSPACE="${WORKSPACE:-/workspace}"
CPLEX_ROOT_DIR="${CPLEX_ROOT_DIR:-/opt/ibm/ILOG/CPLEX_Studio2211}"
DEST="${ARTIFACT_ROOT}/CBRP_original"

echo "Artifact root: ${ARTIFACT_ROOT}"
echo "Workspace (source): ${WORKSPACE}"
echo "Destination: ${DEST}"
echo "CPLEX root: ${CPLEX_ROOT_DIR}"

mkdir -p "${ARTIFACT_ROOT}"

echo "Copying project into artifacts..."
rm -rf "${DEST}"
cp -a "${WORKSPACE}" "${DEST}"
mkdir -p "${DEST}/logs"
# Drop staged artifact tree from the copy (avoids nesting huge .julia trees on re-stage).
rm -rf "${DEST}/build/ip-artifacts" 2>/dev/null || true

echo "Setting up Julia environment..."
cd "${DEST}"

if [ ! -d "${CPLEX_ROOT_DIR}" ]; then
    echo "ERROR: CPLEX Studio root not found at ${CPLEX_ROOT_DIR}"
    echo "Set CPLEX_HOST_PATH (compose) or mount CPLEX to match CPLEX_ROOT_DIR."
    exit 1
fi

if [ -n "${CPLEX_STUDIO_BINARIES:-}" ]; then
    echo "Using preset CPLEX_STUDIO_BINARIES=${CPLEX_STUDIO_BINARIES}"
elif [ -d "${CPLEX_ROOT_DIR}/cplex/bin/x86-64_linux" ]; then
    export CPLEX_STUDIO_BINARIES="${CPLEX_ROOT_DIR}/cplex/bin/x86-64_linux"
    echo "Set CPLEX_STUDIO_BINARIES=${CPLEX_STUDIO_BINARIES}"
elif [ -d "${CPLEX_ROOT_DIR}/cplex/bin/arm64_linux" ]; then
    export CPLEX_STUDIO_BINARIES="${CPLEX_ROOT_DIR}/cplex/bin/arm64_linux"
    echo "Set CPLEX_STUDIO_BINARIES=${CPLEX_STUDIO_BINARIES}"
else
    echo "ERROR: No CPLEX binaries under ${CPLEX_ROOT_DIR}/cplex/bin/"
    echo "Expected x86-64_linux or arm64_linux. Is CPLEX_ROOT_DIR the Studio install root?"
    exit 1
fi

export JULIA_DEPOT_PATH="${ARTIFACT_ROOT}/.julia"
mkdir -p "${JULIA_DEPOT_PATH}"
echo "JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH}"

echo "Installing Julia packages (IP + CPLEX)..."
julia --threads=1 --project=. -e "
using Pkg
Pkg.instantiate()
Pkg.precompile()
"

echo "=== Staging complete ==="
echo "Artifacts at: ${ARTIFACT_ROOT}"
