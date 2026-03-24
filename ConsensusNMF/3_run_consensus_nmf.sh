#!/bin/bash
# 3_run_consensus_nmf.sh

# Step 3 of the ConsensusNMF pipeline.
# Runs a final NMF on the consensus co-assignment matrix produced in Step 2.
# This extracts stable factor matrices (U, V) and cluster assignments
# (U_assign, V_assign) from the aggregated signal across all initialisations.
#
# Can be run on the cpp output, the matlab output.
# Results are written back into the same consensus/ subdirectory as the
# consensus matrix itself.

#
# Usage:
#   bash 3_run_consensus_nmf.sh <cpp|matlab> <k>
#
# Arguments:
#   cpp|matlab   NMF implementation to use.
#                  cpp    — compiled C++ binary (faster; no MATLAB licence required)
#                  matlab — MATLAB implementation (requires MATLAB on PATH)
#   k            Number of NMF factors (latent clusters).
#
# Output:
#   example/example_out/<prefix>/cpp/consensus/   U.txt, V.txt, U_assign.txt, V_assign.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Arguments ─────────────────────────────────────────────────────────────────

CHOOSE=${1:?"ERROR: Missing required argument. Usage: bash 3_run_consensus_nmf.sh <cpp|matlab> <k>"}
k=${2:?"ERROR: Missing required argument. Usage: bash 3_run_consensus_nmf.sh <cpp|matlab> <k>"} 

if [[ "$CHOOSE" != "cpp" && "$CHOOSE" != "matlab" ]]; then
    echo "ERROR: First argument must be 'cpp' or 'matlab', got: '$CHOOSE'"
    exit 1
fi

# ── Configuration ─────────────────────────────────────────────────────────────
# Paths to executables / source
EXE="${SCRIPT_DIR}/../NMF/nmf_sparse_bp/run_nmf_bp"
MATLAB_DIR="${SCRIPT_DIR}/../NMF/"

# Input data
N_ROWS=4988        # number of cells (rows) in DATA
N_COLS=5000        # number of bins  (cols) in DATA

# NMF regularisation parameters
ALPHA=0.1
BETA=0.1
REG="regularized"  # regularisation mode: plain | regularized | sparse

# Output root — cpp and matlab each get their own subfolder under this prefix
PREFIX="chr6_k${k}_alpha_${ALPHA}_beta_${BETA}"
OUTDIR="${SCRIPT_DIR}/example/example_out/${PREFIX}/${CHOOSE}/consensus"
DATA=$OUTDIR/consensus_matrix.txt
# ── Preflight checks ──────────────────────────────────────────────────────────

echo "=== Step 3: NMF on consensus matrix ==="
echo "  Implementation : $CHOOSE"
echo "  k (factors)    : $k"
echo "  Dimensions     : ${N_ROWS} rows x ${N_COLS} cols"
echo ""

if [[ ! -f "$EXE" ]]; then
    echo "ERROR: C++ binary not found: $EXE"
    echo "  Compile it first — see NMF/nmf_sparse_bp/README.md"
    exit 1
fi

# ── Run NMF Implementation ─────────────────────────────────────────────────────
if [[ "$CHOOSE" == "cpp" ]]; then
    "$EXE" --data "$DATA"     \
           --nRows "$N_ROWS"  \
           --nCols "$N_COLS"  \
           -k      "$k"       \
           -o      "$OUTDIR"  \
           -R      "normal"     \
           -a      0   \
           -b      0   

elif [[ "$CHOOSE" == "matlab" ]]; then

        matlab -nodisplay -r \
            "addpath('${MATLAB_DIR}'); run_nmf('${DATA}', ${k}, '${OUTDIR}', 'ALPHA', 0, 'BETA', 0); exit;"

fi

echo ""
echo "=== Done. Consensus matrix and factorization written to: $OUTDIR"
echo ""
echo "    Next: run Step 2 (makeConsensusMatrix) across the I_* directories, then:"
echo "    bash run_consensus_nmf.sh $CHOOSE $k $N_ROWS $N_COLS"
