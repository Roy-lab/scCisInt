#!/bin/bash
# run_consensus_nmf.sh
#
# Step 3 of the ConsensusNMF pipeline.
# Runs a final NMF on the consensus co-assignment matrix produced in Step 2.
# This extracts stable factor matrices (U, V) and cluster assignments
# (U_assign, V_assign) from the aggregated signal across all initialisations.
#
# Can be run on the cpp output, the matlab output, or both in sequence.
# Results are written back into the same consensus/ subdirectory as the
# consensus matrix itself.
#
# Usage:
#   bash run_consensus_nmf.sh <cpp|matlab|both> <k> <nRows> <nCols>
#
# Arguments:
#   cpp|matlab|both   Which consensus matrix to run NMF on.
#   k                 Number of NMF factors (must match Steps 1 and 2).
#   nRows             Number of cells in the consensus matrix.
#   nCols             Number of bins  in the consensus matrix.
#
# Output layout (written into the existing consensus/ directory):
#   example/example_out/<prefix>/cpp/consensus/U.txt
#   example/example_out/<prefix>/cpp/consensus/V.txt
#   example/example_out/<prefix>/cpp/consensus/U_assign.txt
#   example/example_out/<prefix>/cpp/consensus/V_assign.txt
#   example/example_out/<prefix>/matlab/consensus/  (same, if matlab or both)
#
# Next step:
#   bash run_make_pseudobulk_matrix.sh <cpp|matlab|both> <k>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Arguments ─────────────────────────────────────────────────────────────────

CHOOSE=${1:?"ERROR: Missing required argument. Usage: bash run_consensus_nmf.sh <cpp|matlab|both> <k> <nRows> <nCols>"}
k=${2:?"ERROR: Missing required argument. Usage: bash run_consensus_nmf.sh <cpp|matlab|both> <k> <nRows> <nCols>"}
N_ROWS=${3:?"ERROR: Missing required argument. Usage: bash run_consensus_nmf.sh <cpp|matlab|both> <k> <nRows> <nCols>"}
N_COLS=${4:?"ERROR: Missing required argument. Usage: bash run_consensus_nmf.sh <cpp|matlab|both> <k> <nRows> <nCols>"}

if [[ "$CHOOSE" != "cpp" && "$CHOOSE" != "matlab" && "$CHOOSE" != "both" ]]; then
    echo "ERROR: First argument must be 'cpp', 'matlab', or 'both', got: '$CHOOSE'"
    exit 1
fi

# ── Configuration ─────────────────────────────────────────────────────────────

EXE="${SCRIPT_DIR}/../NMF/nmf_sparse_bp/run_nmf_bp"

ALPHA=0.1
BETA=0.1
PREFIX="chr6_k${k}_alpha_${ALPHA}_beta_${BETA}"

# ── Helper: run consensus NMF for one implementation ─────────────────────────
# Takes the implementation name (cpp or matlab) and runs NMF on its consensus
# matrix, writing stable factor outputs back into the same consensus/ directory.

run_consensus_for_impl() {
    local impl="$1"
    local consensus_dir="${SCRIPT_DIR}/example/example_out/${PREFIX}/${impl}/consensus"
    local consensus_matrix="${consensus_dir}/consensus_matrix.txt"

    echo "  [${impl}] Consensus matrix : $consensus_matrix"
    echo "  [${impl}] Output directory : $consensus_dir"

    if [[ ! -f "$consensus_matrix" ]]; then
        echo "  ERROR: Consensus matrix not found: $consensus_matrix"
        echo "    Has Step 2 (makeConsensusMatrix) completed for '$impl'?"
        return 1
    fi

    "$EXE" --data   "$consensus_matrix" \
           --nRows  "$N_ROWS"            \
           --nCols  "$N_COLS"            \
           -k       "$k"                 \
           -o       "$consensus_dir"     \
           --alpha  0                    \
           --beta   0                    \
           --reg    normal

    echo "  [${impl}] Done."
}

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

# ── Run ───────────────────────────────────────────────────────────────────────

if [[ "$CHOOSE" == "cpp" || "$CHOOSE" == "both" ]]; then
    run_consensus_for_impl "cpp"
    echo ""
fi

if [[ "$CHOOSE" == "matlab" || "$CHOOSE" == "both" ]]; then
    run_consensus_for_impl "matlab"
    echo ""
fi

echo "=== Done. Stable NMF factors written to consensus/ directories."
echo ""
echo "    Next: generate pseudobulk matrix for scCisInt:"
echo "    bash run_make_pseudobulk_matrix.sh $CHOOSE $k"
