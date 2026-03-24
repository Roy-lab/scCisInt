#!/bin/bash
# run_make_pseudobulk_matrix.sh
#
# Step 4 of the ConsensusNMF pipeline.
# Converts the stable cluster assignments from Step 3 into a pseudobulk
# cluster-mean accessibility matrix. This is the final output of ConsensusNMF
# and the direct input to the scCisInt EP-prediction module.
#
# For each cell, the cluster assignment from U_assign.txt is used to group
# cells, and the mean scATAC-seq accessibility across cells in each cluster
# is computed for every genomic bin, producing a clusters x bins matrix.
#
# Usage:
#   bash run_make_pseudobulk_matrix.sh <cpp|matlab|both> <k>
#
# Arguments:
#   cpp|matlab|both   Which consensus output to use (must match Steps 1-3).
#   k                 Number of NMF factors.
#
# Output:
#   example/example_out/<prefix>/cpp/pseudobulk_<prefix>.txt    (cpp or both)
#   example/example_out/<prefix>/matlab/pseudobulk_<prefix>.txt (matlab or both)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Arguments ─────────────────────────────────────────────────────────────────

CHOOSE=${1:?"ERROR: Missing required argument. Usage: bash run_make_pseudobulk_matrix.sh <cpp|matlab|both> <k>"}
k=${2:?"ERROR: Missing required argument. Usage: bash run_make_pseudobulk_matrix.sh <cpp|matlab|both> <k>"}

if [[ "$CHOOSE" != "cpp" && "$CHOOSE" != "matlab" && "$CHOOSE" != "both" ]]; then
    echo "ERROR: First argument must be 'cpp', 'matlab', or 'both', got: '$CHOOSE'"
    exit 1
fi

# ── Configuration ─────────────────────────────────────────────────────────────

ALPHA=0.1
BETA=0.1
PREFIX="chr6_k${k}_alpha_${ALPHA}_beta_${BETA}"

# Original scATAC-seq matrix (cells x bins, comma-separated, no header)
ATAC_MATRIX="${SCRIPT_DIR}/example/example_in/chr6_subsample.txt"

# One genomic bin name per line (chr_start_end format)
BINS_FILE="${SCRIPT_DIR}/example/example_in/chr6_bins.txt"

# ── Helper: generate pseudobulk matrix for one implementation ─────────────────
# Takes the implementation name (cpp or matlab), reads U_assign.txt from its
# consensus/ directory, and writes a clusters x bins TSV ready for scCisInt.

run_pseudobulk_for_impl() {
    local impl="$1"
    local u_assign="${SCRIPT_DIR}/example/example_out/${PREFIX}/${impl}/consensus/U_assign.txt"
    local outfile="${SCRIPT_DIR}/example/example_out/${PREFIX}/${impl}/pseudobulk_${PREFIX}.txt"

    echo "  [${impl}] U_assign file : $u_assign"
    echo "  [${impl}] Output file   : $outfile"

    for f in "$ATAC_MATRIX" "$BINS_FILE" "$u_assign"; do
        if [[ ! -f "$f" ]]; then
            echo "  ERROR: Required file not found: $f"
            echo "    Has Step 3 (run_consensus_nmf.sh) completed for '$impl'?"
            return 1
        fi
    done

    Rscript "${SCRIPT_DIR}/make_pseudobulk_matrix.R" \
        "$ATAC_MATRIX" \
        "$BINS_FILE"   \
        "$u_assign"    \
        "$outfile"

    echo "  [${impl}] Done."
    echo "  [${impl}] Pass to scCisInt with:"
    echo "    bash ../scCisInt/run_scATAC_EPpredict_stability_RF.sh <gene> <ntrees> $outfile <chr>"
}

# ── Preflight checks ──────────────────────────────────────────────────────────

echo "=== Step 4: Generating pseudobulk matrix ==="
echo "  Implementation : $CHOOSE"
echo "  k (factors)    : $k"
echo "  ATAC matrix    : $ATAC_MATRIX"
echo "  Bins file      : $BINS_FILE"
echo ""

for f in "$ATAC_MATRIX" "$BINS_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Shared input file not found: $f"
        exit 1
    fi
done

# ── Run ───────────────────────────────────────────────────────────────────────

if [[ "$CHOOSE" == "cpp" || "$CHOOSE" == "both" ]]; then
    run_pseudobulk_for_impl "cpp"
    echo ""
fi

if [[ "$CHOOSE" == "matlab" || "$CHOOSE" == "both" ]]; then
    run_pseudobulk_for_impl "matlab"
    echo ""
fi

echo "=== All done. ==="
