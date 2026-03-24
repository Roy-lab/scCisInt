#!/bin/bash
# 2_run_make_consensus_matrix.sh
#
# Step 2 of the ConsensusNMF pipeline.
# Combines N independent NMF results with different random initilaizations. 
# Produces a N-by-M co-assignment matrix where cell and peaks which are assigned to same 
# clusters have non-zero value. The value is incremented per each observation co-assignment
# obvservation. The results of the method will be placed in a consensus directory. 
#
# Usage:
#   bash run_make_consensus_matrix.sh <cpp|matlab> <k>
#
# Arguments:
#   cpp|matlab   NMF implementation to use.
#                  cpp    — compiled C++ binary (faster; no MATLAB licence required)
#                  matlab — MATLAB implementation (requires MATLAB on PATH)
#   k            Number of NMF factors (latent clusters).
#   n_init       Number of random initialisations to run (default: 10).
#
# Output:
#  example/example_out/<prefix>/cpp/consensus 
#
# Next step:
#   bash run_consensus_nmf.sh <cpp|matlab|both> <k> 

set -euo pipefail 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Arguments ─────────────────────────────────────────────────────────────────

CHOOSE=${1:?"ERROR: Missing required argument. Usage: run_make_consensus_matrix.sh <cpp|matlab> <k> [n_init]"}
k=${2:?"ERROR: Missing required argument. Usage: bash run_make_consensus_matrix.sh <cpp|matlab> <k> [n_init]"}
N_INIT=${3:-10}

if [[ "$CHOOSE" != "cpp" && "$CHOOSE" != "matlab" ]]; then
    echo "ERROR: First argument must be 'cpp' or 'matlab', got: '$CHOOSE'"
    exit 1
fi

# ── Configuration ─────────────────────────────────────────────────────────────
# NMF regularisation parameters
ALPHA=0.1
BETA=0.1
REG="regularized"  # regularisation mode: plain | regularized | sparse

# Output root — cpp and matlab each get their own subfolder under this prefix
PREFIX="chr6_k${k}_alpha_${ALPHA}_beta_${BETA}"
OUTDIR_HEAD="${SCRIPT_DIR}/example/example_out/${PREFIX}/${CHOOSE}"

# ── Preflight checks ──────────────────────────────────────────────────────────

echo "=== Step 2: Make Consensus Matrix  ===" 
echo "  Implementation   : $CHOOSE"
echo "  k (factors)      : $k"
echo "  N initializations: $N_INIT"
echo "  Output root      : $OUTDIR_HEAD"

OUTDIR="$OUTDIR_HEAD/consensus"
if [[ -z "$OUTDIR" ]]; then
    echo "Error: OUTDIR is not set"
    exit 1
fi

if [[ ! -d "$OUTDIR" ]]; then
    mkdir -p "$OUTDIR"
fi


U_file=U_assign.txt
V_file=V_assign.txt


matlab -batch "makeConsensusMatrix('$OUTDIR_HEAD', '$U_file', '$V_file', $N_INIT, '$OUTDIR');"

echo ""
echo "=== Done. Consensus matrix written to $OUTDIR"
echo ""
echo "    Next: run Step 3 run_consensus_nmf.sh"
echo "    bash run_consensus_nmf.sh $CHOOSE $k"
