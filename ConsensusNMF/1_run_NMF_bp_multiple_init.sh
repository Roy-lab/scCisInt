#!/bin/bash
# 1_run_NMF_bp_multiple_init.sh
#
# Step 1 of the ConsensusNMF pipeline.
# Runs NMF from N independent random initializations to sample the solution space.
# Each run uses a different random seed and writes to its own numbered subdirectory.
# The resulting U_assign.txt / V_assign.txt files are aggregated in Step 2
# (makeConsensusMatrix) to build a stable co-assignment consensus matrix.
#
# Usage:
#   bash run_NMF_bp_multiple_init.sh <cpp|matlab> <k> [n_init]
#
# Arguments:
#   cpp|matlab   NMF implementation to use.
#                  cpp    — compiled C++ binary (faster; no MATLAB licence required)
#                  matlab — MATLAB implementation (requires MATLAB on PATH)
#   k            Number of NMF factors (latent clusters).
#   n_init       Number of random initializations to run (default: 10).
#
# Output layout:
#   example/example_out/<prefix>/cpp/I_1/   U.txt, V.txt, U_assign.txt, V_assign.txt
#   example/example_out/<prefix>/cpp/I_2/   ...
#   ...
#   example/example_out/<prefix>/matlab/I_1/
#   ...
#
# Next step:
#   Run makeConsensusMatrix (Step 2) across the I_* directories, then:
#   bash run_consensus_nmf.sh <cpp|matlab|both> <k> <nRows> <nCols>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Arguments ─────────────────────────────────────────────────────────────────

CHOOSE=${1:?"ERROR: Missing required argument. Usage: bash run_NMF_bp_multiple_init.sh <cpp|matlab> <k> [n_init]"}
k=${2:?"ERROR: Missing required argument. Usage: bash run_NMF_bp_multiple_init.sh <cpp|matlab> <k> [n_init]"}
N_INIT=${3:-10}

if [[ "$CHOOSE" != "cpp" && "$CHOOSE" != "matlab" ]]; then
    echo "ERROR: First argument must be 'cpp' or 'matlab', got: '$CHOOSE'"
    exit 1
fi

# ── Configuration ─────────────────────────────────────────────────────────────
# Paths to executables / source
EXE="${SCRIPT_DIR}/../NMF/nmf_sparse_bp/run_nmf_bp"
MATLAB_DIR="${SCRIPT_DIR}/../NMF/"

# Input data
DATA="${SCRIPT_DIR}/example/example_in/chr6_subsample.txt"
N_ROWS=4988        # number of cells (rows) in DATA
N_COLS=5000        # number of bins  (cols) in DATA

# NMF regularisation parameters
ALPHA=0.1
BETA=0.1
REG="regularized"  # regularisation mode: plain | regularized | sparse

# Output root — cpp and matlab each get their own subfolder under this prefix
PREFIX="chr6_k${k}_alpha_${ALPHA}_beta_${BETA}"
OUTDIR_HEAD="${SCRIPT_DIR}/example/example_out/${PREFIX}/${CHOOSE}"

# ── Preflight checks ──────────────────────────────────────────────────────────

echo "=== Step 1: NMF multiple initializations ==="
echo "  Implementation   : $CHOOSE"
echo "  k (factors)      : $k"
echo "  N initializations: $N_INIT"
echo "  Input data       : $DATA  (${N_ROWS} x ${N_COLS})"
echo "  Output root      : $OUTDIR_HEAD"
echo ""

if [[ ! -f "$DATA" ]]; then
    echo "ERROR: Input data file not found: $DATA"
    exit 1
fi

if [[ "$CHOOSE" == "cpp" && ! -f "$EXE" ]]; then
    echo "ERROR: C++ binary not found: $EXE"
    echo "  Compile it first — see NMF/nmf_sparse_bp/README.md"
    exit 1
fi

mkdir -p "$OUTDIR_HEAD"

# ── Run N initializations ─────────────────────────────────────────────────────

if [[ "$CHOOSE" == "cpp" ]]; then

    for i in $(seq 1 "$N_INIT"); do
        outdir="${OUTDIR_HEAD}/I_${i}"
        echo "  [cpp] Init ${i}/${N_INIT}  seed=${i}  →  $outdir"
        "$EXE" --data "$DATA"     \
               --nRows "$N_ROWS"  \
               --nCols "$N_COLS"  \
               -k      "$k"       \
               -o      "$outdir"  \
               -R      "$REG"     \
               -a      "$ALPHA"   \
               -b      "$BETA"    \
               -r      "$i"
    done

elif [[ "$CHOOSE" == "matlab" ]]; then

    for i in $(seq 1 "$N_INIT"); do
        outdir="${OUTDIR_HEAD}/I_${i}"
	mkdir -p $outdir
        echo "  [matlab] Init ${i}/${N_INIT}  seed=${i}  →  $outdir"
        matlab -nodisplay -r \
            "addpath('${MATLAB_DIR}'); rng(${i}); run_nmf('${DATA}', ${k}, '${outdir}', 'ALPHA', $ALPHA, 'BETA', $BETA); exit;"
    done

fi

echo ""
echo "=== Done. ${N_INIT} initializations written to: $OUTDIR_HEAD"
echo ""
echo "    Next: run Step 2 (makeConsensusMatrix) across the I_* directories."
echo "    Wrapper: ./2_run_make_consensus_matrix.sh $CHOOSE $k $N_INIT
