#!/bin/bash
# run_scATAC_EPpredict_stability_RF.sh
#
# Runs the full scCisInt EP-prediction pipeline for a single gene on one chromosome.
#
# Usage:
#   bash run_scATAC_EPpredict_stability_RF.sh <gene_name> <ntrees> <pseudobulk_matrix> <chr> [outdir]
#
# Arguments:
#   gene_name         Gene symbol (e.g. Sox5)
#   ntrees            Number of trees for the Random Forest (e.g. 500)
#   pseudobulk_matrix Path to the cluster-mean accessibility matrix (clusters x bins, TSV)
#   chr               Chromosome (e.g. chr6)
#   outdir            Output directory (optional; default: Results)
#
# Dependencies: R (with data.table, dplyr), MATLAB

set -euo pipefail

GENE=$1
NTREES=$2
NMFDATA=$3
CHR=$4
OUTDIR=${5:-Results}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"

###############################################################################
## Step 1: Prepare per-gene enhancer/promoter data
###############################################################################
echo "=== Step 1: Preparing data for gene: $GENE ==="

GENEBINFILE="${SCRIPT_DIR}/example/example_in/Mus_musculus.GRCm38.74.TSS_trans2gene2bin_morebins_${CHR}.txt"

Rscript --vanilla "${SCRIPT_DIR}/scCisInt_preprocessing_prep_data_by_gene.R" \
    "$GENE" "$NMFDATA" "$GENEBINFILE" 1000000 "$OUTDIR"

###############################################################################
## Step 2: Run Random Forest predictions for each region associated with the gene
###############################################################################
echo "=== Step 2: Running RF predictions ==="

# Each region produced by step 1 has a corresponding _enhancer_data.txt file.
# Loop over all such files and call scCisInt_EP_predict for each region.
shopt -s nullglob
ENHANCER_FILES=("${OUTDIR}"/*_enhancer_data.txt)

if [ ${#ENHANCER_FILES[@]} -eq 0 ]; then
    echo "ERROR: No enhancer data files found in ${OUTDIR}. Did step 1 succeed?"
    exit 1
fi

for enhfile in "${ENHANCER_FILES[@]}"; do
    region=$(basename "$enhfile" "_enhancer_data.txt")
    echo "  Predicting region: $region"
    matlab -batch \
        "addpath('${SCRIPT_DIR}'); scCisInt_EP_predict('${region}', ${NTREES}, '${OUTDIR}/', '${OUTDIR}/', true)"
done

###############################################################################
## Step 3: Map predictions to peaks of interest
###############################################################################
echo "=== Step 3: Merging predictions ==="

Rscript --vanilla "${SCRIPT_DIR}/scATAC_EPpredict_mergePred.R" "$GENE" "$CHR"

echo "=== Done. Results written to ${OUTDIR}/ ==="
