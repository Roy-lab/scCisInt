# scCisInt

**scCisInt** (single-cell Cis-regulatory Interaction) is a computational pipeline for predicting enhancer–promoter (EP) interactions from single-cell ATAC-seq data. It uses Non-negative Matrix Factorization (NMF) to identify latent chromatin accessibility patterns across cells, then trains a Random Forest model to predict which candidate enhancers regulate each gene.

## Overview

The pipeline consists of three main stages, each in its own subdirectory:

```
scCisInt/
├── NMF/                  # NMF algorithm implementations (MATLAB + C++)
├── ConsensusNMF/         # Runs NMF multiple times and builds a consensus matrix
└── scCisInt/             # Main EP prediction pipeline (R + MATLAB)
```

### Stage 1 — NMF (`NMF/`)

Factorizes a cells × genomic-bins matrix **X ≈ U · Vᵀ** using block-pivot NMF with optional sparsity regularization. Two implementations are provided:

- **MATLAB** (`nmf.m`, `run_nmf.m`) — requires MATLAB
- **C++** (`nmf_sparse_bp/`) — for users without a MATLAB license; compiled separately

See [`NMF/README.md`](NMF/README.md) and [`NMF/nmf_sparse_bp/README.md`](NMF/nmf_sparse_bp/README.md) for build instructions and usage.

### Stage 2 — Consensus NMF (`ConsensusNMF/`)

Runs NMF from multiple random initializations and aggregates the results into a single consensus matrix, improving stability. Scripts in this directory orchestrate the multi-run workflow using either the C++ or MATLAB NMF implementation.

### Stage 3 — EP Prediction (`scCisInt/`)

Uses the consensus NMF output (pseudobulk accessibility profiles per cluster) plus Hi-C interaction data to train a per-gene Random Forest model that scores candidate enhancer–promoter pairs. Includes preprocessing, prediction, and postprocessing scripts in R and MATLAB.

---

## Dependencies

| Tool | Version | Used in |
|------|---------|---------|
| R | ≥ 3.5 | scCisInt preprocessing & postprocessing |
| MATLAB | any recent | NMF (optional), consensus matrix, prediction |
| C++11 + CMake | — | C++ NMF build |
| [GSL](https://www.gnu.org/software/gsl/) | ≥ 2.6 | C++ NMF build |

**R packages:** `data.table`, `dplyr`, `tidyr`, `caret`, `cluster`, `GenomicRanges`

---

## Quick Start

### 1. Build the C++ NMF executable (if not using MATLAB)

```bash
cd NMF/nmf_sparse_bp
mkdir -p build && cd build
cmake ..
make
mv run_nmf_bp ..
```

### 2. Run Consensus NMF

```bash
cd ConsensusNMF

# Run NMF with 10 random initializations using the C++ binary
bash run_NMF_bp_multiple_init.sh cpp <k>

# Build the consensus matrix (requires MATLAB)
bash run_makeConsensusMatrix.sh

# Run a final NMF on the consensus matrix
bash run_consensus_nmf.sh
```

Replace `<k>` with your desired number of factors.

### 3. Run the EP prediction pipeline

```bash
cd scCisInt

# Preprocessing: prepare per-gene enhancer/promoter data
Rscript scCisInt_preprocessing_prep_data_by_gene.R \
    <gene_name> <pseudobulk_file> <gene_bin_file>

# Prediction: run Random Forest for each gene/region
bash run_scATAC_EPpredict_stability_RF.sh <MCR_root> <gene> <nseed> <nmfdata> <chr>

# Postprocessing: merge and score predictions
Rscript scCisInt_postprocessing_combine_all_predictions.R <indir> <outfile>
```

---

## Example Data

Each subdirectory contains an `example/` folder with small input datasets you can use to test the pipeline:

- `ConsensusNMF/example/example_in/` — subsampled PBMC scATAC-seq matrix
- `scCisInt/example/example_in/` — mouse chr6 ATAC-seq and Hi-C data

---

## Repository Structure

```
scCisInt/
│
├── NMF/
│   ├── nmf.m                      # MATLAB NMF entry point
│   ├── run_nmf.m                  # MATLAB wrapper
│   ├── nnlsm_blockpivot.m         # Block-pivot NNLS solver
│   ├── nnlsm_activeset.m          # Active-set NNLS solver
│   ├── solveNormalEqComb.m        # Normal equation solver
│   ├── README.md                  # NMF citations and notes
│   └── nmf_sparse_bp/             # C++ NMF implementation
│       ├── modules/               # Core C++ source files
│       ├── CMakeLists.txt
│       ├── run_nmf.cpp
│       └── README.md
│
├── ConsensusNMF/
│   ├── run_NMF_bp_multiple_init.sh   # Run NMF across multiple seeds
│   ├── run_makeConsensusMatrix.sh    # Build consensus matrix
│   ├── run_consensus_nmf.sh          # NMF on consensus matrix
│   ├── makeConsensusMatrix.m         # MATLAB consensus function
│   └── example/
│       ├── example_in/               # Input data
│       └── example_out/              # Pre-computed output
│
└── scCisInt/
    ├── scCisInt_preprocessing_prep_data_by_gene.R
    ├── scCisInt_preprocessing_prep_data_by_region.R
    ├── scCisInt_predict.m
    ├── scCisInt_postprocessing_combine_all_predictions.R
    ├── scCisInt_postprocessing_enhancer_scores.R
    ├── scCisInt_postprocessing_integrate_aux_data.R
    ├── scCisInt_postprocessing_kmeans.R
    ├── scCisInt_aux_functions.R
    ├── scATAC_EPpredict_mergePred.R
    ├── run_scATAC_EPpredict_stability_RF.sh
    └── example/
        └── example_in/               # Mouse chr6 example inputs
```

---

## Citation

If you use scCisInt, please also cite the NMF algorithm it depends on:

> Jingu Kim and Haesun Park. *Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons.* ICDM 2008, pp. 353–362.

If using the active-set solver (`nnls_solver = "as"`):

> Hyunsoo Kim and Haesun Park. *Nonnegative Matrix Factorization Based on Alternating Nonnegativity Constrained Least Squares and Active Set Method.* SIAM J. Matrix Anal. Appl., 30:713–730, 2008.

---

## License

See [LICENSE](LICENSE) for details.
