# ConsensusNMF

This module runs Non-negative Matrix Factorization (NMF) across multiple random initializations and aggregates the results into a single **consensus matrix**, improving the stability of cluster assignments before passing results to the scCisInt prediction module.

---

## How it works

NMF results vary depending on random initialization. ConsensusNMF addresses this by:

1. Running NMF `N` times from different random seeds, producing a cluster assignment per cell for each run.
2. Building a **consensus co-assignment matrix** — an cells × bins matrix where each entry counts how often two items were assigned to the same cluster across all runs. More co-assignments = more stable pairing.
3. Running one final NMF on the consensus matrix to extract stable cluster assignments.
4. Converting those assignments into a **pseudobulk accessibility matrix** (clusters × bins) for use by the scCisInt EP prediction module.

---

## Pipeline steps

```
Step 1  run_NMF_bp_multiple_init.sh
        Runs NMF N times (default 10), each from a different random seed.
        Output: example_out/<prefix>/I_1/ ... I_N/  (U.txt, V.txt, U_assign.txt, V_assign.txt)

Step 2  run_makeConsensusMatrix.sh  (requires MATLAB)
        Aggregates the N runs into a consensus co-assignment matrix.
        Output: example_out/<prefix>/consensus/consensus_matrix.txt

Step 3  run_consensus_nmf.sh
        Runs one final NMF on the consensus matrix.
        Output: example_out/<prefix>/consensus/U.txt, V.txt, U_assign.txt, V_assign.txt

Step 4  make_pseudobulk_matrix.R  (requires R)
        Converts the final U_assign.txt + original ATAC matrix into a
        clusters × bins pseudobulk mean accessibility matrix for scCisInt.
        Output: <pseudobulk_matrix>.txt
```

---

## Dependencies

| Tool | Used in |
|------|---------|
| C++11 + GSL (or MATLAB) | Step 1 — NMF runs |
| MATLAB | Step 2 — consensus matrix |
| R + `data.table` | Step 4 — pseudobulk matrix |

Build the C++ binary first if not using MATLAB (see [`../NMF/nmf_sparse_bp/README.md`](../NMF/nmf_sparse_bp/README.md)):

```bash
cd ../NMF/nmf_sparse_bp
mkdir -p build && cd build
cmake .. && make
mv run_nmf_bp ..
```

---

## Usage

Run all steps from inside the `ConsensusNMF/` directory.

### Step 1 — Multiple NMF runs

```bash
# Using the C++ binary (recommended if no MATLAB license)
bash run_NMF_bp_multiple_init.sh cpp <k>

# Using MATLAB
bash run_NMF_bp_multiple_init.sh matlab <k>
```

Replace `<k>` with your desired number of factors (e.g. `20`).

This runs NMF 10 times and writes results to `example/example_out/<prefix>/I_1/` through `I_10/`.

### Step 2 — Build consensus matrix

```bash
bash run_makeConsensusMatrix.sh
```

Requires MATLAB. Reads the `U_assign.txt` and `V_assign.txt` from each `I_*/` directory and writes `consensus_matrix.txt` to the `consensus/` subdirectory.

### Step 3 — NMF on consensus matrix

```bash
bash run_consensus_nmf.sh
```

Runs a final NMF on `consensus_matrix.txt` to produce stable factor matrices and cluster assignments in `consensus/`.

### Step 4 — Generate pseudobulk matrix

```bash
Rscript make_pseudobulk_matrix.R \
    <atac_matrix>    \   # original cells x bins matrix (comma-separated, no header)
    <bins_file>      \   # one bin name per line (e.g. example/example_in/chr6_bins.txt)
    <u_assign_file>  \   # consensus/U_assign.txt from Step 3
    <outfile>            # output path for the pseudobulk TSV
```

#### Example (chr6 data)

```bash
Rscript make_pseudobulk_matrix.R \
    example/example_in/chr6_subsample.txt \
    example/example_in/chr6_bins.txt \
    example/example_out/pbmc_k20_alpha_0.1_beta_0.1/consensus/U_assign.txt \
    pseudobulk_chr6_k20.txt
```

The output file is the `<pseudobulk_matrix>` argument for the scCisInt module.

---

## Output files

| File | Description |
|------|-------------|
| `I_*/U.txt` | Cell factor matrix from run `*` (cells × k) |
| `I_*/V.txt` | Bin factor matrix from run `*` (bins × k) |
| `I_*/U_assign.txt` | Cell cluster assignments from run `*` (one integer per line) |
| `I_*/V_assign.txt` | Bin cluster assignments from run `*` (one integer per line) |
| `consensus/consensus_matrix.txt` | Co-assignment matrix (cells × bins, tab-separated) |
| `consensus/consensus_matrix_sorted.png` | Heatmap of the consensus matrix sorted by cluster |
| `consensus/U.txt` | Final stable cell factor matrix |
| `consensus/U_assign.txt` | Final stable cell cluster assignments |
| `<outfile>.txt` | Pseudobulk mean accessibility matrix (clusters × bins) for scCisInt |

---

## Example data

| File | Description |
|------|-------------|
| `example/example_in/chr6_subsample.txt` | 4,988 cells × 5,000 bins (comma-separated, no header) |
| `example/example_in/chr6_bins.txt` | 5,000 bin names in `chr_start_end` format |
| `example/example_out/` | Pre-computed output (requires rerunning Step 1–3 to regenerate on your platform) |

> **Note:** The pre-computed `U_assign.txt` files in `example/example_out/` were generated on macOS and the binary must be recompiled for Linux/HPC before running. `make_pseudobulk_matrix.R` will print a clear error if it detects uninitialized assignments.

---

## File descriptions

| Script | Language | Role |
|--------|----------|------|
| `run_NMF_bp_multiple_init.sh` | Bash | Step 1 — run NMF N times |
| `makeConsensusMatrix.m` | MATLAB | Step 2 — build consensus matrix |
| `run_makeConsensusMatrix.sh` | Bash | Step 2 — wrapper to call MATLAB |
| `run_consensus_nmf.sh` | Bash | Step 3 — NMF on consensus matrix |
| `make_pseudobulk_matrix.R` | R | Step 4 — generate pseudobulk matrix for scCisInt |
