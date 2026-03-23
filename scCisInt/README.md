# scCisInt — EP Prediction Module

This directory contains the main scCisInt pipeline for predicting enhancer–promoter (EP) interactions from single-cell ATAC-seq data using Random Forest regression.

---

## Prerequisites

Before running this module you must have:

1. **A pseudobulk accessibility matrix** — a clusters × bins TSV where each row is a cluster and each column is a genomic bin. This is produced by the ConsensusNMF step:

   ```bash
   Rscript ../ConsensusNMF/make_pseudobulk_matrix.R \
       <atac_matrix> <bins_file> <u_assign_file> <pseudobulk_out.txt>
   ```

2. **A gene/bin annotation file** — maps each gene to its promoter genomic bin. The file expected follows the format of `example/example_in/Mus_musculus.GRCm38.74.TSS_trans2gene2bin_morebins_chr6.txt`.

3. **Software** — R (≥ 3.5) and MATLAB installed and on your `PATH`.

---

## Pipeline overview

```
run_scATAC_EPpredict_stability_RF.sh
 │
 ├─ Step 1  scCisInt_preprocessing_prep_data_by_gene.R
 │          Selects candidate enhancer bins within a genomic window of each
 │          gene's promoter and writes per-region input files to Results/.
 │
 ├─ Step 2  scCisInt_predict.m  (called via matlab -batch)
 │          Trains a 5-fold cross-validated Random Forest on each region and
 │          scores every candidate enhancer by variable importance.
 │
 └─ Step 3  scATAC_EPpredict_mergePred.R
            Maps raw predictions back to genomic peaks and integrates Hi-C
            contact data to produce the final scored EP interaction table.
```

---

## Running the pipeline

```bash
bash run_scATAC_EPpredict_stability_RF.sh \
    <gene_name>         \   # e.g. Sox5
    <ntrees>            \   # number of RF trees, e.g. 500
    <pseudobulk_matrix> \   # TSV from ConsensusNMF/make_pseudobulk_matrix.R
    <chr>               \   # chromosome, e.g. chr6
    [outdir]                # output directory (default: Results)
```

### Example (using provided chr6 data)

```bash
bash run_scATAC_EPpredict_stability_RF.sh \
    Sox5 \
    500 \
    example/example_in/cluster_mean_matrix_K70_chr6.txt \
    chr6 \
    Results
```

---

## Input files

| File | Description |
|------|-------------|
| `<pseudobulk_matrix>` | Clusters × bins accessibility matrix (TSV with header). Produced by `make_pseudobulk_matrix.R`. |
| `Mus_musculus.GRCm38.74.TSS_trans2gene2bin_morebins_<chr>.txt` | Gene-to-promoter-bin annotation. |
| `marker_peaks_ds_<chr>.txt` | Marker peaks for the chromosome (used in Step 3). |
| `<chr>_5kb_pairs_duan_scaled_1.txt` | Hi-C 5 kb contact pairs (used in Step 3). |

---

## Output files

After a successful run, `Results/` contains:

| File pattern | Description |
|---|---|
| `{region}_bins.txt` | Candidate bins within the search window of a promoter |
| `{region}_enhancer_data.txt` | Accessibility profiles for candidate enhancers (cells × enhancers) |
| `{region}_promoter_data.txt` | Accessibility profile for the promoter (cells × 1) |
| `{region}_consensus_EP_predictions.txt` | RF variable-importance scores for each enhancer–promoter pair |
| `{gene}_hicCount_confidence_predictions_beforemerging.txt` | Final merged predictions with Hi-C support |

---

## File descriptions

| Script | Language | Role |
|--------|----------|------|
| `run_scATAC_EPpredict_stability_RF.sh` | Bash | Master pipeline script |
| `scCisInt_preprocessing_prep_data_by_gene.R` | R | Step 1 — data preparation |
| `scCisInt_predict.m` | MATLAB | Step 2 — RF prediction per region |
| `scATAC_EPpredict_mergePred.R` | R | Step 3 — merge predictions with Hi-C |
| `scCisInt_postprocessing_combine_all_predictions.R` | R | Postprocessing — combine results across genes |
| `scCisInt_postprocessing_enhancer_scores.R` | R | Postprocessing — score enhancers |
| `scCisInt_postprocessing_kmeans.R` | R | Postprocessing — cluster enhancer programs |
| `scCisInt_postprocessing_integrate_aux_data.R` | R | Postprocessing — integrate auxiliary data |
| `scCisInt_aux_functions.R` | R | Shared helper functions |

---

## R package dependencies

```r
install.packages(c("data.table", "dplyr", "tidyr", "caret", "cluster"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```
