#!/usr/bin/env Rscript
# make_pseudobulk_matrix.R
#
# Bridge script: converts Consensus NMF output into the cluster-mean
# accessibility matrix expected by the scCisInt module.
#
# NMF factorizes  X (cells x bins) ≈ U * V'.
# After consensus NMF, each cell is hard-assigned to the cluster with the
# highest loading in U (stored in U_assign.txt).  This script uses those
# assignments to compute the mean accessibility of each bin per cluster,
# producing the pseudobulk matrix that scCisInt_preprocessing_prep_data_by_gene.R
# expects as its <pseudobulk_matrix> argument.
#
# Usage:
#   Rscript make_pseudobulk_matrix.R \
#       <atac_matrix>   \  # cells x bins matrix (whitespace/comma-separated, no header)
#       <bins_file>     \  # one bin name per line (e.g. chr6_bins.txt)
#       <u_assign_file> \  # U_assign.txt from ConsensusNMF consensus/ output
#       <outfile>          # output TSV path (clusters x bins, with header)
#
# Output format (tab-separated, matches cluster_mean_matrix_K70_chr6.txt):
#   cluster  <bin1>  <bin2>  ...
#   1        0.12    0.00    ...
#   2        0.00    0.45    ...
#   ...

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: Rscript make_pseudobulk_matrix.R",
    "<atac_matrix> <bins_file> <u_assign_file> <outfile>"
  ), call. = FALSE)
}

atac_file   <- args[1]
bins_file   <- args[2]
assign_file <- args[3]
outfile     <- args[4]

# ---- Read inputs ------------------------------------------------------------

message("Reading ATAC-seq matrix: ", atac_file)
# Auto-detect separator (comma for the example data, whitespace otherwise)
first_line <- readLines(atac_file, n = 1)
sep <- if (grepl(",", first_line)) "," else ""
X <- fread(atac_file, header = FALSE, sep = sep)

message("  Matrix dimensions: ", nrow(X), " cells x ", ncol(X), " bins")

message("Reading bin names: ", bins_file)
bins <- readLines(bins_file)
bins <- bins[nzchar(trimws(bins))]   # drop blank lines

if (ncol(X) != length(bins)) {
  stop(sprintf(
    "Column count mismatch: matrix has %d columns but bins file has %d entries.",
    ncol(X), length(bins)
  ))
}
setnames(X, bins)

message("Reading cluster assignments: ", assign_file)
assignments <- scan(assign_file, what = integer(), quiet = TRUE)

if (nrow(X) != length(assignments)) {
  stop(sprintf(
    "Row count mismatch: matrix has %d rows (cells) but U_assign has %d entries.",
    nrow(X), length(assignments)
  ))
}

# Guard: INT_MIN (-2147483648) is the sentinel value written when NMF did not
# converge or the binary was compiled for a different platform and never ran.
INT_MIN <- -2147483648L
if (all(assignments == INT_MIN)) {
  stop(paste(
    "U_assign.txt contains only -2147483648 (INT_MIN) — the consensus NMF has",
    "not been run successfully on this machine yet.",
    "Please compile run_nmf_bp for your platform and re-run run_NMF_bp_multiple_init.sh",
    "followed by run_makeConsensusMatrix.sh before calling this script."
  ), call. = FALSE)
}
n_invalid <- sum(assignments == INT_MIN)
if (n_invalid > 0) {
  warning(sprintf(
    "%d of %d cells have invalid assignments (INT_MIN). They will be excluded.",
    n_invalid, length(assignments)
  ))
  keep <- assignments != INT_MIN
  X <- X[keep, ]
  assignments <- assignments[keep]
}

# ---- Compute cluster means --------------------------------------------------

clusters <- sort(unique(assignments))
message(sprintf("Computing pseudobulk means for %d clusters across %d bins...",
                length(clusters), ncol(X)))

result_list <- lapply(clusters, function(cl) {
  idx <- which(assignments == cl)
  if (length(idx) == 1) {
    as.numeric(X[idx, ])
  } else {
    colMeans(X[idx, ])
  }
})

result <- as.data.frame(do.call(rbind, result_list))
colnames(result) <- bins
result <- cbind(cluster = clusters, result)

# ---- Write output -----------------------------------------------------------

message("Writing pseudobulk matrix to: ", outfile)
fwrite(result, outfile, sep = "\t", quote = FALSE,
       row.names = FALSE, col.names = TRUE)

message(sprintf("Done. Output: %d clusters x %d bins (+ cluster column)",
                length(clusters), length(bins)))
