?EXE="../NMF/nmf_sparse_bp/run_nmf_bp"
DATA="example/example_out/pbmc_k20_alpha_0.1_beta_0.1/consensus/consensus_matrix.txt"

$EXE --data $DATA --nRows 2638 --nCols 6566 -k 20 -o example/example_out/pbmc_k20_alpha_0.1_beta_0.1/consensus -s 1 --alpha 0 --beta 0
 

