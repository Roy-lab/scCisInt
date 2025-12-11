EXE="../NMF/nmf_sparse_bp/run_nmf_bp"
DATA="example/example_in/pbmc_10x_transpose.txt"

for i in {1..10}
do
$EXE --data $DATA --nRows 2638 --nCols 6566 -k 20 -o example/example_out/pbmc_k20_alpha_0.1_beta_0.1/I_${i}/ -s 1 -r $i/
done 

