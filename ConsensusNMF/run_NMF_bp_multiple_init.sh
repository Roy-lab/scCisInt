EXE="../NMF/nmf_sparse_bp/run_nmf_bp"
MATLAB_DIR="../NMF/"
DATA="example/example_in/chr6_subsample.txt"

CHOOSE=$1
k=$2

outdir_head="example/example_out/chr6_k${k}_alpha_0.1_beta_0.1"

if [[ $CHOOSE == "cpp" ]]; then

    for i in {1..10}; do
        outdir="${outdir_head}/I_${i}"
        echo "$EXE --data $DATA --nRows 4988 --nCols 5000 -k $k -o $outdir -R regularized -s 1 -r $i"
	$EXE --data $DATA --nRows 4988 --nCols 5000 -k $k -o $outdir -R regularized -s 1 -r $i
    done 

elif [[ $CHOOSE == "matlab" ]]; then

    for i in {1..10}; do
        outdir="${outdir_head}/I_${i}"
        # Proper quoting for matlab command
        echo "matlab -nodisplay -r \"addpath('$MATLAB_DIR'); run_nmf('$DATA', $k, '$outdir'); exit;\""
	matlab -nodisplay -r "addpath('$MATLAB_DIR'); rng($i); run_nmf('$DATA', $k, '$outdir'); exit;"
    done

fi


