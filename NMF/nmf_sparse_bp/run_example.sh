export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/mnt/dv/wid/projects2/Roy-common/programs/thirdparty/gsl-2.6/lib

# NMF sparse block-pivot with alpha 0.1 beta 0.1
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --output output/sparse/

# NMF normal block-pivot
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --reg normal --output output/norm/

# NMF regularized block-pivot with alpha 0.1 beta 0.1
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --reg regularized --output output/reg/

# NMF with stop method pgrad
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --stop pgrad --output output/stop_meth_grad/

# NMF sparse block-pivot with alpha=1 beta=1
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --alpha 1 --beta 1 --output output/sparse_1_1/

# NMF sparse block-pivot with verbose
./run_nmf_bp --data input/toy/A.txt --nRows 95 --nCols 120 --nFactors 3 --output output/sparse/ --verbose
