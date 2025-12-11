export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/mnt/dv/wid/projects2/Roy-common/programs/thirdparty/gsl-2.6/lib
DATA=/Volumes/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/datasets/liger/pbmc_alignment_srprocess/data/pbmc_10x_transpose.txt

mkdir output/pbmc/

./run_nmf_bp --data $DATA --nRows 2638 --nCols 6566 -k 8 -o output/pbmc/ -s 1  ## NMF sparse bp with alpha 0.1 beta 0.1
