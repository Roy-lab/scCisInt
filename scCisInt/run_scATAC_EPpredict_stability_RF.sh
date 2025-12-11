#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
# Add these lines to run_foo.sh
tar xzf r2014b.tar.gz
mkdir cache
export MCR_CACHE_ROOT=`pwd`/cache

# untar your R installation. Make sure you are using the right version!
if [[ ! -f /usr/bin/R ]];then
echo "R not installed"
tar -xzf R351.tar.gz

# make sure the script will use your R installation, 
# and the working directory as its home location
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
fi
## R packages:
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf packages.tar.gz
export R_LIBS=$PWD/packages

genename=$2
nseed=$3
nmfdata=$4
chr=$5

###################################################################################
## Step 1: prepare data
echo Rscript --vanilla scATAC_EPpredict_prepData_forgenes.R $genename $nmfdata $chr
Rscript --vanilla scATAC_EPpredict_prepData_forgenes.R $genename $nmfdata $chr


###################################################################################
## Step 2: make predictions using RF
exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} \"${token}\"" 
      shift
  done
  echo "\"${exe_dir}/scATAC_EPpredict_stability_RF\"" $genename $nseed
  eval "\"${exe_dir}/scATAC_EPpredict_stability_RF\"" $genename $nseed
fi

###################################################################################
## Step 3: map predictions to peaks of interest
Rscript --vanilla scATAC_EPpredict_mergePred.R $genename $chr

#outpath=Results/
#bin1kb5kbfile=${outpath}/${genename}_bin1kbs_bin5kbs.txt
#cat $bin5kbfile | while read bin1kb bin5kb; 
#do
  #echo $bin1kb $bin5kb
  #predfile=${outpath}/${genename}_${bin1kb}_bins_prediction.txt
  #hicfile=${outpath}/${bin5kb}.txt
  #$bedtools/bedtools intersect -a $predfile -b $hicfile -wo | awk '$9>0' > $outpath/${genename}_${bin1kb}_enhancers_map.txt
#done


tar cvzf Results.tar.gz Results/
exit

