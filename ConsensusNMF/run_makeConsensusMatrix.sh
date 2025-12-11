#!/bin/bash

dir="example/example_out/pbmc_k20_alpha_0.1_beta_0.1"
U_file="U_assign.txt"
V_file="V_assign.txt"
num_iterations=10
outdir="$dir/consensus"

matlab -batch "makeConsensusMatrix('$dir', '$U_file', '$V_file', $num_iterations, '$outdir')"
