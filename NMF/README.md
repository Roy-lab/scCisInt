Non-Negative Matrix Factorization (NMF) Components

This subdirectory contains resources related to the non-negative matrix factorization (NMF) algorithms originally developed by Jingu Kim and collaborators. These implementations are provided for use in this project, and users are expected to cite the appropriate publications listed below.

References

If you use any part of this software, please cite:

[1] Core NMF Algorithm
Jingu Kim and Haesun Park.
Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons.
In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM'08), pp. 353–362, 2008.

If you use the nnls_solver = "as" option, please also cite:

[2] Alternating Non-Negativity Constrained Least Squares (Active Set Method)
Hyunsoo Kim and Haesun Park.
Nonnegative Matrix Factorization Based on Alternating Nonnegativity Constrained Least Squares and Active Set Method.
SIAM Journal on Matrix Analysis and Applications, 30, 713–730, 2008.

Included Implementations

The original block-pivot and active-set NMF algorithms were implemented in MATLAB. For users who do not have a MATLAB license, this project also provides a C++ implementation of the block-pivot algorithm.

If you use the C++ code, please cite the original publications listed above as well as the repository below.

C++ Implementation: nmf_sparse_bp

The C++ implementation of NMF_sparse_bp is also available as an independent GitHub repository:

https://github.com/Roy-lab/nmf_sparse_bp

Please refer to the README inside the nmf_sparse_bp subdirectory for build instructions, usage examples, and additional documentation.
