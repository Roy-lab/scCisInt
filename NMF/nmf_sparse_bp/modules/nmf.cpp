#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <cmath>
#include <iostream>
#include <numeric>
#include "initialization.h"
#include "nnlsm_blockpivot.h"
#include "nmf.h"
#include "utils.h"
#include "io.h"
#include "stdio.h"

//SH has hadded functionality from the alternating nonnegative constrained least squares
// using block principal pivoting/active Set method

//The modification are generated using the script
//mnt/dv/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/nmf_blockpivot/nmf_bpas/nmf.m and
//mnt/dv/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/nmf_blockpivot/nmf_bpas/nnlsm_blockpivot.m

//If utilized use site:
// Jingu Kim and Haesun Park, Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons,
//      In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM'08), 353-362, 200


// Using formulation:
// (1) minimize 1/2 * || A-WH ||_F^2
// (2) minimize 1/2 * ( || A-WH ||_F^2 + alpha * || W ||_F^2 + beta * || H ||_F^2 )
// (3) minimize 1/2 * ( || A-WH ||_F^2 + alpha * || W ||_F^2 + beta * (sum_(i=1)^n || H(:,i) ||_1^2 ) )


NMF::NMF(int k, reg_type regType, method_type methodType, stop_rule_type stopType, int maxIter, int minIter, bool verb,
         double termTol, double alpha_in, double beta_in)
{
        n_components = k;
        max_iter = maxIter;
        min_iter = minIter;
        type = regType;
        method = methodType;
        stop_rule = stopType;
        verbose = verb;
        tol = termTol;
        alpha = alpha_in;
        beta = beta_in;
        list<double> reconstruction_err_;
}

NMF::~NMF()
{
};

int
NMF::appendRows(gsl_matrix *append, gsl_matrix *Block1, gsl_matrix *Block2)
{
        gsl_matrix_view append_view1 = gsl_matrix_submatrix(append, 0, 0, Block1->size1, Block1->size2);
        gsl_matrix_memcpy(&append_view1.matrix, Block1);
        gsl_matrix_view append_view2 = gsl_matrix_submatrix(append, Block1->size1, 0, Block2->size1, Block2->size2);
        gsl_matrix_memcpy(&append_view2.matrix, Block2);
        return 0;
}

int
NMF::appendRows(gsl_matrix *append, gsl_matrix *Block1, gsl_vector *Block2)
{
        gsl_matrix_view append_view1 = gsl_matrix_submatrix(append, 0, 0, Block1->size1, Block1->size2);
        gsl_matrix_memcpy(&append_view1.matrix, Block1);
        gsl_vector_view append_view2 = gsl_matrix_row(append, Block1->size1);
        gsl_vector_memcpy(&append_view2.vector, Block2);
        return 0;
}


int
NMF::getGradient()
{
        int n = A->size1;
        int m = A->size2;
        //Compute the base of all three gradiants. This is the plain version.
        // Compute gradW: (original code) W (H * H') - A * H'

        gsl_matrix *HHt = gsl_matrix_calloc(n_components, n_components);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, H, H, 0, HHt);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, A, H, 0, gradW);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, W, HHt, -1, gradW);

        gsl_matrix_free(HHt);

        //Compute gradH: (W' * W)* H - W' * A
        gsl_matrix *WtW = gsl_matrix_calloc(n_components, n_components);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, W, W, 0, WtW);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, W, A, 0, gradH);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, WtW, H, -1, gradH);

        gsl_matrix_free(WtW);

        //No additional process needed
        if (type == regularized)
        {
                //gradU_tot = GradW + alpha * W = gradU + alpha * U
                gsl_matrix_scale(W, alpha);
                gsl_matrix_add(gradW, W);
                gsl_matrix_scale(W, 1.0 / alpha);


                //gradV_tot = Grad
                gsl_matrix_scale(H, beta);
                gsl_matrix_add(gradH, H);
                gsl_matrix_scale(H, 1.0 / beta);

        } else if (type == sparse)
        {
                //gradU_tot = GradW + alpha * W = gradU + alpha * U
                gsl_matrix_scale(W, alpha);
                gsl_matrix_add(gradW, W);
                gsl_matrix_scale(W, 1.0 / alpha);

                //gradV_tot = GradV to  beta * [ones matrix] * H
                gsl_matrix *ones = gsl_matrix_calloc(n_components, n_components);
                gsl_matrix_set_all(ones, 1.0);
                gsl_matrix_scale(ones, beta);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ones, H, 1, gradH);
                gsl_matrix_free(ones);
        }
	
	
	for (int i = 0; i < gradW->size1; i++)
        {
        	for (int j = 0; j < gradW->size2; j++)
                {
                	if (abs(gsl_matrix_get(gradW, i, j)) < 1e-12)
                        { // For numerical stability
                        	gsl_matrix_set(gradW, i, j, 0);
                        }
                }
        }

	for (int i = 0; i < gradH->size1; i++)
        {
                for (int j = 0; j < gradH->size2; j++)
                {
                        if (abs(gsl_matrix_get(gradH, i, j)) < 1e-12)
                        { // For numerical stability
                                gsl_matrix_set(gradH, i, j, 0);
                        }
                }
        }

        return 0;
}

int
NMF::nnlsm(gsl_matrix *A_mat, gsl_matrix *B_mat, gsl_matrix *X, gsl_matrix *Y)
{
        if (method == block_pivot)
        {
                nnlsm_bp::nnlsmBlockPivot(A_mat, B_mat, X, Y);
        } else if (method == active_set)
        {
                //nlsm_as::nnlsmActiveSet(A, B,X, Y);
                cout << "Active set method is not implemented. Coming soon." << endl;
        }
}

int
NMF::fit(gsl_matrix *inputmat, gsl_matrix *U, gsl_matrix *V)
{
        A = inputmat;
        W = U;
        H = V;
        gsl_matrix *At = gsl_matrix_alloc(A->size2, A->size1);
        gsl_matrix_transpose_memcpy(At, A);
        gsl_matrix *Wt = gsl_matrix_alloc(W->size2, W->size1);
        gsl_matrix *Ht = gsl_matrix_alloc(H->size2, H->size1);
        int m = A->size1;
        int n = A->size2;
        int k = W->size2;

        if (W->size2 != H->size1)
        {
                cout << "Inner dimension of W and H must be the same." << endl;
                return 1;
        }

        if (W->size2 != n_components || H->size1 != n_components)
        {
                cout << "Inner dimension of W and H must be equal to k" << endl;
                return 1;
        }

        if (W->size1 != m || H->size2 != n)
        {
                cout << "outer dimensions of W (m x k) and H (k * n) must be the same as A( m * n)." << endl;
                return 1;
        }

        if (verbose)
        {
                cout << "Dimension of input: " << m << " x " << n << "." << endl;
                cout << "Lower dimension embedding: " << k << '.' << endl;
        }

        gradW = gsl_matrix_calloc(W->size1, W->size2); //Same size
        gsl_matrix *gradWt = gsl_matrix_calloc(W->size2, W->size1);
        gradH = gsl_matrix_calloc(H->size1, H->size2);

        getGradient();
        io::write_dense_matrix("tmp/gradW.txt", gradW);
        io::write_dense_matrix("tmp/gradH.txt", gradH);
        io::write_dense_matrix("tmp/W.txt", W);
        io::write_dense_matrix("tmp/H.txt", H);
        io::write_dense_matrix("tmp/data.txt", A);

        double initSC = getStartValue();
        double SC = numeric_limits<double>::infinity();
        int SCconv = 0;
        int SC_COUNT = 3;

        //Make additional append matrices for adding reg rows.
        if (type == normal)
        {
                for (int n_iter = 0; n_iter < max_iter; n_iter++)
                {
                        //Run to Optimize
                        nnlsm(W, A, H, gradH);
                        if (utils::has_nan(H, "H iter " + to_string(n_iter))) { return -1; }
                        gsl_matrix_transpose_memcpy(Wt, W);
                        gsl_matrix_transpose_memcpy(Ht, H);
                        nnlsm(Ht, At, Wt, gradWt);
                        gsl_matrix_transpose_memcpy(W, Wt);
                        gsl_matrix_transpose_memcpy(gradW, gradWt);
                        if (utils::has_nan(W, "W iter " + to_string(n_iter))) { return -1; }

                        //Check conditions
                        getGradient();
                        SC = getStopValue();


                        if (n_iter > min_iter)
                        {
                                if (SC / initSC <= tol)
                                {
                                        SCconv = SCconv + 1;
                                        if (SCconv >= SC_COUNT)
                                        {
                                                break;
                                        }
                                } else
                                {
                                        SCconv = 0;
                                }
                        }

                        if (verbose)
                        {
                                char buffer[75];
                                sprintf(buffer,
                                        "iter: %i\t stop value: %.02e\t relative stop value: %.02e\t # stop achieved: %i",
                                        n_iter, SC, SC / initSC, SCconv);
                                cout << buffer << endl;
                        }
                }
        } else if (type == regularized)
        {
                gsl_matrix *salphaI = gsl_matrix_calloc(n_components, n_components);
                gsl_matrix_add_diagonal(salphaI, sqrt(alpha));
                gsl_matrix *sbetaI = gsl_matrix_calloc(n_components, n_components);
                gsl_matrix_add_diagonal(sbetaI, sqrt(beta));
                gsl_matrix *zerokn = gsl_matrix_calloc(n_components, n);
                gsl_matrix *zerokm = gsl_matrix_calloc(n_components, m);
                gsl_matrix *W_append = gsl_matrix_alloc(W->size1 + n_components, W->size2);
                gsl_matrix *Ht_append = gsl_matrix_alloc(H->size2 + n_components, H->size1);
                gsl_matrix *A_append = gsl_matrix_alloc(A->size1 + n_components, A->size2);
                gsl_matrix *At_append = gsl_matrix_alloc(A->size2 + n_components, A->size1);
                appendRows(A_append, A, zerokn);
                appendRows(At_append, At, zerokm);

                bool nan_found = false;
                for (int n_iter = 0; n_iter < max_iter; n_iter++)
                {
                        appendRows(W_append, W, sbetaI);
                        nnlsm(W_append, A_append, H, gradH);
                        if (utils::has_nan(H, "H iter " + to_string(n_iter))) { nan_found = true; break; }

                        gsl_matrix_transpose_memcpy(Wt, W);
                        gsl_matrix_transpose_memcpy(Ht, H);

                        appendRows(Ht_append, Ht, salphaI);
                        nnlsm(Ht_append, At_append, Wt, gradWt);
                        gsl_matrix_transpose_memcpy(W, Wt);
                        gsl_matrix_transpose_memcpy(gradW, gradWt);
                        if (utils::has_nan(W, "W iter " + to_string(n_iter))) { nan_found = true; break; }

                        //Check conditions
                        getGradient();
                        SC = getStopValue();

                        if (n_iter > min_iter)
                        {
                                if (SC / initSC <= tol)
                                {
                                        SCconv = SCconv + 1;
                                        if (SCconv >= SC_COUNT)
                                        {
                                                break;
                                        }
                                } else
                                {
                                        SCconv = 0;
                                }
                        }
                        if (verbose)
                        {
                                char buffer[75];
                                sprintf(buffer,
                                        "iter: %i\t stop value: %.02e\t relative stop value: %.02e\t # stop achieved: %i",
                                        n_iter, SC, SC / initSC, SCconv);
                                cout << buffer << endl;
                        }
                }


                gsl_matrix_free(salphaI);
                gsl_matrix_free(sbetaI);
                gsl_matrix_free(zerokn);
                gsl_matrix_free(zerokm);
                gsl_matrix_free(W_append);
                gsl_matrix_free(Ht_append);
                gsl_matrix_free(A_append);
                gsl_matrix_free(At_append);
                if (nan_found) { return -1; }

        } else if (type == sparse)
        {
                gsl_matrix *salphaI = gsl_matrix_calloc(n_components, n_components);
                gsl_matrix_add_diagonal(salphaI, sqrt(alpha));
                gsl_vector *sbetaE = gsl_vector_calloc(n_components);
                gsl_vector_set_all(sbetaE, sqrt(beta));
                gsl_vector *zero1n = gsl_vector_calloc(n);
                gsl_matrix *zerokm = gsl_matrix_calloc(n_components, m);
                gsl_matrix *W_append = gsl_matrix_alloc(W->size1 + 1, W->size2);
                gsl_matrix *Ht_append = gsl_matrix_alloc(H->size2 + n_components, H->size1);
                gsl_matrix *A_append = gsl_matrix_alloc(A->size1 + 1, A->size2);
                gsl_matrix *At_append = gsl_matrix_alloc(A->size2 + n_components, A->size1);

                appendRows(A_append, A, zero1n);
                appendRows(At_append, At, zerokm);


                bool nan_found = false;
                for (int n_iter = 0; n_iter < max_iter; n_iter++)
                {
                        appendRows(W_append, W, sbetaE);
                        nnlsm(W_append, A_append, H, gradH);
                        if (utils::has_nan(H, "H iter " + to_string(n_iter))) { nan_found = true; break; }
                        //io::write_dense_matrix("tmp/H.txt", H);

                        gsl_matrix_transpose_memcpy(Wt, W);
                        gsl_matrix_transpose_memcpy(Ht, H);

                        appendRows(Ht_append, Ht, salphaI);
                        nnlsm(Ht_append, At_append, Wt, gradWt);
                        gsl_matrix_transpose_memcpy(W, Wt);
                        gsl_matrix_transpose_memcpy(gradW, gradWt);
                        if (utils::has_nan(W, "W iter " + to_string(n_iter))) { nan_found = true; break; }
                        //io::write_dense_matrix("tmp/W.txt", W);

                        //Check conditions
                        getGradient();
                        SC = getStopValue();
			//io::write_dense_matrix("tmp/gradH.txt", gradH);
			//io::write_dense_matrix("tmp/gradW.txt", gradW);
                        if (n_iter > min_iter)
                        {
                                if (SC / initSC <= tol)
                                {
                                        SCconv = SCconv + 1;
                                        if (SCconv >= SC_COUNT)
                                        {
                                                break;
                                        }
                                } else
                                {
                                        SCconv = 0;
                                }
                        }
                        if (verbose)
                        {
                                char buffer[75];
                                sprintf(buffer,
                                        "iter: %i\t stop value: %.02e\t relative stop value: %.02e\t # stop achieved: %i",
                                        n_iter, SC, SC / initSC, SCconv);
                                cout << buffer << endl;
                        }
                }

                gsl_matrix_free(salphaI);
                gsl_vector_free(sbetaE);
                gsl_vector_free(zero1n);
                gsl_matrix_free(zerokm);
                gsl_matrix_free(W_append);
                gsl_matrix_free(Ht_append);
                gsl_matrix_free(A_append);
                gsl_matrix_free(At_append);
                if (nan_found) { return -1; }


        } else
        {
                cout << "Incorrect type provided. Use type normal, regularized, or sparse." << endl;
                return -999;
        }

	rescaleVectors();	
        gsl_matrix_free(gradW);
        gsl_matrix_free(gradH);


        return 0;
}

double
NMF::getStartValue()
{
        int m = W->size1;
        int k = W->size2;
        int n = H->size2;
        gsl_matrix *gradHt = gsl_matrix_alloc(n, k);
        gsl_matrix_transpose_memcpy(gradHt, gradH);
        gsl_matrix *gradMerge = gsl_matrix_alloc(m + n, k);
        appendRows(gradMerge, gradW, gradHt);
        double retVal = sqrt(utils::get_frobenius_norm(gradMerge));

        int numAll = (m * k) + (n * k);
        if (stop_rule == normalized_pgrad)
        {
                return retVal / numAll;
        } else if (stop_rule == pgrad)
        {
                return retVal;
        } else
        {
                cout << "Unsupported stop rule." << endl;
                return -999;
        }
}

double
NMF::getStopValue()
{
        /*getGradient();
        gsl_vector *grad_vector = gsl_vector_alloc(W->size1 * W->size2 + H->size1 * H->size2);
        gsl_vector_view W_vector = gsl_vector_subvector(grad_vector, 0, W->size1 * W->size2);
        gsl_vector_view H_vector = gsl_vector_subvector(grad_vector, W->size1 * W->size2, H->size1 * H->size2);
        gslMatrixRowsFlatVect(&W_vector.vector, W);
        gslMatrixRowsFlatVect(&H_vector.vector, H);
        double norm = gsl_blas_dnrm2(grad_vector);
        */

	double ij_val;
	double ij_grad;
        list<double> pGrad;
        for (int i = 0; i < W->size1; i++)
        {
                for (int j = 0; j < W->size2; j++)
                {
			ij_val = gsl_matrix_get(W, i, j);
			ij_grad = gsl_matrix_get(gradW, i, j);
                        if (ij_val > 0 ||ij_grad < 0)
                        {
                                pGrad.push_back(ij_grad);
                        }
                }

        }

        for (int i = 0; i < H->size1; i++)
        {
                for (int j = 0; j < H->size2; j++)
                {
			ij_val = gsl_matrix_get(H, i, j);
			ij_grad = gsl_matrix_get(gradH, i, j);
                        if (ij_val > 0 || ij_grad < 0)
                        {
                                pGrad.push_back(ij_grad);
                        }
                }

        }

        double norm = sqrt(inner_product(pGrad.begin(), pGrad.end(), pGrad.begin(), 0));

        if (stop_rule == normalized_pgrad)
        {
                int nonzero_length = pGrad.size();
                if (nonzero_length > 0)
                {
                        norm = norm / nonzero_length;
                } else
                {
                        norm = 0; //for stability
                }

        } else if (stop_rule == pgrad)
        {
                //Nothing to do here.
        } else
        {
                cout << "Unselected stop rule." << endl;
                return -999;
        }
        return norm;
}

int
NMF::gslMatrixRowsFlatVect(gsl_vector *Vector, gsl_matrix *Matrix)
{
        for (int i = 0; i < Matrix->size1; i++)
        {
                gsl_vector_view i_row = gsl_matrix_row(Matrix, i);
                gsl_vector_view sub_vector = gsl_vector_subvector(Vector, i * Matrix->size2, Matrix->size2);
                gsl_vector_memcpy(&sub_vector.vector, &i_row.vector);
        }

        return 0;
}

int
NMF::gslMatrixColFlatVect(gsl_vector *Vector, gsl_matrix *Matrix)
{
        for (int j = 0; j < Matrix->size2; j++)
        {
                gsl_vector_view i_row = gsl_matrix_column(Matrix, j);
                gsl_vector_view sub_vector = gsl_vector_subvector(Vector, j * Matrix->size1, Matrix->size1);
                gsl_vector_memcpy(&sub_vector.vector, &i_row.vector);
        }

        return 0;
}

int
NMF::setSize(int NSamp, int NFeat)
{
        nSamples = NSamp;
        nFeatures = NFeat;
        return 0;
}

int 
NMF::rescaleVectors() //Rescale W to unit norm 
{

	int k = W->size2;


        for(int i = 0; i < k; i++)
	{
		gsl_matrix_view W_col_mat = gsl_matrix_submatrix(W, 0, i, W->size1, 1);
		if (utils::is_all_zeros(&W_col_mat.matrix, "W column " + to_string(i)))
		{
			continue;
		}

		gsl_vector_view W_i_view = gsl_matrix_column(W, i);
		gsl_vector_view H_i_view = gsl_matrix_row(H, i);
		double norm = gsl_blas_dnrm2(&W_i_view.vector);
		gsl_vector_scale(&W_i_view.vector, 1/norm);
		gsl_vector_scale(&H_i_view.vector, norm);
	}
	return 0;
}


int
NMF::assignClusters()
{
        row_clusters.clear();
        column_clusters.clear();
        for (int i = 0; i < nSamples; i++)
        {
                int max_idx = std::numeric_limits<int>::min();
                double max_val = std::numeric_limits<double>::min();
                double val;

                for (int j = 0; j < n_components; j++)
                {
                        val = gsl_matrix_get(W, i, j);
                        if (val > max_val)
                        {
                                max_val = val;
                                max_idx = j;
                        }
                }
                row_clusters.push_back(max_idx);
        }


        for (int j = 0; j < nFeatures; j++)
        {
                int max_idx = std::numeric_limits<int>::min();
                double max_val = std::numeric_limits<double>::min();
                double val;

                for (int i = 0; i < n_components; i++)
                {
                        val = gsl_matrix_get(H, i, j);
                        if (val > max_val)
                        {
                                max_val = val;
                                max_idx = i;
                        }
                }
                column_clusters.push_back(max_idx);
        }
        return 0;
}


vector<int>*
NMF::getRowClusters()
{
        return &row_clusters;
}


vector<int>*
NMF::getColumnClusters()
{
        return &column_clusters;
}

