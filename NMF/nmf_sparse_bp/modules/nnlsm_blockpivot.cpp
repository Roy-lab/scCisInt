//
// Created by Spencer Halberg on 6/6/23.
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include "nnlsm_blockpivot.h"
#include <iostream>
#include "io.h"


using namespace std;


//This is based on the the function nnlsm_blockpivot from
//Nonnegativity Constrained Least Squares with Multiple Righthand Sides
//using Block Principal Pivoting method by Jingu Kim and Haesun Park


//If using cite: Jingu Kim and Haesun Park, Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons,
//In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM'08), 353-362, 2008

//A: input matrix (n * k) for W or (m * k) for H. 
//B: input matrix (n * m) data matrix X or (m*n) data transpose X' for updates of W or H respectively. 
//X: output solution
//Y: Gradient



int
nnlsm_bp::nnlsmBlockPivot(gsl_matrix *A, gsl_matrix *B, gsl_matrix *X, gsl_matrix *Y)
{
        //Initialize Key Variables.
        int m = A->size1;
        int k = A->size2;
        int n = B->size2;
        int iter = 0;

        gsl_matrix *AtA = gsl_matrix_calloc(k, k);
        gsl_matrix *AtB = gsl_matrix_calloc(k, n);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, AtA);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, B, 0, AtB);

        int MAX_ITER = k * 5; // Default in code

        //io::write_dense_matrix("tmp/X.txt", X);
        //io::write_dense_matrix("tmp/Y.txt", Y);
        //io::write_dense_matrix("tmp/A.txt", A);
        //io::write_dense_matrix("tmp/B.txt", B);
        //io::write_dense_matrix("tmp/AtA.txt", AtA);
        //io::write_dense_matrix("tmp/AtB.txt", AtB);

        multimap<vector<bool>, int> PassSetMap;
        vector<vector<bool>> PassSet;
        initializeVectorVector<bool>(PassSet, n, k, false);
        vector<vector<bool>> negX;
        initializeVectorVector<bool>(negX, n, k, false);
        vector<vector<bool>> negY;
        initializeVectorVector<bool>(negY, n, k, false);
        vector<vector<bool>> NonOptSet;
        initializeVectorVector<bool>(NonOptSet, n, k, false);
        vector<vector<bool>> InfeaSet;
        initializeVectorVector<bool>(InfeaSet, n, k, false);

        vector<int> NotGood;
        initializeVector<int>(NotGood, n, 0);
        set<int> NonOptCol;

        vector<int> P_cols;
        initializeVector(P_cols, n, 3);
        vector<int> T_cols;
        initializeVector(T_cols, n, k + 1);


        //Initialize the original Matrix Note that PassSet is empty.

        //if (gsl_matrix_ispos(X) == 0)
        //{
        //        gsl_matrix_memcpy(Y, AtB);
        //        gsl_matrix_scale(Y, -1.0);
        //} else
        //{
        // We are always passing an init X so lets remove the if statement
	checkPositivity(X, PassSet, false);
        for (int i = 0; i < PassSet.size(); i++)
        {
        	NonOptCol.insert(i);
        }
        solveNormalEqComb(AtA, AtB, NonOptCol, PassSet, X, iter);
        gsl_matrix_memcpy(Y, AtB);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, AtA, X, -1, Y);
        //}
        NonOptCol.erase(NonOptCol.begin(), NonOptCol.end());

        //io::write_dense_matrix("tmp/X.txt", X);
        //io::write_dense_matrix("tmp/Y.txt", Y);

        // Redefine Sets
        checkPositivity(Y, negY, true);
        checkPositivity(X, negX, true);
        intersectBoolMatrix(negX, PassSet, InfeaSet);
        negateBoolMatrix(PassSet); //Convert Passive set to active set
        intersectBoolMatrix(negY, PassSet, NonOptSet); //Find non-opt set of active set
        negateBoolMatrix(PassSet); // Convert activeSet back to Passive Set

        sumMatrixIntColumn(NonOptSet, NotGood); //Count the number of passive set elements
        sumMatrixIntColumn(InfeaSet, NotGood);

        findNonOptCols(NotGood, NonOptCol); // find the non_opt set;
	int bigIter = 0;
        while (!NonOptCol.empty())
        {
		bigIter = bigIter + 1;
                if (bigIter > MAX_ITER)
                {
                        cout << "WARNING: Reached max iterations without converging." << endl;
                        return 0;
                }
                // Break columns into three sets based on the paper description.
                updatePassiveSet(NonOptCol, NotGood, T_cols, P_cols, PassSet, NonOptSet, InfeaSet);

                // Update X
                solveNormalEqComb(AtA, AtB, NonOptCol, PassSet, X, iter);


                for (int i = 0; i < X->size1; i++)
                {
                        for (int j = 0; j < X->size2; j++)
                        {
                                if (abs(gsl_matrix_get(X, i, j)) < 1e-12)
                                { // For numerical stability
                                        gsl_matrix_set(X, i, j, 0);
                                }
                        }
                }

                gsl_matrix_memcpy(Y, AtB);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, AtA, X, -1, Y);

                for (int i = 0; i < X->size1; i++)
                {
                        for (int j = 0; j < X->size2; j++)
                        {
                                if (abs(gsl_matrix_get(X, i, j)) < 1e-12)
                                { // For numerical stability
                                        gsl_matrix_set(X, i, j, 0);
                                }
                        }
                }
		//io::write_dense_matrix("tmp/X.txt", X);
        	//io::write_dense_matrix("tmp/Y.txt", Y);

                checkPositivity(Y, negY, true);
                checkPositivity(X, negX, true);
                intersectBoolMatrix(negX, PassSet, InfeaSet);
                negateBoolMatrix(PassSet); //Convert Passive set to active set
                intersectBoolMatrix(negY, PassSet, NonOptSet); //Find non-opt set of active set
                negateBoolMatrix(PassSet); // Convert activeSet back to Passive Set
                initializeVector<int>(NotGood, n, 0);
		sumMatrixIntColumn(NonOptSet, NotGood);
                sumMatrixIntColumn(InfeaSet, NotGood);
                findNonOptCols(NotGood, NonOptCol); // find the non_opt set;

        }

        return 0;
}

int
nnlsm_bp::checkPositivity(gsl_matrix *X, vector<vector<bool>> &Set, bool negate)
{
        int m = X->size1;
        int n = X->size2;
        if (Set.size() != n || Set[1].size() != m)
        {
                cout << "Passive set of inconsistent size." << endl;
                return 1;
        }

        double val = gsl_matrix_get(X, 0, 0);  // Dummy value.
        for (int i = 0; i < m; i++)
        {
                for (int j = 0; j < n; j++)
                {
                        val = gsl_matrix_get(X, i, j);
                        if (!negate)
                        {
                                Set[j][i] = val > 0;
                        } else
                        {
                                Set[j][i] = val < 0;
                        }

                }
        }
        return 0;
}

int
nnlsm_bp::negateBoolMatrix(vector<vector<bool>> &M)
{
        int n = M.size();
        int k = M[0].size();

        for (int i = 0; i < n; i++)
        {
                for (int j = 0; j < k; j++)
                {
                        M[i][j] = !M[i][j];
                }
        }

        return 0;
}

int
nnlsm_bp::intersectBoolMatrix(vector<vector<bool>> &A, vector<vector<bool>> &B, vector<vector<bool>> &Intersect)
{
        //check sizes
        if (A.size() != B.size())
        {
                cout << "Matrix A and B must be the same size to intersect." << endl;
                return -999;
        }

        if (A.size() != Intersect.size())
        {
                cout << "Matrix A and result matrix must be the same size to intersect." << endl;
                return -999;
        }

        for (int i = 0; i < A.size(); i++)
        {
                int k = A[i].size();
                for (int j = 0; j < k; j++)
                {
                        if (A[i].size() != B[i].size())
                        {
                                cout << "Matrix A and B must be the same size to intersect." << endl;
                                return -999;
                        }
                        if (A[i].size() != Intersect[i].size())
                        {
                                cout << "Matrix A and result matrix  must be the same size to intersect." << endl;
                                return -999;
                        }
                        Intersect[i][j] = A[i][j] & B[i][j];

                }
        }
        return 0;
}


int
nnlsm_bp::sumMatrixIntColumn(vector<vector<bool>> &X, vector<int> &sum)
{
        int n = X.size();
        for (int i = 0; i < n; i++)
        {
                sum[i] = sum[i] + count_if(X[i].begin(), X[i].end(), [](bool i) { return i; });
        }
        return 0;
}

int
nnlsm_bp::findNonOptCols(vector<int> &NotGood, set<int> &nonOptCols)
{
        int n = NotGood.size();
        for (int i = 0; i < n; i++)
        {
                if (NotGood[i] > 0)
                {
                        nonOptCols.insert(i);
                } else
                {
                        nonOptCols.erase(i);
                }
        }
        return 0;
}

int
nnlsm_bp::updatePassiveSet(set<int> &NonOptCols, vector<int> &NotGood, vector<int> &Tcol, vector<int> &Pcol,
                           vector<vector<bool>> &passSet, vector<vector<bool>> &NonOptSet,
                           vector<vector<bool>> &InfeaSet)
{
        //Three partition of NonOptCols to derive passive set mask. Col1 (first Update), Col2 (second or third update), Col3 (ill-conditioned so pivot).
        set<int> Col1;
        set<int> Col2;
        set<int> Col3;
        int k = passSet[0].size();

        for (int i:NonOptCols)
        {
                if (NotGood[i] < Tcol[i])
                {
                        Col1.insert(i);
                        Col2.erase(i);
                        Col3.erase(i);
                        Pcol[i] = 3;
                        Tcol[i] = NotGood[i];

                } else if (NotGood[i] >= Tcol[i] && Pcol[i] >= 1)
                {
                        Col1.erase(i);
                        Col2.insert(i);
                        Col3.erase(i);
                        Pcol[i] = Pcol[i] - 1;
                } else
                {
                        Col1.erase(i);
                        Col2.erase(i);
                        Col3.insert(i);
                }
        }

        if (!Col3.empty())
        {
                for (int j:Col3)
                {
                        int max_idx = 0;
                        for (int i = 0; i < k; i++)
                        {
                                if (NonOptSet[j][i] || InfeaSet[j][i])
                                {
                                        max_idx = i;
                                }
                        }
                        passSet[j][max_idx] = !passSet[j][max_idx]; //Flip the one element. For extreme cases!
                }
        } else
        {
                set<int> U;
                set_union(Col1.begin(), Col1.end(), Col2.begin(), Col2.end(), std::inserter(U, U.end()));
                for (int j:U)
                {
                        for (int i = 0; i < k; i++)
                        {
                                if (NonOptSet[j][i] == 1)
                                {
                                        passSet[j][i] = true;
                                } else if (InfeaSet[j][i] == 1)
                                {
                                        passSet[j][i] = false;
                                }
                        }
                }
        }
        return 0;
}


int
nnlsm_bp::inRowColumnMaskedMatrix(set<int> &row_mask, set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask)
{
        int mask_j = 0;
        for (int j:col_mask)
        {
                gsl_vector_view X_mask_j = gsl_matrix_column(X_mask, mask_j);
                gsl_vector_view X_j = gsl_matrix_column(X, j);

                int mask_i = 0;
                for (int i:row_mask)
                {
                        gsl_vector_set(&X_mask_j.vector, mask_i, gsl_vector_get(&X_j.vector, i));
                        mask_i++;
                }
                mask_j++;
        }
        return 0;
}

int
nnlsm_bp::inColumnMaskedMatrix(set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask)
{
        int mask_j = 0;
        for (int j:col_mask)
        {
                gsl_vector_view X_mask_j = gsl_matrix_column(X_mask, mask_j);
                gsl_vector_view X_j = gsl_matrix_column(X, j);
                gsl_vector_memcpy(&X_mask_j.vector, &X_j.vector);
                mask_j++;
        }
}


int
nnlsm_bp::outRowColumnMaskedMatrix(set<int> &row_mask, set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask)
{
        int mask_j = 0;
        for (int j:col_mask)
        {
                gsl_vector_view X_mask_j = gsl_matrix_column(X_mask, mask_j);
                gsl_vector_view X_j = gsl_matrix_column(X, j);
                int mask_i = 0;
                for (int i:row_mask)
                {
                        gsl_vector_set(&X_j.vector, i, gsl_vector_get(&X_mask_j.vector, mask_i));
                        mask_i++;
                }
                mask_j++;
        }
        return 0;
}

int
nnlsm_bp::outColumnMaskedMatrix(set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask)
{
        int mask_j = 0;
        for (int j:col_mask)
        {
                gsl_vector_view X_mask_j = gsl_matrix_column(X_mask, mask_j);
                gsl_vector_view X_j = gsl_matrix_column(X, j);
                gsl_vector_memcpy(&X_j.vector, &X_mask_j.vector);
                mask_j++;
        }
}

int
nnlsm_bp::outRowMaskedMatrix(set<int> &row_mask, gsl_matrix *X, gsl_matrix *X_mask)
{
        int mask_i = 0;
        for (int i:row_mask)
        {
                gsl_vector_view X_mask_i = gsl_matrix_row(X_mask, mask_i);
                gsl_vector_view X_i = gsl_matrix_row(X, i);
                gsl_vector_memcpy(&X_i.vector, &X_mask_i.vector);
                mask_i++;
        }
}


int
nnlsm_bp::makePassiveSetMap(vector<vector<bool>> &passSet, multimap<vector<bool>, int> &passSetMap)
{
        passSetMap.erase(passSetMap.begin(), passSetMap.end()); //Clear map
        int n = passSet.size();
        for (int i = 0; i < n; i++)
        {
                passSetMap.insert({passSet[i], i});
        }
        return 0;
}

int
nnlsm_bp::solveNormalEqComb(gsl_matrix *AtA, gsl_matrix *AtB, set<int> &nonOptSet, vector<vector<bool>> &passSet,
                            gsl_matrix *Z, int &iter)
{

        multimap<vector<bool>, int> passSetMap;
        makePassiveSetMap(passSet, passSetMap);
        set<vector<bool>> passSetKey;
        for (pair<vector<bool>, int> element:passSetMap)
        {
                passSetKey.insert(element.first);
        }


        //For each type of vector in passive set, we need to subset AtA, and AtB.
        set<int> vars;
        set<int> cols;

        for (vector<bool> key:passSetKey)
        {
                vars.erase(vars.begin(), vars.end());
                //select the variables in passive set vector that are 1.
                for (int i = 0; i < key.size(); i++)
                {
                        if (key[i])
                        {
                                vars.insert(i);
                        }
                }


                //select columns in the passive set whose value matches key
                cols.erase(cols.begin(), cols.end());
                for (auto iter = passSetMap.lower_bound(key); iter != passSetMap.upper_bound(key); iter++)
                {
                        if (nonOptSet.find(iter->second) != nonOptSet.end())
                        { //Double check the col needs updated.
                                cols.insert(iter->second);
                        }
                }

                if (!cols.empty() && !vars.empty())
                {

                        gsl_matrix *AtA_mask = gsl_matrix_alloc(vars.size(), vars.size());
                        inRowColumnMaskedMatrix(vars, vars, AtA, AtA_mask);
                        gsl_matrix *AtB_mask = gsl_matrix_alloc(vars.size(), cols.size());
                        inRowColumnMaskedMatrix(vars, cols, AtB, AtB_mask);
                        gsl_matrix *Z_mask_col = gsl_matrix_calloc(Z->size1, cols.size());
                        gsl_matrix *Z_mask_col_row = gsl_matrix_calloc(vars.size(), cols.size());


                        gsl_permutation *p = gsl_permutation_alloc(AtA_mask->size1);
                        int s;
                        gsl_linalg_LU_decomp(AtA_mask, p, &s);
                        for (int j = 0; j < AtB_mask->size2; j++)
                        {
                                gsl_vector_view Z_mask_j = gsl_matrix_column(Z_mask_col_row, j);
                                gsl_vector_view AtB_mask_j = gsl_matrix_column(AtB_mask, j);
                                gsl_linalg_LU_solve(AtA_mask, p, &AtB_mask_j.vector, &Z_mask_j.vector);
                        }
                        outRowMaskedMatrix(vars, Z_mask_col, Z_mask_col_row);
                        outColumnMaskedMatrix(cols, Z, Z_mask_col);

                        iter++;
                        gsl_matrix_free(AtA_mask);
                        gsl_matrix_free(AtB_mask);
                        gsl_matrix_free(Z_mask_col);
                        gsl_matrix_free(Z_mask_col_row);
                }
        }
}






