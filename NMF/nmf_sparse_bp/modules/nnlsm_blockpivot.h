//
// Created by Spencer Halberg on 6/8/23.
//

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include "io.h"

#ifndef TMF_NNLSM_BLOCKPIVOT_H
#define TMF_NNLSM_BLOCKPIVOT_H

using namespace std;


typedef pair<int, int> matrixIndex;
namespace nnlsm_bp {
        int
        nnlsmBlockPivot(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);

        int
        solveNormalEqComb(gsl_matrix *, gsl_matrix *, set<int> &nonOptSet, vector<vector<bool>> &passSet, gsl_matrix *Z,
                          int &iter);


        int
        checkPositivity(gsl_matrix *X, vector<vector<bool>> &Set, bool negate);

        int
        negateBoolMatrix(vector<vector<bool>> &M);

        int
        intersectBoolMatrix(vector<vector<bool>> &, vector<vector<bool>> &, vector<vector<bool>> &);

        int
        sumMatrixIntColumn(vector<vector<bool>> &X, vector<int> &sum);

        int
        makePassiveSetMap(vector<vector<bool>> &passSet, multimap<vector<bool>, int> &passSetMap);


        int
        findNonOptCols(vector<int> &NotGood, set<int> &nonOptCols);

        int
        updatePassiveSet(set<int> &NonOptCols, vector<int> &NotGood, vector<int> &Tcol, vector<int> &Pcol,
                         vector<vector<bool>> &passSet, vector<vector<bool>> &NonOptSet,
                         vector<vector<bool>> &InfeaSet);


        int
        inRowColumnMaskedMatrix(set<int> &row_mask, set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask);

        int
        outRowColumnMaskedMatrix(set<int> &row_mask, set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask);


        template<class T>
        int
        initializeVectorVector(vector<vector<T>> &Vect, int size1, int size2, T val)
        {
                vector<T> v;
                for (int j = 0; j < size2; j++)
                {
                        v.push_back(val);
                }
                Vect.erase(Vect.begin(), Vect.end()); //Ensure the vector is empty and then fill with value val.
                for (int i = 0; i < size1; i++)
                {
                        Vect.push_back(v);
                }
                return 0;
        }

        template<class T>
        int
        initializeVector(vector<T> &Vect, int size, T val)
        {
                Vect.erase(Vect.begin(), Vect.end());
                for (int i = 0; i < size; i++)
                {
                        Vect.push_back(val);
                }
                return 0;
        }


        int
        inColumnMaskedMatrix(set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask);

        int
        outColumnMaskedMatrix(set<int> &col_mask, gsl_matrix *X, gsl_matrix *X_mask);

        int
        outRowMaskedMatrix(set<int> &row_mask, gsl_matrix *X, gsl_matrix *X_mask);
}
#endif //TMF_NNLSM_BLOCKPIVOT_H



