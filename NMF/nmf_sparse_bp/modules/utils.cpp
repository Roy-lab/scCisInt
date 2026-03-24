#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <list>
#include <vector>
#include <math.h>
#include <iostream>
#include <string>
#include "utils.h"

double
utils::get_frobenius_norm(gsl_matrix *X)
{
        int n = X->size1;
        int m = X->size2;
        double sum = 0;
        gsl_vector *row_copy = gsl_vector_alloc(m);
        for (int i = 0; i < n; i++)
        {
                gsl_vector_view row = gsl_matrix_row(X, i);
                gsl_vector_memcpy(row_copy, &row.vector);
                gsl_vector_mul(row_copy, &row.vector);
                double rowsum = gsl_blas_dasum(row_copy);
                sum += rowsum;
        }
        gsl_vector_free(row_copy);
        return sum;
}

int
utils::get_inverse(gsl_matrix *A)
{
        gsl_linalg_cholesky_decomp1(A);
        gsl_linalg_cholesky_invert(A);
        return 0;
}

// Returns true if every element of X is exactly 0.0.
// Prints a warning to stderr if so, using label to identify the matrix.
bool
utils::is_all_zeros(gsl_matrix *X, string label)
{
        int n = X->size1;
        int m = X->size2;
        for (int i = 0; i < n; i++)
        {
                for (int j = 0; j < m; j++)
                {
                        if (gsl_matrix_get(X, i, j) != 0.0)
                                return false;
                }
        }
        if (!label.empty())
        {
                cerr << "Warning: " << label << " is an all-zero matrix ("
                     << n << " x " << m << ")." << endl;
        }
        return true;
}

// Returns true if X contains one or more NaN values.
// Prints the total NaN count and the location of the first one to stderr.
bool
utils::has_nan(gsl_matrix *X, string label)
{
        int n = X->size1;
        int m = X->size2;
        int count    = 0;
        int first_row = -1;
        int first_col = -1;
        for (int i = 0; i < n; i++)
        {
                for (int j = 0; j < m; j++)
                {
                        if (isnan(gsl_matrix_get(X, i, j)))
                        {
                                if (count == 0)
                                {
                                        first_row = i;
                                        first_col = j;
                                }
                                count++;
                        }
                }
        }
        if (count > 0)
        {
                string prefix = label.empty() ? "Matrix" : label;
                cerr << "Warning: " << prefix << " contains " << count
                     << " NaN value(s). First at ["
                     << first_row << ", " << first_col << "]." << endl;
        }
        return count > 0;
}
