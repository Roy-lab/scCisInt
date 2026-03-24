#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <string>

#ifndef _utils_
#define _utils_
using namespace std;

namespace utils {
        double
        get_frobenius_norm(gsl_matrix *);

        int
        get_inverse(gsl_matrix *);

        // Returns true if every element of X is exactly 0.0.
        // Optionally prints a labelled warning to stderr.
        bool
        is_all_zeros(gsl_matrix *X, string label = "");

        // Returns true if X contains any NaN values.
        // Optionally prints count and first location to stderr.
        bool
        has_nan(gsl_matrix *X, string label = "");
};
#endif
