#include <gsl/gsl_matrix.h>
#include <list>
#include "nnlsm_blockpivot.h"

#ifndef _nmf_
#define _nmf_
using namespace std;


enum reg_type
{
        normal, regularized, sparse
};
enum method_type
{
        block_pivot, active_set
};
enum stop_rule_type
{
        normalized_pgrad, pgrad
};

class NMF
{
public:
        NMF(int, reg_type, method_type, stop_rule_type, int, int, bool, double, double, double);
        int setSize(int, int);
        int fit(gsl_matrix *, gsl_matrix *, gsl_matrix *);
        int assignClusters();

        vector<int>* getRowClusters();
        vector<int>* getColumnClusters();

        ~NMF();


private:
        gsl_matrix *A;
        gsl_matrix *W;
        gsl_matrix *H;
        gsl_matrix *gradW;
        gsl_matrix *gradH;

        double tol;
        double alpha;
        double beta;
        int nSamples;
        int nFeatures;
        int max_iter;
        int min_iter;
        int random_state;
        bool verbose;
        int n_components;

        vector<int> row_clusters;
        vector<int> column_clusters;
        reg_type type;
        method_type method;
        stop_rule_type stop_rule;

        int nnlsm(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);

        int getGradient();

        int appendRows(gsl_matrix *, gsl_matrix *, gsl_matrix *);

        int appendRows(gsl_matrix *append, gsl_matrix *Block1, gsl_vector *Block2);

        int gslMatrixRowsFlatVect(gsl_vector *, gsl_matrix *);

        int gslMatrixColFlatVect(gsl_vector *Vector, gsl_matrix *Matrix);

        double getStopValue();

        double getStartValue();
	
	int rescaleVectors();
};

#endif
