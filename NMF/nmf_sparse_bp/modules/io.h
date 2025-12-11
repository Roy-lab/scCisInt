#include <gsl/gsl_matrix.h>
#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#ifndef _io_
#define _io_
using namespace std;
namespace io {
        int
        read_sparse_matrix(string, gsl_matrix *);

        //SR added this
        int
        read_dense_matrix(string, gsl_matrix *);

        int
        write_dense_matrix(string, gsl_matrix *);

        int
        write_dense_matrix_trans(string, gsl_matrix *);

        int
        write_list(string, list<double> &);

        int
        read_tree(string, vector<int> &, vector<string> &, vector<string> &, vector<int> &);

        int
        print_usage(string);

        int
        read_prev_results(string &in_dir, gsl_matrix *U, gsl_matrix *V);

        int
        write_mem_and_time(string out_file, unsigned long time_diff, unsigned long mem_diff);

        int
        write_nmf_output(gsl_matrix *U, gsl_matrix *V, string &out_dir);

        template<typename T>
        int write_vector (string outfile, vector<T> * vector)
        {
                ofstream fout(outfile);
                if (!fout.is_open())
                {
                        cerr << "Cannot open file " << outfile << " for writing." << endl;
                        return -1;
                }

                for (int i = 0; i < vector->size(); i++)
                {
                        fout << vector->at(i) << endl;
                }

                fout.close();
                return 0;
        }

};
#endif
