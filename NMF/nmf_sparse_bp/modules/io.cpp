#include <iostream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "io.h"
#include <sstream>


int
io::print_usage(string inputFile)
{
        ifstream f(inputFile.c_str());
        string line;
        while (getline(f, line))
        {
                cout << line << endl;
        }
        f.close();
        return 0;
}

int
io::write_dense_matrix(string outputFile, gsl_matrix *X)
{
        int rowNum = X->size1;
        int colNum = X->size2;
        ofstream ofs;
        ofs.open(outputFile.c_str());
        for (int i = 0; i < rowNum; i++)
        {
                for (int j = 0; j < colNum; j++)
                {
                        ofs << X->data[i * X->tda + j] << "\t";
                }
                ofs << endl;
        }
        ofs.close();
        return 0;
}


int
io::write_dense_matrix_trans(string outputFile, gsl_matrix *X)
{
        int rowNum = X->size1;
        int colNum = X->size2;
        ofstream ofs;
        ofs.open(outputFile.c_str());
        for (int j = 0; j < colNum; j++)
        {
                for (int i = 0; i < rowNum; i++)
                {
                        ofs << X->data[i * X->tda + j] << "\t";
                }
                ofs << endl;
        }
        ofs.close();
        return 0;
}

int
io::read_sparse_matrix(string inputFile, gsl_matrix *X)
{
        int rowNum = X->size1;
        int colNum = X->size2;
        ifstream input(inputFile.c_str());
        int i, j;
        double val;
        while (input >> i >> j >> val)
        {
                //val = log2(val+1);
                gsl_matrix_set(X, i, j, val);
                gsl_matrix_set(X, j, i, val);
        }
        input.close();
        return 0;
}


int
io::read_dense_matrix(string inputFile, gsl_matrix *X)
{
        int rowNum = X->size1;
        int colNum = X->size2;
        ifstream input(inputFile.c_str());
        string buffstr;
        string tok;
        double val;

        if (!input.is_open())
        {
                cerr << "Error opening file: " << inputFile << endl;
                return -1;
        }

        // --- Check file is tab-delimited ---
        // Peek at the first line, then rewind so the main loop sees it too.
        if (getline(input, buffstr))
        {
                if (buffstr.find('\t') == string::npos)
                {
                        cerr << "Error: file does not appear to be tab-delimited: "
                             << inputFile << endl;
                        cerr << "  First line preview: "
                             << buffstr.substr(0, 80) << endl;
                        input.close();
                        return -1;
                }
                input.seekg(0);
        }

        int rowid = 0;
        while (getline(input, buffstr))
        {
                // --- Row bounds check ---
                if (rowid >= rowNum)
                {
                        cerr << "Error: " << inputFile << " has more rows than expected ("
                             << rowNum << ")." << endl;
                        input.close();
                        return -1;
                }

                stringstream ss(buffstr);
                int colid = 0;
                while (getline(ss, tok, '\t'))
                {
                        // Skip empty tokens produced by trailing tabs.
                        if (tok.empty()) { continue; }

                        // --- Column bounds check ---
                        if (colid >= colNum)
                        {
                                cerr << "Error: row " << rowid << " of " << inputFile
                                     << " has more columns than expected (" << colNum << ")." << endl;
                                input.close();
                                return -1;
                        }

                        // --- Safe numeric parse ---
                        try
                        {
                                val = stod(tok);
                        }
                        catch (const invalid_argument &)
                        {
                                cerr << "Error: non-numeric value \"" << tok
                                     << "\" at row " << rowid << ", col " << colid
                                     << " in " << inputFile << endl;
                                input.close();
                                return -1;
                        }
                        catch (const out_of_range &)
                        {
                                cerr << "Error: value out of double range \"" << tok
                                     << "\" at row " << rowid << ", col " << colid
                                     << " in " << inputFile << endl;
                                input.close();
                                return -1;
                        }

                        gsl_matrix_set(X, rowid, colid, val);
                        colid++;
                }

                // --- Column count check ---
                if (colid != colNum)
                {
                        cerr << "Error: row " << rowid << " of " << inputFile
                             << " has " << colid << " column(s) but expected "
                             << colNum << "." << endl;
                        input.close();
                        return -1;
                }

                rowid++;
        }

        // --- Row count check ---
        if (rowid != rowNum)
        {
                cerr << "Error: " << inputFile << " has " << rowid
                     << " row(s) but expected " << rowNum << "." << endl;
                input.close();
                return -1;
        }

        input.close();
        return 0;
}


int
io::write_list(string outputFile, list<double> &err)
{
        ofstream ofs;
        ofs.open(outputFile.c_str());
        for (list<double>::iterator itr = err.begin(); itr != err.end(); ++itr)
        {
                ofs << *itr << endl;
        }
        ofs.close();
        return 0;
}

int
io::read_tree(string inputFile,
              vector<int> &parentIds, vector<string> &aliases, vector<string> &fileNames, vector<int> &numSamples)
{
        ifstream input(inputFile.c_str());
        int id, pid;
        string alias, filename, n;
        while (input >> id >> pid >> alias >> filename >> n)
        {
                parentIds.push_back(pid);
                aliases.push_back(alias);
                if (filename != "N/A")
                {
                        fileNames.push_back(filename);
                        stringstream nTemp(n);
                        int numSample = 0;
                        nTemp >> numSample;
                        numSamples.push_back(numSample);
                }
        }
        input.close();
        return 0;
}

int
io::read_prev_results(string &in_dir, gsl_matrix *U, gsl_matrix *V)
{
        io::read_dense_matrix(in_dir + "U.txt", U);
        io::read_dense_matrix(in_dir + "V.txt", V);
        return 0;
}

int
io::write_mem_and_time(string out_file, unsigned long int time_diff, unsigned long int mem_diff)
{
        ofstream ofs;
        ofs.open(out_file.c_str());
        ofs << "time:\t" << time_diff << endl;
        ofs << "mem:\t" << mem_diff << endl;
        ofs.close();
        return 0;
}

int
io::write_nmf_output(gsl_matrix *U, gsl_matrix *V, string &out_dir)
{
        io::write_dense_matrix(out_dir + "U.txt", U);
        io::write_dense_matrix_trans(out_dir + "V.txt", V);
        //io::write_dense_matrix(out_dir + "R.txt", R);
        return 0;
}





