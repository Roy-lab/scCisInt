//
// Created by Spencer Halberg on 6/6/23.
//

#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <list>
#include <string>
#include <sstream>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include "modules/io.h"
#include "modules/nmf.h"
#include "modules/utils.h"
#include "modules/initialization.h"
#include <string.h>


int mkdir_recursive(const std::string &dir) {
	std::stringstream path;
	std::string token;
	std::istringstream ss(dir);

	while (std::getline(ss, token, '/')) {
		path << token << "/";
		std::string p = path.str();
		if (mkdir(p.c_str(), 0766) != 0 && errno != EEXIST) {
			return false;
		}
	}
	return 0;
}


int
main(int argc, char **argv)
{
	struct timeval beginTime;
	struct timeval factorTime;
	struct timeval endTime;
	gettimeofday(&beginTime, NULL);

	struct rusage bUsage;
	struct rusage eUsage;
	getrusage(RUSAGE_SELF, &bUsage);

	string matrixFile;
	int nSamples = -1;
	int nFeatures = -1;
	int nComponents = -1;

	string outputPrefix = string(".");
	string inPrefix = string("");
	int seed = 1010;
	int verbose = false;
	int maxIter = 100;
	int minIter = 20;
	double tol = 1e-3;
	double alpha = 0.1;
	double beta = 0.1;

	method_type method = block_pivot;
	reg_type reg = sparse;
	stop_rule_type stop = normalized_pgrad;

	//Get option_longs
	static struct option long_options[] = {
		{"data", required_argument, 0, 'x'},
		{"nRows", required_argument, 0, 'n'},
		{"nCols", required_argument, 0, 'N'},
		{"nFactors", required_argument, 0, 'k'},
		{"output", required_argument, 0, 'o'},
		{"inputPrefix", required_argument, 0, 'p'},
		{"seed", required_argument, 0, 'r'},
		{"verbose", no_argument, 0, 's'},
		{"maxIter", required_argument, 0, 'M'},
		{"minIter", required_argument, 0, 'm'},
		{"tolerance", required_argument, 0, 't'},
		{"alpha", required_argument, 0, 'a'},
		{"beta", required_argument, 0, 'b'},
		{"reg", required_argument, 0, 'R'},
		{"stop", required_argument, 0, 'S'},
		{"help", no_argument, 0, 'h'},
		{0,0,0,0}
	};

	string usage = string("usage.txt");
	string reg_in = "";
	string stop_in = "";
	int c;
	int option_index = 0;
	while ((c = getopt_long(argc, argv, "x:n:N:k:o:p:r:sM:m:t:a:b:R:S:h",
				long_options, &option_index)) != -1)
	{
		switch(c) {
			case 'x': matrixFile = optarg; break;
			case 'n': nSamples = atoi(optarg); break;
			case 'N': nFeatures = atoi(optarg); break;
			case 'k': nComponents = atoi(optarg); break;
			case 'o': outputPrefix = optarg; break;
			case 'p': inPrefix = optarg; break;
			case 'r': seed = atoi(optarg); break;
			case 's': verbose = true; break;
			case 'M': maxIter = atoi(optarg); break;
			case 'm': minIter = atoi(optarg); break;
			case 't': tol = atof(optarg); break;
			case 'a': alpha = atof(optarg); break;
			case 'b': beta = atof(optarg); break;
			case 'R':
				if (string(optarg) == "normal") reg = normal;
				else if (string(optarg) == "sparse") reg = sparse;
				else if (string(optarg) == "regularized") reg = regularized;
				else { cerr << "Invalid reg type" << endl; return -1; }
				break;
			case 'S':
				if (string(optarg) == "normalized_pgrad") stop = normalized_pgrad;
				else if (string(optarg) == "pgrad") stop = pgrad;
				else { cerr << "Invalid stop rule" << endl; return -1; }
				break;
			case 'h':
				io::print_usage(usage);
				return 0;
			default:
				io::print_usage(usage);
				return 1;
		}
	}

	if (matrixFile.empty() || nSamples <= 0 || nFeatures <= 0 || nComponents <= 0) {
		cerr << "Error: --data, --nRows, --nCols, --nFactors are required\n";
		io::print_usage(usage);
		return 1;
	}


	const gsl_rng_type *T;
	gsl_rng *ri;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, seed);



	//Read Initial Matrix
	gsl_matrix *X = gsl_matrix_calloc(nSamples, nFeatures);
	int read_fail = io::read_dense_matrix(matrixFile, X);
	if (read_fail != 0)
	{
		cout << "Failed to read input data file." << endl;
		return -1;
	}
	if (utils::is_all_zeros(X, "Data"))
	{
		return -1;
	}
	if (utils::has_nan(X, "Data"))
	{
		return -1;
	}



        gsl_matrix *U = gsl_matrix_calloc(nSamples, nComponents);
        gsl_matrix *V = gsl_matrix_calloc(nComponents, nFeatures);


        if (inPrefix == "")
        {
                init::random(U, ri);
                init::random(V, ri);
        } else
        {
                io::read_prev_results(inPrefix, U, V);
        	if (utils::is_all_zeros(U, "Previous U"))
        	{
        		return -1;
        	}
        	if (utils::is_all_zeros(V, "Previous V"))
        	{
        		return -1;
        	}
        }

        stringstream out_dir_str;
        out_dir_str << outputPrefix << "/"; //<< "/k_" << nComponents << "/";
        string out_dir = out_dir_str.str();

	if (verbose)
	{
		cout << "Making " << out_dir_str.str() << endl;
	}
	mkdir_recursive(out_dir);

        //io::write_dense_matrix(out_dir + "U_init.txt", U);
        //io::write_dense_matrix(out_dir + "V_init.txt", V);


        NMF nmf = NMF(nComponents, reg, method, stop, maxIter, minIter, verbose, tol, alpha, beta);
	nmf.setSize(nSamples, nFeatures);


        nmf.fit(X, U, V);
	gettimeofday(&factorTime, NULL);

	nmf.assignClusters();

        gettimeofday(&endTime, NULL);

        getrusage(RUSAGE_SELF, &eUsage);

        unsigned long int bt = beginTime.tv_sec;
        unsigned long int ft = factorTime.tv_sec;
        unsigned long int et = endTime.tv_sec;

	unsigned long int totalSeconds  = et - bt;
	unsigned long int hours         = totalSeconds / 3600;
	unsigned long int minutes       = (totalSeconds % 3600) / 60;
	unsigned long int seconds       = totalSeconds % 60;

	cout << "Total time elapsed: ";
	if (hours > 0)
		cout << hours << "h ";
	if (hours > 0 || minutes > 0)
		cout << minutes << "m ";
	cout << seconds << "s" << endl;

        cout << "Total time elapsed: " << et - bt << " seconds" << endl;

        unsigned long int bu = bUsage.ru_maxrss;
        unsigned long int eu = eUsage.ru_maxrss;

        cout << "Memory usage: " << (eu - bu) / 1000 << "MB" << endl;

        io::write_mem_and_time(out_dir + "usage.txt", et - bt, (eu - bu) / 1000);
        io::write_nmf_output(U, V, out_dir);

	vector<int>* U_clusters = nmf.getRowClusters();
	vector<int>* V_clusters = nmf.getColumnClusters();

	io::write_vector(out_dir + "/U_assign.txt", U_clusters);
	io::write_vector(out_dir + "/V_assign.txt", V_clusters);


        return 0;
}
