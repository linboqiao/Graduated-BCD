//============================================================================
// Name        : ABCDEF_asynchronous
// Author      : linboqiao
// Version     :
// Copyright   : Your copyright notice
// Description : main function
//============================================================================

#include <fstream>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <sys/time.h>
#include "ulti.h"

int main() {
	printf("start!\n");
	int d = 1000, n = 20, m = 5;

	//input data
	MatrixXd X = MatrixXd::Zero(n,d);		// array of size nxd
	MatrixXd C = MatrixXd::Zero(m,d);		// array of size mxd
	MatrixXd Y = MatrixXd::Zero(n,1);		// array of size n
	MatrixXd b = MatrixXd::Zero(m,1);		// array of size m
	MatrixXd betaTrue = MatrixXd::Zero(d,1);	// array of size d
	double L   = 0;
	random_init(d, n, m, &X, &C, &Y, &b, &betaTrue, &L);
	//read_from_csv(d, n, m, &X, &C, &Y, &b, &betaTrue, &L);

	int maxiter   = 1000;
	int regtype   = 3;
	double theta  = 1.4;
	double lambda = 2.0;
	char funcname[255] = {0};

	int num_iters = 100;

	//output results
	MatrixXd iters      = MatrixXd::Zero(num_iters,1);
	MatrixXd trace_beta = MatrixXd::Zero(num_iters,d);
	MatrixXd trace_obj  = MatrixXd::Zero(num_iters,1);
	MatrixXd trace_time = MatrixXd::Zero(num_iters,1);

	for (int idx_method =0; idx_method<8;idx_method++)
	{
		int n_threads[] = {1,2,4,8};
		int n_thread = 4;
		for(int i = n_thread-1; i>=0;i--)
		{
			int num_threads = n_threads[i];
			printf("num_threads: %d\n",num_threads);
			double time_start = timing();
			switch(idx_method)
			{
				case 0:
				ABCDEF(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 1:
				ABCDEF_coord_cyclic(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 2:
				ABCDEF_coord_cyclic_asynchronous(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 3:
				ABCDEF_coord_cyclic_distributed(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 4:
				ABCDEF_coord_cyclic_distributed_asynchronous(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 5:
				ABCDEF_coord_rand_distributed(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 6:
				ABCDEF_coord_rand_distributed_asynchronous(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
				case 7:
				ABCDEF_coord_rand_asynchronous(&X, &C, &Y, &b, lambda, regtype, theta, maxiter,\
				&iters, &trace_beta, &trace_obj, &trace_time, m, n, d, num_threads, funcname);
				break;
			}
			double time_end = timing();

			char resultFile[255]= {0};
			sprintf(resultFile, "result_%s_%d_mu%d_d%d_n%d_m%d.csv",\
					funcname, num_threads, (int)(lambda*10), d, n, m);
			//printf("result file name %s\n", resultFile);
			//write_to_csv(resultFile, &iters, &trace_beta, &trace_obj, &trace_time);

			double elapsetime = time_end - time_start;
			printf("time start  is \t%g\n", time_start);
			printf("time end    is \t%g\n", time_end);
			printf("time elapse is \t%g\n\n", elapsetime);
		}
	}
	return 0;
}

