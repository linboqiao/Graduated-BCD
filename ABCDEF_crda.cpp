//============================================================================
// Name        : ABCDEF
// Author      : linboqiao & qianglan
// Version     :
// Copyright   : Your copyright notice
// Description : ABCDEF coordinate randomized asynchronous
//============================================================================

#include <fstream>
#include "omp.h"
#include "ulti.h"

int ABCDEF_coord_rand_distributed_asynchronous(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads, char *funcname)
{
	sprintf(funcname,"%s",__func__);
	printf("enter %s\n",funcname);

	MatrixXd x     = 0.5*MatrixXd::Ones(d, 1);
	MatrixXd x_pre = 0.5*MatrixXd::Ones(d, 1);
	MatrixXd s   = MatrixXd::Zero(m, 1);
	double H     = 0;
	double nu    = 1;
	double tolnu = 1e-7;
	double tol   = 1e-3;
	double gamma = 0.8;
	MatrixXd AT  = A->adjoint();
	//L          = eigs(A'*A, 1);
	//EigenSolver<MatrixXd> eigvals(ATA, false);
	//double L   = eigvals.pseudoEigenvalueMatrix().diagonal().maxCoeff();
	double L     = 1056.56;


	MatrixXd CT  = C->adjoint();
	MatrixXd CCT = (*C)*(C->adjoint());
	MatrixXd terror = MatrixXd::Zero(n, 1);
	MatrixXd tobj   = MatrixXd::Zero(5, 1);

	int k = 0;
	int j = 0;
	int iter_outer = 0;
	int	inneriter = 0.01*d;
	double x_temp = 0;

	printf("start iterate\n");

	// The counter
	int checki = 1;					// the frequency of recording.
	int Num_i = 0;
	double time_start = timing();	// starting point of time.

	// main loop
	srand (time(NULL));
	while (nu>tolnu  &&  k<maxiter)
	{
		for(j = 0;j<inneriter;j++)  // inner loop
		{
			x_pre = x;
			//s = 1./(b-C*x);
			s = ((*b)-(*C)*x).cwiseInverse();

			//H = 2*(L(idx) + nu*s'*CCT*s);
			H = 2*(L + (nu*(s.adjoint())*CCT*s)(0,0));

			omp_set_num_threads(n_threads);
            #pragma omp parallel proc_bind(spread) num_threads(n_threads) //shared(x, AT, A, Y, nu, CT, s, H, lambda, theta,x_temp)
			{
				int icur = omp_get_thread_num();
				int numT = floor(d/n_threads)+1;
				int lb = icur*numT;
				int ub = (icur+1)*numT-1;
				if (ub>d-1)
					ub = d-1;

				//x(idx) = proximalRegC( (A(:,idx)'*(A*x-y) + nu*C(:,idx)'*s)/H, 1, mu/H, theta, regtype);
				#pragma omp for
				for(int idx_d=lb;idx_d<=ub;idx_d++)
				{
					int idx = rand()%d;
					x_temp = x(idx,0)-( ((AT.row(idx))*((*A)*x-(*Y)))(0,0) + nu*((CT.row(idx))*s)(0,0) )/H;
					double xtt = 0;
					proximalSCAD_d(&xtt, &x_temp, 1, lambda/H, theta);
					x(idx) = xtt;
				}
			}

			//if ((x-x_pre).squaredNorm() < 0.01*nu);	break;
		}
		k += 1;
		nu *= gamma;

		// log
		Num_i = Num_i + 1;
		if (Num_i >= checki)
		{
			(*iters)(iter_outer,0)= k;
			(*trace_beta).col(iter_outer) = x.col(0);
			(*trace_time)(iter_outer,0) = timing()-time_start;

			//In Eigen, there is no norm(X)!!!!!
			double err = 0;
			terror = (*Y)-((*A)*x);
			err = normSquared(&terror);
			(*trace_obj)(iter_outer,0) = funSCAD(&x, d, lambda, theta) + 0.5*err;
			Num_i = 0;
		}
		printf("iter_outer, %d, obj, %.7f, time cost, %.7f\n",(int)(*iters)(iter_outer,0), (*trace_obj)(iter_outer,0), (*trace_time)(iter_outer,0));

		iter_outer += 1;

		// Stopping Criterion
		int num_last = 3;
		if (iter_outer>= num_last+1)
		{
			tobj = trace_obj->middleRows(iter_outer-num_last-1, num_last) - trace_obj->middleRows(iter_outer-num_last, num_last);
			double max_err = (tobj.array().abs().maxCoeff());
			if( max_err < tol)
				return 0;
		}
	}
	return 0;
}


