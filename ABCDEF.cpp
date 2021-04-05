//============================================================================
// Name        : ABCDEF
// Author      : linboqiao & qianglan
// Version     :
// Copyright   : Your copyright notice
// Description : ABCDEF batch
//============================================================================

#include <fstream>
#include "omp.h"
#include "ulti.h"

int ABCDEF(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads, char *funcname)
{
	sprintf(funcname,"%s",__func__);
	printf("enter %s\n",funcname);

	MatrixXd x   = 0.5*MatrixXd::Ones(d, 1);
	MatrixXd x_pre = 0.5*MatrixXd::Ones(d, 1);
	MatrixXd x_temp = 0.5*MatrixXd::Ones(d, 1);
	MatrixXd s   = MatrixXd::Zero(m, 1);
	double H     = 0;
	double nu    = 1;
	double tolnu = 1e-7;
	double tol   = 1e-3;
	double gamma = 0.8;
	MatrixXd AT  = A->adjoint();
	MatrixXd ATA = AT*(*A);
	//L          = eigs(A'*A, 1);
	//EigenSolver<MatrixXd> eigvals(ATA, false);
	//double L   = eigvals.pseudoEigenvalueMatrix().diagonal().maxCoeff();
	double L     = 1056.56;

	MatrixXd CT  = C->adjoint();
	MatrixXd CCT = (*C)*(C->adjoint());
	int num_last = 5;
	MatrixXd tobj   = MatrixXd::Zero(num_last, 1);

	int k = 0;
	int j = 0;
	int iter_outer = 0;
	int	inneriter = 20;

	printf("start iterate\n");

	// The counter
	int checki = 1;					// the frequency of recording.
	int Num_i  = 0;
	double time_start = timing();	// starting point of time.

	// main loop
	while (nu>tolnu  &&  k<maxiter)
	{
		for(j = 0;j<inneriter;j++)  // inner loop
		{
			x_pre = x;
			//s = 1./(b-C*x);
			s = ((*b)-(*C)*x).cwiseInverse();

			//H = 2*(L + nu*s'*CCT*s);
			H = 2*(L + (nu*(s.adjoint())*CCT*s)(0,0));

			//x = proximalRegC(x-(A'*(A*x-y) + nu*C'*s)/H, d, mu/H, theta, regtype);
			x_temp = x-(AT*((*A)*x-(*Y)) + nu*CT*s)/H;
			proximalSCAD(&x, &x_temp, d, lambda/H, theta);

			if ((x-x_pre).squaredNorm() < 0.01*nu)
				break;
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
			MatrixXd terror = (*Y)-((*A)*x);
			err = normSquared(&terror);
			(*trace_obj)(iter_outer,0) = funSCAD(&x, d, lambda, theta) + 0.5*err;
			Num_i = 0;
		}
		printf("iter_outer, %d, obj cur, %f, time cost, %f\n",(int)(*iters)(iter_outer,0), (*trace_obj)(iter_outer,0), (*trace_time)(iter_outer,0));

		iter_outer += 1;

		// Stopping Criterion
		if (iter_outer> num_last)
		{
			double max_err=0,u;
			tobj = trace_obj->middleRows(iter_outer-num_last-1,num_last) - trace_obj->middleRows(iter_outer-num_last,num_last);
			for(int i=0;i<num_last;i++)
			{
				u = tobj(i,0)>0?tobj(i,0):-tobj(i,0);
				if(u>tobj(i,0))
					max_err= u;
			}
			if( max_err < tol)
				return 0;
		}
	}
	return 0;
}


