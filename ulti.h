//============================================================================
// Name        : ABCDEF_distributed
// Author      : linboqiao
// Version     :
// Copyright   : Your copyright notice
// Description : subfunctions
//============================================================================


#include "Eigen/Dense"
using namespace Eigen;

double timing();
long mymin(double *x, long n);
double min(double a,double b);
double max(double a,double b);
void proximalCapL1(double *x, double *d, long n, double lambda, double theta);
void proximalLSP(double *x, double *d, long n, double lambda, double theta);
//void proximalSCAD(double *x, double *d, long n, double lambda, double theta);
void proximalSCAD_d(double *x, double *d, int n, double lambda, double theta);
void proximalSCAD(MatrixXd *x, MatrixXd *d, int n, double lambda, double theta);
void proximalMCP(double *x, double *d, long n, double lambda, double theta);
void updateX(double *d,long n,double lambda,double theta,int type,double *x);
double funCapL1( double *x, long n, double lambda, double theta);
double funLSP(double *x, long n, double lambda, double theta);
//double funSCAD(double *x, long n, double lambda, double theta);
double funSCAD(MatrixXd *x, int n, double lambda, double theta);
double funMCP( double *x, long n, double lambda, double theta);
double computeObj(double *x,long n,double lambda,double theta,int type);
void matrixMult(double *MA,double *MB,double *MC,int m,int n,int p);
void transpose(double *MA, double *MB, int m,int n);
void computeL(double *MA,double *ML,int m1,int n1);
void assignArray(double *MA,double *MB,int m,int n);
double normSquared(MatrixXd *var);

double fRand(double fMin, double fMax);
int random_init(int d, int n, int m, \
		MatrixXd *eX, MatrixXd *eC, MatrixXd *eY, MatrixXd *eb, MatrixXd *ebetaTrue, double *eL);

int read_from_csv(int d, int n, int m, \
		MatrixXd *eX, MatrixXd *eC, MatrixXd *eY, MatrixXd *eb, MatrixXd *ebetaTrue, double *eL);

int write_to_csv(char* resultFile, \
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time);


//ABCDEF batch
int ABCDEF(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads, char *funcname);

//ABCDEF coordinate cyclic
int ABCDEF_coord_cyclic(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads,char *funcname);

//ABCDEF coordinate cyclic distributed
int ABCDEF_coord_cyclic_distributed(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads,char *funcname);

//ABCDEF coordinate randomized distributed
int ABCDEF_coord_rand_distributed(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads, char *funcname);

//ABCDEF coordinate cyclic distributed asynchronous
int ABCDEF_coord_cyclic_distributed_asynchronous(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads,char *funcname);

//ABCDEF coordinate randomized distributed asynchronous
int ABCDEF_coord_rand_distributed_asynchronous(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads, char *funcname);

//ABCDEF coordinate cyclic asynchronous
int ABCDEF_coord_cyclic_asynchronous(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads,char *funcname);

//ABCDEF coordinate randomized asynchronous
int ABCDEF_coord_rand_asynchronous(MatrixXd *A, MatrixXd *C, MatrixXd *Y, MatrixXd *b,\
		double lambda, int regtype, double theta, int maxiter,\
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time,\
		int m, int n, int d, int n_threads, char *funcname);
