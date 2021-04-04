//============================================================================
// Name        : ABCDEF_distributed
// Author      : linboqiao
// Version     :
// Copyright   : Your copyright notice
// Description : subfunctions
//============================================================================


#include <fstream>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <sys/time.h>
using namespace std;

#include "ulti.h"

double max(double a,double b) {return (a>b)?a:b;}
double min(double a,double b) {return (a<b)?a:b;}

double timing(){
	double time;
	struct timeval timmer;

	gettimeofday(&timmer,NULL);
	time = timmer.tv_sec + timmer.tv_usec*1e-6;
	return time;
}


long mymin(double *x, long n)
{
	long i,ind;
	double temp;
	temp = x[0];
	ind = 0;
	for (i=1;i<n;i++)
	{
		if (x[i] < temp)
		{
			ind = i;
			temp = x[i];
		}

	}
	return ind;
}


void proximalCapL1(double *x, double *d, long n, double lambda, double theta)
{
	long i;
	double u,x1,x2;
	for(i=0;i<n;i++)
	{   
		u = fabs(d[i]);
		x1 = max(u,theta); 
		x2 = min(theta, max(0, u - lambda));
		if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + lambda*(theta - x2) < 0)
			x[i] = x1;
		else
			x[i] = x2;
		x[i] = d[i]>= 0 ? x[i] : -x[i];

	}	
	return;
}


void proximalLSP(double *x, double *d, long n, double lambda, double theta)
{
	long i;
	double v,u,z,sqrtv;
	double xtemp[3],ytemp[3];
	xtemp[0] = 0.0;

	for(i=0;i<n;i++)
	{ 
		u = fabs(d[i]);
		z = u - theta;
		v = z*z - 4.0*(lambda - u*theta);

		if (v < 0)
			x[i] = 0;
		else
		{
			sqrtv = sqrt(v);
			xtemp[1] = max(0,0.5*(z + sqrtv));
			xtemp[2] = max(0,0.5*(z - sqrtv));

			ytemp[0] = 0.5*u*u;
			ytemp[1] = 0.5*(xtemp[1] - u)*(xtemp[1] - u) + lambda*log(1.0 + xtemp[1]/theta);
			ytemp[2] = 0.5*(xtemp[2] - u)*(xtemp[2] - u) + lambda*log(1.0 + xtemp[2]/theta);

			x[i] = xtemp[mymin(ytemp,3)];
			x[i] = d[i]>= 0 ? x[i] : -x[i];
		}
	}	
	return;
}


void proximalMCP(double *x, double *d, long n, double lambda, double theta)
{
	long i;
	double x1,x2,v,u,z;
	z = theta*lambda;
	if (theta  >  1)
	{
		for(i=0;i<n;i++)
		{ 
			u = fabs(d[i]);
			x1 = min(z,max(0,theta*(u - lambda)/(theta - 1.0)));
			x2 = max(z,u);
			if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + x1*(lambda - 0.5*x1/theta) - 0.5*z*lambda < 0)
				x[i] = x1;
			else
				x[i] = x2;
			x[i] = d[i]>= 0 ? x[i] : -x[i];
		}
	}
	else if (theta < 1)
	{
		for(i=0;i<n;i++)
		{ 
			u = fabs(d[i]);
			v = theta*(u - lambda)/(theta - 1);
			x1 = fabs(v) > fabs(v - z) ? 0.0 : z;
			x2 = max(z,u);
			if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + x1*(lambda - 0.5*x1/theta) - 0.5*z*lambda < 0)
				x[i] = x1;
			else
				x[i] = x2;
			x[i] = d[i]>= 0 ? x[i] : -x[i];
		}
	}
	else
	{
		for (i=0;i<n;i++)
		{
			u = fabs(d[i]);
			x1 = lambda > u ? 0.0 : z;
			x2 = max(z,u);
			if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + x1*(lambda - 0.5*x1/theta) - 0.5*z*lambda < 0)
				x[i] = x1;
			else
				x[i] = x2;
			x[i] = d[i]>= 0 ? x[i] : -x[i];	
		}
	}
	return;
}


void updateX(double *d,long n,double lambda,double theta,int type,double *x)
{
	switch (type)
	{
	case 1:
		proximalCapL1(x, d, n, lambda, theta);
		break;
	case 2:
		proximalLSP(x, d, n, lambda, theta);
		break;
	case 3:
		//proximalSCAD(x, d, n, lambda, theta);
		break;
	case 4:
		proximalMCP(x, d, n, lambda, theta);
		break;
	default:
		proximalCapL1(x, d, n, lambda, theta);
	}
}


double funCapL1( double *x, long n, double lambda, double theta)
{
	long i;
	double u = 0.0;;
	for(i=0;i<n;i++)
	{   
		u += min(fabs(x[i]),theta);
	}
	double f = u*lambda;
	return f;
}


double funLSP(double *x, long n, double lambda, double theta)
{
	long i;
	double u = 0.0;

	for(i=0;i<n;i++)
	{ 
		u += log(1.0 + fabs(x[i])/theta);
	}
	double f = u*lambda;
	return f;
}

/*
double funSCAD(double *x, long n, double lambda, double theta)
{
    long i;
	double u,v,y,z,w;
	y = theta*lambda;
	w = lambda*lambda;
	z = 0.5*(theta+1.0)*w;

	u = 0.0;
	for(i=0;i<n;i++)
	{ 
		v = fabs(x[i]);
		if (v <= lambda) 
			u += lambda*v;
		else if (v > y)
			u += z;
		else
			u += 0.5*(v*(2*y - v) - w)/(theta-1.0);
	}	

	return u;
}
 */


double funMCP( double *x, long n, double lambda, double theta)
{
	long i;
	double v,u,y;
	y = theta*lambda;
	u = 0.0;
	for(i=0;i<n;i++)
	{
		v = fabs(x[i]);
		if (v <= y)
			u += v*(lambda - 0.5*v/theta);
		else
			u += 0.5*y*lambda;
	}

	return u;
}


double computeObj(double *x,long n,double lambda,double theta,int type)
{
	double f = 0;
	switch (type)
	{
	case 1:
		f = funCapL1(x, n, lambda, theta);
		return f;
		break;
	case 2:
		f = funLSP(x, n, lambda, theta);
		return f;
		break;
	case 3:
		//f = funSCAD(x, n, lambda, theta);
		return f;
		break;
	case 4:
		f = funMCP(x, n, lambda, theta);
		return f;
		break;
	default:
		f = funCapL1(x, n, lambda, theta);
		return f;
	}
}


void matrixMult(double *MA,double *MB,double *MC,int m,int n,int p)
{
	/* MA:mxn  MB:nxp   MC:mxp */
	for(int i=0;i<m;i++)
		for(int j=0;j<p;j++){
			MC[i+j*m] = 0.0;
			for(int k=0;k<n;k++)
				MC[i+j*m] += MA[i+k*m]*MB[k+j*n];
		}
	return;
}


void transpose(double *MA, double *MB, int m,int n)
{
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			MB[i+j*n] = MA[j+i*m];
	return;
}


void computeL(double *MA,double *ML,int m1,int n1)
{
	double *MAT = (double *)malloc(m1*n1*sizeof(double));
	double *MATA = (double *)malloc(n1*n1*sizeof(double));
	transpose(MA,MAT,m1,n1);
	matrixMult(MAT,MA,MATA,n1,m1,n1);
	for(int i=0;i<n1;i++)
		ML[i] = MATA[i+i*n1];
	free(MAT);
	free(MATA);
}


void assignArray(double *MA,double *MB,int m,int n)
{
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			MB[i+j*m] = MA[i+j*m];
	return;
}


void proximalSCAD(MatrixXd *x, MatrixXd *d, int n, double lambda, double theta)
{
	long i;
	double u,z,w,td;
	double xtemp[3],ytemp[3];
	z = theta*lambda;
	w = lambda*lambda;
	for(i=0;i<n;i++)
	{
		td = (*d)(i,0);
		u = fabs(td);

		xtemp[0] = min(lambda,max(0,u - lambda));
		xtemp[1] = min(z,max(lambda,(u*(theta-1.0)-z)/(theta-2.0)));
		xtemp[2] = max(z,u);

		ytemp[0] = 0.5*(xtemp[0] - u)*(xtemp[0] - u) + lambda*xtemp[0];
		ytemp[1] = 0.5*(xtemp[1] - u)*(xtemp[1] - u) + 0.5*(xtemp[1]*(-xtemp[1] + 2*z) - w)/(theta-1.0);
		ytemp[2] = 0.5*(xtemp[2] - u)*(xtemp[2] - u) + 0.5*(theta+1.0)*w;

		(*x)(i,0) = xtemp[mymin(ytemp,3)];
		(*x)(i,0) = (*d)(i,0)>= 0 ? ((*x)(i,0)) : -((*x)(i,0));
	}
	return;
}


void proximalSCAD_d(double *x, double *d, int n, double lambda, double theta)
{
	long i;
	double u,z,w;
	double xtemp[3],ytemp[3];
	z = theta*lambda;
	w = lambda*lambda;
	for(i=0;i<n;i++)
	{
		u = fabs(d[i]);
		xtemp[0] = min(lambda,max(0,u - lambda));
		xtemp[1] = min(z,max(lambda,(u*(theta-1.0)-z)/(theta-2.0)));
		xtemp[2] = max(z,u);

		ytemp[0] = 0.5*(xtemp[0] - u)*(xtemp[0] - u) + lambda*xtemp[0];
		ytemp[1] = 0.5*(xtemp[1] - u)*(xtemp[1] - u) + 0.5*(xtemp[1]*(-xtemp[1] + 2*z) - w)/(theta-1.0);
		ytemp[2] = 0.5*(xtemp[2] - u)*(xtemp[2] - u) + 0.5*(theta+1.0)*w;

		x[i] = xtemp[mymin(ytemp,3)];
		x[i] = d[i]>= 0 ? x[i] : -x[i];

	}
	return;
}


double funSCAD(MatrixXd *x, int n, double lambda, double theta)
{
	long i;
	double u,v,y,z,w,tx;
	y = theta*lambda;
	w = lambda*lambda;
	z = 0.5*(theta+1.0)*w;

	u = 0.0;
	for(i=0;i<n;i++)
	{
		tx = (*x)(i,0);
		v = fabs(tx);
		if (v <= lambda)
			u += lambda*v;
		else if (v > y)
			u += z;
		else
			u += 0.5*(v*(2*y - v) - w)/(theta-1.0);
	}
	return u;
}

double normSquared(MatrixXd *var)
{
	double r=0;
	for(int i=0;i<var->rows();i++)
		for(int j=0;j<var->cols();j++)
			r += (*var)(i,j)*(*var)(i,j);
	return r;

}

double fRand(double fMin=-1, double fMax=1)
{
    return  ((double)(fMin)) +  ((double)(fMax - fMin))*((double)rand())/((double)RAND_MAX);
}

int random_init(int d, int n, int m, \
		MatrixXd *eX, MatrixXd *eC, MatrixXd *eY, MatrixXd *eb, MatrixXd *ebetaTrue, double *eL)
{
	unsigned int tseed = 5;
	std::srand(tseed); /* Seed = 5 */
	
	/*

		X         = randn(n, d);                      % generate X.
		X         = X./repmat(sqrt(sum(X.*X)), n,1);  % column normalization.
		betaTrue  = [ones(5,1); zeros(d-5,1)];        % true solution.
		sigma     = 0.1;                              % the level of noise
		Y         = X*betaTrue + sigma^2*rand(n, 1);  % generate Y.
		C         = randn(m, d);                      % generate C.
		b         = C*betaTrue + sigma*rand(m, 1);        % generate b.
		L         = eigs(X'*X, 1);
	*/

	// X Matrix init
	for(int j=0; j < n; j++){
		for(int i=0; i < d; i++){
			(*eX)(j, i) = fRand();
		}
	}
	
	// betaTrue init 
	for(int j=0; j < 5; j++){
		(*ebetaTrue)(j, 0) = 1;
	}
	for(int j=5; j < d; j++){
		(*ebetaTrue)(j, 0) = 0;
	}

	// Y Matrix init
	for(int j=0; j < n; j++){
		double temp = 0;		
		for(int i=0; i < d; i++){
			temp +=(*eX)(j, i)*(*ebetaTrue)(i, 0);
		}
		temp += ((double)rand())/((double)RAND_MAX);
		(*eY)(j, 0) = temp;
	}

	// C Matrix init
	for(int j=0; j < m; j++){
		for(int i=0; i < d; i++){
			(*eC)(j, i) = fRand();
		}
	}
		
	// b init	
	for(int j=0; j < m; j++){
		double temp = 0;		
		for(int i=0; i < d; i++){
			temp +=(*eC)(j, i)*(*ebetaTrue)(i, 0);
		}
		temp += ((double)rand())/((double)RAND_MAX);
		(*eb)(j, 0) = temp;
	}

	// L set
	/*
	cout << "eigenvalues calculating..." << endl;
	MatrixXd XTX = eX->transpose()*(*eX);
	VectorXcd eivals = XTX.eigenvalues();
	max_eigen = eivals.array().max(-999999);
	cout << eivals << "\t" << max_eigen << endl;
	*eL = eivals[0];
	*/
	
	// *eL = fRand();
	// *eL = -211.746863604871;
	*eL = 999;
	
	return 0;
}



int read_from_csv(int d, int n, int m, \
		MatrixXd *eX, MatrixXd *eC, MatrixXd *eY, MatrixXd *eb, MatrixXd *ebetaTrue, double *eL)
{
	char filename[255] = {0};
	char X_file[255] = {0};	
	char C_file[255] = {0};
	char Y_file[255] = {0};
	char b_file[255] = {0};
	char betaTrue_file[255] = {0};
	char L_file[255] = {0};	
	sprintf(filename,"data_d%d_n%d_m%d", d, n, m);
	sprintf(X_file,"./csv/data_d%d_n%d_m%d_X.csv", d, n, m);
	sprintf(C_file,"./csv/data_d%d_n%d_m%d_C.csv", d, n, m);
	sprintf(Y_file,"./csv/data_d%d_n%d_m%d_Y.csv", d, n, m);
	sprintf(b_file,"./csv/data_d%d_n%d_m%d_b.csv", d, n, m);
	sprintf(betaTrue_file,"./csv/data_d%d_n%d_m%d_betaTrue.csv", d, n, m);
	sprintf(L_file,"./csv/data_d%d_n%d_m%d_L.csv", d, n, m);
	//printf("filename %s, %s, %s, %s, %s, %s, %s\n",filename, X_file, C_file, Y_file, b_file, betaTrue_file, L_file);
	
	string line;
	
	// read X file 
	ifstream X_in(X_file, ios::in);
	int x_row_count = 0;
	while(getline(X_in, line)){
		string number;
		istringstream readstr(line);
		for(int i=0; i < d; i++){
			getline(readstr, number, ',');
			(*eX)(x_row_count, i) = atof(number.c_str());
		}
		x_row_count++;
	}

	// read C file
	ifstream C_in(C_file, ios::in);
	int c_row_count = 0;
	while(getline(C_in, line)){
		string number;
		istringstream readstr(line);
		for(int i=0; i<d; i++){
			getline(readstr, number, ',');
			(*eC)(c_row_count, i) = atof(number.c_str());
		}
		c_row_count++;
	}

	// read Y file 
	ifstream Y_in(Y_file, ios::in);
	int y_row_count = 0;
	while(getline(Y_in, line)){
		(*eY)(y_row_count, 0) = atof(line.c_str());
		y_row_count++;
	}
		
	//read b file 
	ifstream b_in(b_file, ios::in);
	int b_row_count = 0;
	while(getline(b_in, line)){
		(*eb)(b_row_count, 0) = atof(line.c_str());
		b_row_count++;
	}

	//read betaTrue file 
	ifstream beta_in(betaTrue_file, ios::in);
	int beta_row_count=0;
	while(getline(beta_in, line)){
		(*ebetaTrue)(beta_row_count, 0) = atoi(line.c_str());
		beta_row_count++;
	}

	ifstream L_in(L_file, ios::in);
	getline(L_in, line);
	*eL = atoi(line.c_str()); 
	
	//pmatFile = matOpen(filename,"r");
	return 0;

}

int write_to_csv(char* resultFile, \
		MatrixXd *iters, MatrixXd *trace_beta, MatrixXd *trace_obj, MatrixXd *trace_time)
{
	string iter = "iteration";
	string beta = "trace_beta";
	string object = "trace_obj";
	string time = "trace_time";
	
	ofstream outFile;
	outFile.open(resultFile, ios::out);
	outFile<<iter<<beta<<object<<time<<"\n";
	
	outFile<<*iters<<*trace_beta<<*trace_obj<<*trace_time<<"\n";

	outFile.close();
	
	
	return 0;

}

