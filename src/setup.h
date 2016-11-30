#ifndef OPTIONS_H
#define OPTIONS_H

#include<gsl/gsl_vector.h>
using namespace std;

struct constants {
	double reyn;
	double Cmu;
	double C1;
	double C2;
	double Cep1;
	double Cep2;
	double Ceta;
	double CL;
	double sigmaEp;
};

int Grvy_Input_Parse(constants * modelConst,string &filename);

void log(int verbose,int level, string msg);

int SolveIC(gsl_vector * U, gsl_vector * k, gsl_vector * ep, gsl_vector * v2,double deltaEta, string file); 
int LinInterp(gsl_vector * Vec,double pt1,double pt2, double tempU1,double tempU2,double gridPt,int i );

int Solve4f0(gsl_vector * k, gsl_vector * ep, gsl_vector * v2, gsl_vector * P, gsl_vector *  T,gsl_vector * L,constants * modelConst, double deltaEta, gsl_vector * f);


#endif

