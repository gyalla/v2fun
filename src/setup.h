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

int SolveIC(gsl_vector* xi,double deltaEta, string file); 
int LinInterp(gsl_vector * Vec,double pt1,double pt2, double tempU1,double tempU2,double gridPt,int i );

int Solve4f0(gsl_vector * xi,constants * modelConst, double deltaEta);


#endif

