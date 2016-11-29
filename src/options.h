#ifndef OPTIONS_H
#define OPTIONS_H

#include<gsl/gsl_vector.h>
using namespace std;

int Grvy_Input_Parse(double &reyn, double & Cmu, double & C1, double & C2, double & sigmaEp, double & CL, double &Cep2, double &Cep1, double &Ceta,string & filename);

void log(int verbose,int level, string msg);

int SolveIC(gsl_vector * U, gsl_vector * k, gsl_vector * ep, gsl_vector * v2,double deltaEta, string file); 
int LinInterp(gsl_vector * Vec,double pt1,double pt2, double tempU1,double tempU2,double gridPt,int i );

#endif

