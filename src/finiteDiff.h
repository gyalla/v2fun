#ifndef FINITEDIFF_H
#define FINITEDIFF_H
#include<gsl/gsl_vector.h>
using namespace std;



double Diff1vT(gsl_vector * x,double deltaEta,int i);
double Diff2(gsl_vector * x,double deltaEta,double bdry, int i);

double Diff1(gsl_vector *x, double deltaEta,double bdry,int i);
double BdryDiff2(gsl_vector * x ,double deltaEta,int i);

#endif
