#ifndef FINITEDIFF_H
#define FINITEDIFF_H
#include<gsl/gsl_vector.h>


double Diff2(gsl_vector * x,double deltaEta,int i);

double Diff1(gsl_vector *x, double deltaEta,int i);
double BdryDiff2(gsl_vector * x ,double deltaEta,int i);

#endif
