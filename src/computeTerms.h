#ifndef COMPUTETERMS_H
#define COMPUTETERMS_H

#include<gsl/gsl_vector.h>

int ComputeT(gsl_vector * k, gsl_vector * ep, int reyn, gsl_vector * T);
int ComputeL(gsl_vector * k,gsl_vector * ep,int reyn, double CL,double Ceta,gsl_vector * L);
int ComputeEddyVisc(gsl_vector * v2, gsl_vector * T, double Cmu, gsl_vector * vT);
int ComputeP(gsl_vector * U,gsl_vector *vT,double deltaEta,gsl_vector *P);

#endif
