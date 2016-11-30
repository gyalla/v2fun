#ifndef SYSTEMSOLVE_H
#define SYSTEMSOLVE_H

#include<gsl/gsl_vector.h>

struct FParams{
	gsl_vector * XiN;
	double deltaT; 
	double deltaEta; 
	constants * modelConst;
}; 


int DeconstructXi(gsl_vector * xi, gsl_vector * U, gsl_vector * k,gsl_vector * ep,gsl_vector * v2, gsl_vector * f);

int ReconstructXi(gsl_vector * xi,gsl_vector * U,gsl_vector *k,gsl_vector * ep,gsl_vector *v2,gsl_vector * f);

int SysF(const gsl_vector * xi, void * p, gsl_vector * sysF);

int SetFirst2(gsl_vector * ep, gsl_vector * f,gsl_vector * k,  gsl_vector * v2, FParams * params,gsl_vector * sysF);

int SetUTerms(gsl_vector * U, gsl_vector * vT, FParams * params,gsl_vector * sysF);

int SetKTerms(gsl_vector * k,gsl_vector * P, gsl_vector * ep, gsl_vector * vT,FParams * params,gsl_vector *sysF);


int SetEpTerms(gsl_vector * ep, gsl_vector * P, gsl_vector * T, gsl_vector * vT, FParams * params, gsl_vector * sysF);
int SetV2Terms(gsl_vector * v2, gsl_vector * k, gsl_vector * f, gsl_vector * vT,gsl_vector * ep,  FParams * params, gsl_vector * sysF);


int SetFTerms(gsl_vector *f, gsl_vector * L, gsl_vector*P, gsl_vector * k, gsl_vector * T, gsl_vector * v2,  FParams * params, gsl_vector * sysF);

#endif
