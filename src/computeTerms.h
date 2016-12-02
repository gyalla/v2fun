#ifndef COMPUTETERMS_H
#define COMPUTETERMS_H

#include<gsl/gsl_vector.h>
#include"setup.h"

int ComputeT(gsl_vector * xi, constants * modelConst,int i);
int ComputeL(gsl_vector * xi, constants * modelConst,int i);
int ComputeEddyVisc(gsl_vector * xi, constants * modelConst,int i);
int ComputeP(gsl_vector * xi,constants * modelConst,double deltaEta,int i );

double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta);
double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) ;
#endif
