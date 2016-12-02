#ifndef COMPUTETERMS_H
#define COMPUTETERMS_H

#include<gsl/gsl_vector.h>
#include"setup.h"
using namespace std;

double ComputeT(gsl_vector * xi, constants * modelConst,int i);
double ComputeL(gsl_vector * xi, constants * modelConst,int i);
double ComputeEddyVisc(gsl_vector * xi, constants * modelConst,int i);
double ComputeP(gsl_vector * xi,constants * modelConst,double deltaEta,int i );

double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta);
double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) ;
#endif
