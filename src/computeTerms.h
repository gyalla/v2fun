/**
 * \file 
 */
#ifndef COMPUTETERMS_H
#define COMPUTETERMS_H

#include<gsl/gsl_vector.h>
#include"setup.h"
using namespace std;

double ComputeT(gsl_vector * xi, constants * modelConst,int i);
/**
 * \brief Compute Turbulent length scale, L. 
 *
 * Compute L using \f[ C_L \max\{\frac{k}{\epsilon}\} \f]
 * \param xi vector of unknowns, U,k,ep,v2,f.
 * \param modelConst model constants
 * \param i position at which to compute L.
 * \return Length scale L at i. 
 */
double ComputeL(gsl_vector * xi, constants * modelConst,int i);
double ComputeEddyVisc(gsl_vector * xi, constants * modelConst,int i);
double ComputeP(gsl_vector * xi,constants * modelConst,double deltaEta,int i );

double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta);
double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) ;
#endif
