//-------------------------------------------------- 
// finiteDiff: Defines finite difference approximations. 
//
// 12/3/2016 - (gry88) Written for CSE380 final project. 
//--------------------------------------------------
#include "finiteDiff.h"
#include<gsl/gsl_vector.h>
#define _USE_MATH_DEFINES // Needed for M_PI
#include<math.h>
#include <iostream>
using namespace std;

inline double Diff2(gsl_vector * x,double deltaEta,double bdry,int i)
{
  // i = xiCounter indices
  //if i < 5 use value at boundary. 5 for the index of xi.
  if (i<5)
    return (gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + bdry)/pow(deltaEta,2);
  else
    return ( gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + gsl_vector_get(x,i-5))/pow(deltaEta,2); //i+5 corresponds to i+1 for single terms. U_2 = xi_1+5 for example.
}

inline double Diff1(gsl_vector *x, double deltaEta,double bdry, int i)
{
  //same structure as above but for first derivative.
  if (i<5)
    return (gsl_vector_get(x,i+5)-bdry)/(2*deltaEta);
  else
    return (gsl_vector_get(x,i+5)-gsl_vector_get(x,i-5))/(2*deltaEta);
}

double Deriv2(gsl_vector * x, double bdry, int i, Grid* grid)
{
  // Since xi is uniformly spaced, \xi(0) = \Delta \xi
  double delta = gsl_vector_get(grid->chi, 0);
  double chi = gsl_vector_get(grid->chi, (i/5));
  return grid->d2ChidY2(chi)*Diff1(x, delta, bdry, i) +
      pow(grid->dChidY(chi),2)*Diff2(x, delta, bdry, i);
}

double Deriv1(gsl_vector *x, double bdry, int i, Grid* grid)
{
  // Since xi is uniformly spaced, \xi(0) = \Delta \xi
  double delta = gsl_vector_get(grid->chi, 0);
  double chi = gsl_vector_get(grid->chi, (i/5));
  return Diff1(x, delta, bdry, i)*grid->dChidY(chi);
}

double BdryDeriv2(gsl_vector *x, int i, Grid* grid)
{
  //Using zero Neumann boundary condition we use ghost points for second derivative.
  double delta = gsl_vector_get(grid->chi, 0);
  return BdryDiff2(x,delta,i);
}

double BdryDiff2(gsl_vector * x ,double deltaEta,int i)
{
	//Using zero Neumann boundary condition we use ghost points for second derivative.
	double val = ( 2*gsl_vector_get(x,i-5) - 2*gsl_vector_get(x,i))/pow(deltaEta,2);
	return val; 
}

double Deriv1vT(gsl_vector * x, int i, Grid* grid)
{
  // for vT we use normal finite difference approximation.
  double delta = gsl_vector_get(grid->chi, 0);
  double chi = gsl_vector_get(grid->chi, (i/5));
  return Diff1vT(x,delta,i)*grid->dChidY(chi);
}

double Diff1vT(gsl_vector * x,double deltaEta,int i)
{
	// for vT we use normal finite difference approximation. 
	double val=  (gsl_vector_get(x,i+1) - gsl_vector_get(x,i-1))/(2*deltaEta);
	return val; 
}

