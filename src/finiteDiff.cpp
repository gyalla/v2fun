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

double Diff2(gsl_vector * x,double deltaEta,double bdry,int i)
{
  // i = xiCounter indices
  //if i < 5 use value at boundary. 5 for the index of xi.
  if (i<5)
    return (gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + bdry)/pow(deltaEta,2);
  else
    return ( gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + gsl_vector_get(x,i-5))/pow(deltaEta,2); //i+5 corresponds to i+1 for single terms. U_2 = xi_1+5 for example.
}

double Diff1(gsl_vector *x, double deltaEta,double bdry, int i)
{
  //same structure as above but for first derivative.
  if (i<5)
    return (gsl_vector_get(x,i+5)-bdry)/(2*deltaEta);
  else
    return (gsl_vector_get(x,i+5)-gsl_vector_get(x,i-5))/(2*deltaEta);
}

double Deriv2(gsl_vector * x,double deltaEta,double bdry,int i, Grid* grid)
{
  double xi = -1.0+(i/5)*deltaEta;
  return grid->d2XidY2(xi)*Diff1(x, deltaEta, bdry, i) +
      pow(grid->dXidY(xi),2)*Diff2(x, deltaEta, bdry, i);
}

double Deriv1(gsl_vector *x, double deltaEta,double bdry, int i, Grid* grid)
{
  double xi = -1.0+(i/5)*deltaEta;
  return Diff1(x, deltaEta, bdry, i)*grid->dXidY(xi);
}

double BdryDeriv2(gsl_vector *x, double deltaEta, int i, Grid* grid)
{
  //Using zero Neumann boundary condition we use ghost points for second derivative.
  double xi = -1.0+(i/5)*deltaEta;
  return (2*gsl_vector_get(x,i-5) - 2*gsl_vector_get(x,i)) *
      pow(grid->dXidY(xi)/deltaEta,2);
}

double BdryDiff2(gsl_vector * x ,double deltaEta,int i)
{
	//Using zero Neumann boundary condition we use ghost points for second derivative.
	double val = ( 2*gsl_vector_get(x,i-5) - 2*gsl_vector_get(x,i))/pow(deltaEta,2);
	return val; 
}

double Deriv1vT(gsl_vector * x,double deltaEta,int i, Grid* grid)
{
  // for vT we use normal finite difference approximation.
  double xi = -1.0+(i/5)*deltaEta;
  return (gsl_vector_get(x,i+1) - gsl_vector_get(x,i-1))/(2*deltaEta)*grid->dXidY(xi);
}

double Diff1vT(gsl_vector * x,double deltaEta,int i)
{
	// for vT we use normal finite difference approximation. 
	double val=  (gsl_vector_get(x,i+1) - gsl_vector_get(x,i-1))/(2*deltaEta);
	return val; 
}

