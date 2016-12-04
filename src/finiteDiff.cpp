//-------------------------------------------------- 
// finiteDiff: Defines finite difference approximations. 
//
// 12/3/2016 - (gry88) Written for CSE380 final project. 
//--------------------------------------------------
#include<iostream>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"computeTerms.h"
#include"setup.h"
using namespace std;

double Diff2(gsl_vector * x,double deltaEta,double bdry,int i)
{

	double val;
	//if i < 5 use value at boundary. 5 for the index of xi. 
	if (i<5)
		val =(gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + bdry)/pow(deltaEta,2);
	else
		val = ( gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + gsl_vector_get(x,i-5))/pow(deltaEta,2); //i+5 corresponds to i+1 for single terms. U_2 = xi_1+5 for example. 
	return val;
}

double Diff1(gsl_vector *x, double deltaEta,double bdry, int i) 
{
	//same structure as above but for first derivative. 
	double val; 
	
	if (i<5)
		val = (gsl_vector_get(x,i+5)-bdry)/(2*deltaEta);
	else
		val = (gsl_vector_get(x,i+5)-gsl_vector_get(x,i-5))/(2*deltaEta);
	return val; 
}

double BdryDiff2(gsl_vector * x ,double deltaEta,int i)
{
	//Using zero Nuemann boundary condition we use ghost points for second derivative.  
	double val = ( 2*gsl_vector_get(x,i-5) - 2*gsl_vector_get(x,i))/pow(deltaEta,2);
	return val; 
}

double Diff1vT(gsl_vector * x,double deltaEta,int i)
{
	// for vT we use normal finite difference approximation. 
	double val=  (gsl_vector_get(x,i+1) - gsl_vector_get(x,i-1))/(2*deltaEta);
	return val; 
}

