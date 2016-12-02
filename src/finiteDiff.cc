#include<iostream>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"computeTerms.h"
using namespace std;

double Diff2(gsl_vector * x,double deltaEta,double bdry,int i)
{

	double val;
	if (i<5)
		val =(gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + bdry)/pow(deltaEta,2);
	else
		val = ( gsl_vector_get(x,i+5) - 2*gsl_vector_get(x,i) + gsl_vector_get(x,i-5))/pow(deltaEta,2);
	return val;
}

double Diff1(gsl_vector *x, double deltaEta,double bdry, int i) 
{
	double val; 
	
	if (i<5)
		val = (gsl_vector_get(x,i+5)-bdry)/(2*deltaEta);
	else
		val = (gsl_vector_get(x,i+5)-gsl_vector_get(x,i-5))/(2*deltaEta);
	return val; 
}

double BdryDiff2(gsl_vector * x ,double deltaEta,int i)
{
	double val = ( 2*gsl_vector_get(x,i-5) - 2*gsl_vector_get(x,i))/pow(deltaEta,2);
	return val; 
}

double Diff1vT(gsl_vector * x,double deltaEta,int i)
{
	double val=  (gsl_vector_get(x,i+1) - gsl_vector_get(x,i-1))/(2*deltaEta);
	return val; 
}

