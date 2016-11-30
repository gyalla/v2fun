#include<gsl/gsl_vector.h>
#include<math.h>

double Diff2(gsl_vector * x,double deltaEta,int i)
{
	double val =  ( gsl_vector_get(x,i+1) - 2*gsl_vector_get(x,i) + gsl_vector_get(x,i-1))/pow(deltaEta,2);
	return val;
}

double Diff1(gsl_vector *x, double deltaEta,int i) 
{
	double val = (gsl_vector_get(x,i+1)-gsl_vector_get(x,i-1))/(2*deltaEta);
	return val; 
}

double BdryDiff2(gsl_vector * x ,double deltaEta,int i)
{
	double val = ( 2*gsl_vector_get(x,i-1) - 2*gsl_vector_get(x,i))/pow(deltaEta,2);
	return val; 
}

