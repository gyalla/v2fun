#include<iostream>
#include<gsl/gsl_vector.h>
#include<math.h>

using namespace std;

int ComputeT(gsl_vector * k, gsl_vector * ep, int reyn, gsl_vector * T)
{
	double firstTerm,secondTerm; 

	for (int i=0; i <T->size; i++)
	{
		firstTerm = (gsl_vector_get(k,i)/gsl_vector_get(ep,i));
		if (!isfinite(firstTerm))
		{
			cerr << "Error: non-finite (" << firstTerm << ")" << endl; 
			return 1; 
		}

		secondTerm = 6*sqrt(1/(reyn*gsl_vector_get(ep,i)));
		if(!isfinite(secondTerm))
		{
			cerr << "Error: non-finite (" << secondTerm << ")" << endl;
			return 1; 
		}

		if (firstTerm >= secondTerm)
			 gsl_vector_set(T,i,firstTerm);
		else
			gsl_vector_set(T,i,secondTerm);

	}
	return 0; 
}

int ComputeL(gsl_vector * k,gsl_vector * ep,int reyn, double CL,double Ceta,gsl_vector * L)
{
	double firstTerm,secondTerm; 
	for (int i=0; i<L->size;i++)
	{
		firstTerm = pow(gsl_vector_get(k,i),1.5)/gsl_vector_get(ep,i);
		if (!isfinite(firstTerm))
		{
			cerr << "Error: non-finite (" << firstTerm << ")" << endl;
			return 1; 
		}
		
		secondTerm = Ceta*pow(1/(pow(reyn,3)*gsl_vector_get(ep,i)),0.25);
		if (!isfinite(secondTerm))
		{
			cerr << "Error: non-finite (" << secondTerm << ")" << endl;
			return 1; 
		}

		if (firstTerm >= secondTerm)
			gsl_vector_set(L,i,CL*firstTerm);
		else
			gsl_vector_set(L,i,CL*secondTerm);
	}

	return 0;

}

int ComputeEddyVisc(gsl_vector * v2, gsl_vector * T, double Cmu, gsl_vector * vT)
{
	double val; 
	for (int i = 0; i<vT->size;i++)
	{
		val = Cmu*gsl_vector_get(v2,i)*gsl_vector_get(T,i);
		if (!isfinite(val))
		{
			cerr << "Error: non-finite (" << val << ")" << endl;
			return 1; 
		}
		gsl_vector_set(vT,i,val);
	}
	return 0; 
}



