#include<iostream>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"setup.h"

using namespace std;

int ComputeT(gsl_vector * k, gsl_vector * ep, constants * modelConst, gsl_vector * T)
{
	double firstTerm,secondTerm; 

	for (int i=0; (unsigned)i <T->size; i++)
	{
		firstTerm = (gsl_vector_get(k,i)/gsl_vector_get(ep,i));
		if (!isfinite(firstTerm))
		{
			cerr << "Error: non-finite (" << firstTerm << ")" << endl; 
			return 0; 
		}

		secondTerm = 6*sqrt(1/(modelConst->reyn*gsl_vector_get(ep,i)));
		if(!isfinite(secondTerm))
		{
			cerr << "Error: non-finite (" << secondTerm << ")" << endl;
			return 0; 
		}

		if (firstTerm >= secondTerm)
			 gsl_vector_set(T,i,firstTerm);
		else
			gsl_vector_set(T,i,secondTerm);

	}
	return 1; 
}

int ComputeL(gsl_vector * k,gsl_vector * ep,constants * modelConst,gsl_vector * L)
{
	double firstTerm,secondTerm; 
	for (int i=0; (unsigned)i<L->size;i++)
	{
		firstTerm = pow(gsl_vector_get(k,i),1.5)/gsl_vector_get(ep,i);
		if (!isfinite(firstTerm))
		{
			cerr << "Error: non-finite (" << firstTerm << ")" << endl;
			return 0; 
		}
		
		secondTerm = modelConst->Ceta*pow(1/(pow(modelConst->reyn,3)*gsl_vector_get(ep,i)),0.25);
		if (!isfinite(secondTerm))
		{
			cerr << "Error: non-finite (" << secondTerm << ")" << endl;
			return 0; 
		}

		if (firstTerm >= secondTerm)
			gsl_vector_set(L,i,modelConst->CL*firstTerm);
		else
			gsl_vector_set(L,i,modelConst->CL*secondTerm);
	}

	return 1;

}

int ComputeEddyVisc(gsl_vector * v2, gsl_vector * T, constants * modelConst, gsl_vector * vT)
{
	double val; 
	for (int i = 0; (unsigned)i<vT->size;i++)
	{
		val = modelConst->Cmu*gsl_vector_get(v2,i)*gsl_vector_get(T,i);
		if (!isfinite(val))
		{
			cerr << "Error: non-finite (" << val << ")" << endl;
			return 1; 
		}
		gsl_vector_set(vT,i,val);
	}
	return 0; 
}

int ComputeP(gsl_vector * U,gsl_vector *vT,double deltaEta,gsl_vector *P)
{
	double val; 
	
	for (int i = 0; (unsigned)i< (P->size-1);i++)
	{
		if (i==0)
			val = gsl_vector_get(vT,i)*(gsl_vector_get(U,i+1)/(deltaEta)); //must use forward difference approximation
		else
			val = gsl_vector_get(vT,i)*( (gsl_vector_get(U,i+1)-gsl_vector_get(U,i-1))/(2*deltaEta));
		if (!isfinite(val))
		{
			cerr << "Error: non-finite (" << val << ")" << endl;
			return 0; 
		}
		gsl_vector_set(P,i,val);
	}
	return 1; 
}


