#include<iostream>
#include<gsl/gsl_vector.h>
#include<math.h>
#include"setup.h"
#include"finiteDiff.h"

using namespace std;

double ComputeT(gsl_vector * xi, constants * modelConst,int i)
{
	double firstTerm,secondTerm; 
	double xiCounter = 5*(i-1); 	
	

	firstTerm = (gsl_vector_get(xi,xiCounter+1)/gsl_vector_get(xi,xiCounter+2));
	if (!isfinite(firstTerm))
	{
		cerr << "Error: T non-finite (" << firstTerm << ")" << endl; 
		return -1; 
	}

	secondTerm = 6*sqrt(1/(modelConst->reyn*gsl_vector_get(xi,xiCounter+2)));
	if(!isfinite(secondTerm))
	{
		cerr << "Error: T non-finite (" << secondTerm << ")" << endl;
		return -1; 
	}

	if (firstTerm >= secondTerm)
			return firstTerm;
	else
			return secondTerm;
	}


double ComputeL(gsl_vector * xi,constants * modelConst,int i)
{
	double firstTerm,secondTerm; 
	double xiCounter = 5*(i-1); 

	firstTerm = pow(gsl_vector_get(xi,xiCounter+1),1.5)/gsl_vector_get(xi,xiCounter+2);
	if (!isfinite(firstTerm))
	{
		cerr << "Error: L non-finite (" << firstTerm << ")" << endl;
		return -1; 
	}
	
	secondTerm = modelConst->Ceta*pow(1/(pow(modelConst->reyn,3)*gsl_vector_get(xi,xiCounter+2)),0.25);
	if (!isfinite(secondTerm))
	{
		cerr << "Error: L non-finite (" << secondTerm << ")" << endl;
		return -1; 
	}

	if (firstTerm >= secondTerm)
		return modelConst->CL*firstTerm;
	else
		return modelConst->CL*secondTerm;

}

double ComputeEddyVisc(gsl_vector * xi, constants * modelConst,int i)
{
	double val; 
	double xiCounter = 5*(i-1);
	val = modelConst->Cmu*gsl_vector_get(xi,xiCounter+3)*ComputeT(xi,modelConst,i);
	if (!isfinite(val))
	{
		cerr << "Error: vT non-finite (" << val << ")" << endl;
		return -1; 
	}
	return val; 
}

double ComputeP(gsl_vector * xi,constants * modelConst,double deltaEta,int i)
{
	double val; 
	double xiCounter = 5*(i-1);

	if (i==0)
	{
		cerr << "Error: i = 0" << endl; 
		return -1; 
	}
	val = ComputeEddyVisc(xi,modelConst,i)*pow(Diff1(xi,deltaEta,0,i),2);
	if (!isfinite(val))
	{
		cerr << "Error: P non-finite (" << val << ")" << endl;
		return -1; 
	}
	return val; 
}

double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) 
{
	double ep0 = ((2*gsl_vector_get(xi,1))/(modelConst->reyn*pow(deltaEta,2)));
	if (!isfinite(ep0))
	{
		cerr << "Error: ep0 non-finite (" << ep0 << ")" << endl; 
		return -1; 
	}
	return ep0; 
}

double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta)
{
	double f0  = (( (20*gsl_vector_get(xi,3))/( pow(modelConst->reyn,3)*ComputeEp0(xi,modelConst,deltaEta)*pow(deltaEta,4))));
	if(!isfinite(f0))
	{
		cerr << "Error: f0 non-finite (" << f0 << ")" << endl; 
		return -1; 
	}
	return f0; 
}
