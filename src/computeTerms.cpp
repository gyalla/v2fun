//--------------------------------------------------
// computeTerms: Compute terms T,L,P,vT and f(0),ep(0).
// 
// 12/3/2016 - (gry88) Writen for final project CSE380.  
//-------------------------------------------------- 
#include<gsl/gsl_vector.h>
#include<math.h>
#include"setup.h"
#include"finiteDiff.h"
#include"loglevel.h"
using namespace std;

double ComputeT(gsl_vector * xi, constants * modelConst,int i)
{
	double firstTerm,secondTerm; //1st and 2nd term as in documentation. 
	double xiCounter = 5*(i-1);  //counter relative to xi. 	
	
	Log(logDEBUG1) << "Computing T";	
	firstTerm = (gsl_vector_get(xi,xiCounter+1)/gsl_vector_get(xi,xiCounter+2));
	if (!isfinite(firstTerm))
	{
		Log(logERROR) << "Error: T non-finite (" << firstTerm << ")";
		Log(logERROR) << "-Note ep = " << gsl_vector_get(xi,xiCounter+2) << " at " << i;
		exit(1);
	}

	secondTerm = 6*sqrt(1/(modelConst->reyn*gsl_vector_get(xi,xiCounter+2)));
	if(!isfinite(secondTerm))
	{
		Log(logERROR) << "Error: T non-finite (" << secondTerm << ")";
		Log(logERROR) << "Note ep = " << gsl_vector_get(xi,xiCounter+2) << " at " << i;
		exit(1);
	}
	if (firstTerm >= secondTerm)
			return firstTerm;
	else
			return secondTerm;
}


double ComputeL(gsl_vector * xi,constants * modelConst,int i)
{
	double firstTerm,secondTerm; //see doc.  
	double xiCounter = 5*(i-1); //counter relative to xi.  

	Log(logDEBUG1) << "Computing L";
	firstTerm = pow(gsl_vector_get(xi,xiCounter+1),1.5)/gsl_vector_get(xi,xiCounter+2);
	if (!isfinite(firstTerm))
	{
		Log(logERROR) << "Error: L non-finite (" << firstTerm << ")";
		exit(1);
	}
		
		
	secondTerm = modelConst->Ceta*pow(1/(pow(modelConst->reyn,3)*gsl_vector_get(xi,xiCounter+2)),0.25);
	if (!isfinite(secondTerm))
	{
		Log(logERROR) << "Error: L non-finite (" << secondTerm << ")";
		exit(1);
	}

	if (firstTerm >= secondTerm)
		return modelConst->CL*firstTerm;
	else
		return modelConst->CL*secondTerm;

}

double ComputeEddyVisc(gsl_vector * xi, gsl_vector * T, constants * modelConst,int i)
{
	double val; 
	double xiCounter = 5*(i-1); //counter relative to xi. -1 since U starts a 0. 

	Log(logDEBUG1) << "Computing Eddy Viscosity";
	val = modelConst->Cmu*gsl_vector_get(xi,xiCounter+3)*gsl_vector_get(T,i);
	if (!isfinite(val))
	{
		Log(logERROR) << "Error: vT non-finite (" << val << ")";
		exit(1);
	}
	return val; 
}

double ComputeP(gsl_vector * xi,gsl_vector* vT,double deltaEta,int i)
{
	double val; 

	Log(logDEBUG1) << "Computing P";

	//note: Diff1 takes xicounter indices
	val = gsl_vector_get(vT,i)*pow(Diff1(xi,deltaEta,0,5*(i-1)),2);
	if (!isfinite(val))
	{
		Log(logERROR) << "Error: P non-finite (" << val << ")";
		exit(1);
	}
	return val; 
}

double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) 
{
	Log(logDEBUG1) << "Compute dissipation at wall boundary";
	double ep0 = ((2*gsl_vector_get(xi,1))/(modelConst->reyn*pow(deltaEta,2)));
	if (!isfinite(ep0) || ep0 < 0)
	{
		Log(logERROR) << "Error: unacceptable ep0 (" << ep0 << ")";
		exit(1); 
	}
	return ep0; 
}

double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta)
{
	Log(logDEBUG1)<<"Compute f at wall boundary";
	double f0  = -(( (20*gsl_vector_get(xi,3))/( pow(modelConst->reyn,3)*ComputeEp0(xi,modelConst,deltaEta)*pow(deltaEta,4))));
	if(!isfinite(f0))
	{
		Log(logERROR) << "Error: f0 non-finite (" << f0 << ")";
		exit(1);
	}
	return f0; 
}
