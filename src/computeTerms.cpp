//--------------------------------------------------
// computeTerms: Compute terms T,L,P,vT and f(0),ep(0).
// 
// 12/3/2016 - (gry88) Writen for final project CSE380.  
//-------------------------------------------------- 
#include<gsl/gsl_vector.h>
#include<math.h>
#include"setup.h"
#include"finiteDiff.h"
using namespace std;

#define EP_MIN 1.0e-7
#define K_MIN  1.0e-7
#define V2_MIN 1.0e-12
#define T_MIN  1.0e-7
#define L_MIN  1.0e-5
#define F_MIN  1.0e-8

double ComputeT(gsl_vector * xi, constants * modelConst,int i)
{
	double firstTerm,secondTerm; //1st and 2nd term as in documentation. 
	double xiCounter = 5*(i-1);  //counter relative to xi. 	
	
	Log(logDEBUG1) << "Computing T";	

	double k = fmax(gsl_vector_get(xi,xiCounter+1),K_MIN);
	double ep = fmax(gsl_vector_get(xi,xiCounter+2),EP_MIN);

	firstTerm = k/ep;
	if (!isfinite(firstTerm))
	{
		Log(logERROR) << "Error: T non-finite (" << firstTerm << ")";
		Log(logERROR) << "-Note ep = " << gsl_vector_get(xi,xiCounter+2) << " at " << i;
		exit(1);
	}

	secondTerm = 6*sqrt(1/(modelConst->reyn*ep));
	if(!isfinite(secondTerm))
	{
		Log(logERROR) << "Error: T non-finite (" << secondTerm << ")";
		Log(logERROR) << "Note ep = " << gsl_vector_get(xi,xiCounter+2) << " at " << i;
		exit(1);
	}

	return fmax(fmax(firstTerm,secondTerm),T_MIN);
}

double ComputeL(gsl_vector * xi,constants * modelConst,int i)
{
	double firstTerm,secondTerm; //see doc.  
	double xiCounter = 5*(i-1); //counter relative to xi.  

	double k = fmax(gsl_vector_get(xi,xiCounter+1),K_MIN);
	double ep = fmax(gsl_vector_get(xi,xiCounter+2),EP_MIN);

	Log(logDEBUG1) << "Computing L";
	firstTerm = pow(k,1.5)/ep;
	if (!isfinite(firstTerm))
	{
		Log(logERROR) << "Error: L non-finite (" << firstTerm << ")";
		exit(1);
	}
		
		
	secondTerm = modelConst->Ceta*pow(1/(pow(modelConst->reyn,3)*ep),0.25);
	if (!isfinite(secondTerm))
	{
		Log(logERROR) << "Error: L non-finite (" << secondTerm << ")";
		exit(1);
	}

	return fmax(modelConst->CL*fmax(firstTerm,secondTerm),L_MIN);
}

double ComputeEddyVisc(gsl_vector * xi, gsl_vector * T, constants * modelConst,int i)
{
	double val; 
	double xiCounter = 5*(i-1); //counter relative to xi. -1 since U starts a 0. 

	double v2 = fmax(gsl_vector_get(xi,xiCounter+3),V2_MIN);

	Log(logDEBUG1) << "Computing Eddy Viscosity";
	val = modelConst->Cmu*v2*gsl_vector_get(T,i);
	if (!isfinite(val))
	{
		Log(logERROR) << "Error: vT non-finite (" << val << ")";
		exit(1);
	}
	return val; 
}

double ComputeP(gsl_vector * xi, gsl_vector* vT, Grid* grid, int i)
{
	double val; 

	Log(logDEBUG1) << "Computing P";

	//note: Diff1 takes xi-counter indices
	val = gsl_vector_get(vT,i)*pow(Deriv1(xi,0.0,5*(i-1),grid),2);
	if (!isfinite(val))
	{
		Log(logERROR) << "Error: P non-finite (" << val << ")";
		exit(1);
	}
	return val; 
}

double ComputeEp0(gsl_vector * xi,constants * modelConst, Grid* grid)
{
	Log(logDEBUG1) << "Compute dissipation at wall boundary";
	double delta_y_0 = gsl_vector_get(grid->y, 0);
	double ep0 = ((2*gsl_vector_get(xi,1))/(modelConst->reyn*pow(delta_y_0,2)));
	if (!isfinite(ep0)) //|| ep0 < 0)
	{
		Log(logERROR) << "Error: unacceptable ep0 (" << ep0 << ")";
		exit(1); 
	}
	return ep0; 
}

double Computef0(gsl_vector * xi,constants * modelConst, Grid* grid)
{
	Log(logDEBUG1)<<"Compute f at wall boundary";
  double delta_y_0 = gsl_vector_get(grid->y, 0);
	double f0  = -(( (20*gsl_vector_get(xi,3))/( pow(modelConst->reyn,3) *
	    ComputeEp0(xi, modelConst, grid) * pow(delta_y_0, 4))));
	if(!isfinite(f0))
	{
		Log(logERROR) << "Error: f0 non-finite (" << f0 << ")";
		exit(1);
	}
	return f0; 
}
