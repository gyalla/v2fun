#include<cstdlib>
#include<iostream>
#include"setup.h"
#include<grvy.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<fstream>
#include"computeTerms.h"
#include"systemSolve.h"
using namespace GRVY;
using namespace std; 

int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaX);

int main(int argc, char ** argv)
{
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename; 
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	Grvy_Input_Parse(modelConst,filename); 
	gt.EndTimer("Getting Inputs");

	cout << filename << endl;

	double deltaEta = 1/modelConst->reyn;
	gsl_vector * U = gsl_vector_calloc(1/deltaEta + 1); //mean velocity
	gsl_vector * k = gsl_vector_calloc(1/deltaEta + 1); //turbulent kinetic energy
	gsl_vector * ep = gsl_vector_calloc(1/deltaEta + 1); //energy dissipation
	gsl_vector * v2 = gsl_vector_calloc(1/deltaEta + 1); //turbulent velocity scale
	gsl_vector * f  = gsl_vector_calloc(1/deltaEta + 1); // redistribution term 
	gsl_vector * P  = gsl_vector_calloc(1/deltaEta + 1); // production rate
	gsl_vector * T  = gsl_vector_calloc(1/deltaEta + 1); // turbulent time scale
	gsl_vector * L  = gsl_vector_calloc(1/deltaEta + 1); // turbulent length scale 
	gsl_vector * vT = gsl_vector_calloc(1/deltaEta + 1); // eddy viscosity 
	gsl_vector * xi = gsl_vector_calloc(5*(U->size-1));

	SolveIC(U,k,ep,v2,deltaEta,filename);	

	ComputeT(k,ep,modelConst,T); 
	ComputeEddyVisc(v2,T,modelConst,vT); 
	ComputeP(U,vT,deltaEta,P);
	ComputeL(k,ep,modelConst,L); 
	Solve4f0(k,ep,v2,P,T,L,modelConst,deltaEta,f); 

	ReconstructXi(xi,U,k,ep,v2,f); 	


	NewtonSolve(xi,modelConst,deltaEta);

	
	gsl_vector_free(P);
	gsl_vector_free(T);
	gsl_vector_free(L);
	gsl_vector_free(vT);
	gsl_vector_free(U);
	gsl_vector_free(ep);
	gsl_vector_free(v2);
	gsl_vector_free(f); 

}

int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaEta)
{
	gsl_vector *xin;
	gsl_vector_memcpy(xin,xi); 

	double deltaT = 1; 

	struct FParams p = {xin,deltaT,deltaEta,modelConst}; 
	FParams * params = &p; 
	gsl_multiroot_function FDF = {&SysF,xi->size,params};

	const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_dnewton;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(T,xi->size);
	gsl_multiroot_fsolver_set(s,&FDF,xi); 

	int iter = 0; 
	int status; 
	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);
		if(status)
			break;
		status = gsl_multiroot_test_residual(s->f,1);
	}
	while(status==GSL_CONTINUE && iter < 1);  
	cout << gsl_strerror(status) << endl; 
	gsl_multiroot_fsolver_free(s); 


	return 0; 
}
