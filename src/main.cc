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
int print_state(size_t iter,gsl_multiroot_fsolver * s);

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
	double deltaEta = 0.01;
	double I = 1/deltaEta; 

	gsl_vector * xi = gsl_vector_calloc(5*(I));

	SolveIC(xi,deltaEta,filename);	
	Solve4f0(xi,modelConst,deltaEta); 

	for(int i=0;i<xi->size;i++)
	{
		cout << gsl_vector_get(xi,i) << endl; 
	}

	NewtonSolve(xi,modelConst,deltaEta);
	return 0; 

}

int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaEta)
{
	gsl_vector *xin = gsl_vector_calloc(xi->size);
	gsl_vector_memcpy(xin,xi); 

	double deltaT = 100; 

	struct FParams p = {xin,deltaT,deltaEta,modelConst}; 
	FParams * params = &p; 
	gsl_multiroot_function F = {&SysF,xi->size,params};

	const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid; //dnewton;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(T,xi->size);
	gsl_multiroot_fsolver_set(s,&F,xi); 

	int iter = 0; 
	int status; 
	do
	{

		iter++;
		status = gsl_multiroot_fsolver_iterate(s);
		print_state(iter,s); 
		if(status)
			break;
		status = gsl_multiroot_test_residual(s->f,0.001);
	}
	while(status==GSL_CONTINUE && iter < 1);  
	cout << gsl_strerror(status) << endl; 
	gsl_multiroot_fsolver_free(s); 
	for(int i=0;i<xi->size;i++)
	{
		cout << gsl_vector_get(s->x,i) << endl; 
	}
	return 0; 
}

int print_state(size_t iter,gsl_multiroot_fsolver * s)
{
	cout << "iter = " << iter << " Uend = " << gsl_vector_get(s->x,s->x->size-5); 
	return 0; 
}

