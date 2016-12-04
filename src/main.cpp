//--------------------------------------------------
// main.cpp Main program for running v2f code. 
//
// 12/3/2016 - (gyalla) Written for CSE380 final project. 
//--------------------------------------------------
#include<cstdlib>
#include<iostream>
#include"setup.h"
#include<grvy.h>
#include<string>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<fstream>
#include"computeTerms.h"
#include"systemSolve.h"
using namespace GRVY;
using namespace std; 

//function declarations. 
int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaX);
int print_state(size_t iter,gsl_multiroot_fsolver * s);

int main(int argc, char ** argv)
{
	// Parse inputs 
	log(verbose,0,"Parsing inputs.\n"); 
	double deltaEta;
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename; 
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	if(Grvy_Input_Parse(modelConst,filename,deltaEta))
	{
		cerr << "Error parsing inputs" << endl; 
		return 1; 
	}
	gt.EndTimer("Getting Inputs");

	// Solving for initial conditions 
	log(verbose,0,"Solving initial conditions for U,k,ep,v2\n");
	gt.BeginTimer("Solving Initial Conditions");
	double I = 1/deltaEta; 
	gsl_vector * xi = gsl_vector_calloc(5*(I));
	if(SolveIC(xi,deltaEta,filename))
	{
		cerr << "Error interpolating initial conditions."<<endl; 
		return 1; 
	}
	log(verbose,0,"Solving initial conditions for f\n");

	if(Solve4f0(xi,modelConst,deltaEta))
	{
		cerr << "Error initializing f" << endl; 
		return 1; 
	}
	gt.EndTimer("Solving Initial Conditions");

	// Newton Solve. 
	log(verbose,0,"Solving system\n"); 
	gt.BeginTimer("Newton Solve + Time Marching");
	NewtonSolve(xi,modelConst,deltaEta);
	gt.EndTimer("Newton Solve + Time Marching");

	// summarize timiing.
	gt.Summarize();

	return 0; 
}

int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaEta)
{
	//for(unsigned int i=0;i<xi->size;i++)
	//{
//		cout << gsl_vector_get(xi,i) << endl; 
//	}

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
	log(verbose,1," " + string(gsl_strerror(status)) +  "\n");
	gsl_multiroot_fsolver_free(s); 
	//for(unsigned int i=0;i<xi->size;i++)
	//{
	//	cout << gsl_vector_get(s->x,i) << endl; 
	//}
	return 0; 
}

int print_state(size_t iter,gsl_multiroot_fsolver * s)
{
	log(verbose,1," iter = " + num2st(iter) + " Uend = " + num2st(gsl_vector_get(s->x,s->x->size-5)) + "\n"); 
	return 0; 
}

