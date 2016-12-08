//--------------------------------------------------
// main.cpp Main program for running v2f code. 
//
// 12/3/2016 - (gyalla) Written for CSE380 final project. 
//--------------------------------------------------
#include<cstdlib>
#include<iostream>
#include<iomanip>
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
#include"loglevel.h"
using namespace GRVY;
using namespace std; 

//function declarations. 
int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaX);
int print_state(int i,gsl_multiroot_fsolver * s);

int main(int argc, char ** argv)
{
	// Parse inputs 
	Log(logINFO) << "Parsing inputs";
	double deltaEta;
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename, outFile;
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	if(Grvy_Input_Parse(modelConst,filename,outFile,deltaEta))
	{
		Log(logERROR) << "Error parsing inputs";
		return 1; 
	}
	gt.EndTimer("Getting Inputs");

	// Solving for initial conditions 
	Log(logINFO) << "Solving initial conditions for U,k,ep,v2";
	gt.BeginTimer("Solving Initial Conditions");
	double I = 1/deltaEta; 
	gsl_vector * xi = gsl_vector_calloc(5*(I));
	if(SolveIC(xi,deltaEta,filename))
	{
		Log(logERROR) << "Error interpolating initial conditions.";
		return 1; 
	}
	Log(logINFO) << "Solving initial conditions for f";

	if(Solve4f0(xi,modelConst,deltaEta))
	{
		Log(logERROR) << "Error initializing f";
		return 1; 
	}
	gt.EndTimer("Solving Initial Conditions");

	// Newton Solve. 
	Log(logINFO) << "Solving system";
	gt.BeginTimer("Newton Solve + Time Marching");
	NewtonSolve(xi,modelConst,deltaEta);
	gt.EndTimer("Newton Solve + Time Marching");

	//writing data to output
	Log(logINFO) << "Writing results to output file";
	gt.BeginTimer("Writing results to output"); 
	SaveResults(xi,outFile,deltaEta,modelConst); 
	gt.EndTimer("Writing results to output"); 

	// summarize timiing.
	gt.Summarize();

	return 0; 
}

int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaEta)
{
	double deltaT;
	int status;  // status of solver
	int power = -30;
	int iter = 0; 
	bool converge = false; 
	//set up solver
	Log(logINFO) <<"Setting up Solver...";
	const gsl_multiroot_fsolver_type * Type = gsl_multiroot_fsolver_hybrids; //dnewton for newton solver;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(Type,xi->size);
	//for time marching, starting small and getting bigger works best. 
	//power here represents powers of 2 for deltaT.  
	do
	{
		iter++;
		deltaT = pow(2,power);
		//deltaT = 1/modelConst->reyn + pow(2,power); // start off 1/modelConst->reyn; 
		struct FParams p = {xi,deltaT,deltaEta,modelConst}; 
		FParams * params = &p; 
		Log(logINFO) << "Setting up System...";
		gsl_multiroot_function F = {&SysF,xi->size,params};
		gsl_multiroot_fsolver_set(s,&F,xi); 
		
		//only need one interation per deltaT since we don't care about temporal accuracy. 
		//We are just trying to get to the fully developed region of flow. 
		Log(logINFO) << "Iterating...";
		status = gsl_multiroot_fsolver_iterate(s);
		print_state(iter,s); 
		if(status)
			break;
		//status = gsl_multiroot_test_residual(s->f,0.001);
		Log(logINFO) << string(gsl_strerror(status));
		xi = s->x; 
		power +=5; 
		
		//check if we are in fully developed region
		for(unsigned int i=0; i<xi->size;i++)
		{
			if(fabs(gsl_vector_get(xi,i)-gsl_vector_get(params->XiN,i))<0.001)
				converge = true;
		}
	}while(!converge);

	gsl_multiroot_fsolver_free(s); 
	return 0; 
}

int print_state(int i,gsl_multiroot_fsolver * s)
{
	Log(logINFO) << "iteration = "<< i << ", Uend = " << gsl_vector_get(s->x,s->x->size-5);
	return 0; 
}

