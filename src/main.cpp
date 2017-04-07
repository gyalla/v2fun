//--------------------------------------------------
// main.cpp Main program for running v2f code. 
//
// 12/3/2016 - (gyalla) Written for CSE380 final project. 
//--------------------------------------------------
#include<iostream>
#include<grvy.h>
#include<gsl/gsl_multiroots.h>
#include"systemSolve.h"
#include "Grid.h"
using namespace GRVY;
using namespace std; 

//function declarations. 
int NewtonSolve(gsl_vector * xi,constants * modelConst, Grid* grid);
int print_state(int i,gsl_multiroot_fsolver * s);

int main(int argc, char ** argv)
{
	// Parse inputs 
	Log(logINFO) << "Parsing inputs";
	bool   uniform_grid;
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename, outFile;
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	if(Grvy_Input_Parse(modelConst,filename,outFile, uniform_grid))
	{
		Log(logERROR) << "Error parsing inputs";
		return 1; 
	}
	gt.EndTimer("Getting Inputs");

	// Make a new grid object
	Grid grid(uniform_grid, 1.0, 1.0/Const.reyn);
  double deltaEta = 1.0/Const.reyn; // XXX: Fix this for nonuniform flow.

	// Solving for initial conditions 
	Log(logINFO) << "Solving initial conditions for U,k,ep,v2";
	gt.BeginTimer("Solving Initial Conditions");
	double I = grid.getSize();
	gsl_vector * xi = gsl_vector_calloc(5*(I));
	if(SolveIC(xi,&grid,filename))
	{
		Log(logERROR) << "Error interpolating initial conditions.";
		return 1; 
	}
	Log(logINFO) << "Solving initial conditions for f";

	if(Solve4f0(xi,modelConst,&grid))
	{
		Log(logERROR) << "Error initializing f";
		return 1; 
	}
	gt.EndTimer("Solving Initial Conditions");

	// Newton Solve. 
	Log(logINFO) << "Solving system...";
	gt.BeginTimer("Newton Solve + Time Marching");
	NewtonSolve(xi,modelConst,&grid);
	gt.EndTimer("Newton Solve + Time Marching");

	//writing data to output
	Log(logINFO) << "Writing results to " << outFile;
	gt.BeginTimer("Writing results to output"); 
	SaveResults(xi,outFile,&grid,modelConst);
	gt.EndTimer("Writing results to output"); 

	// summarize timiing.
	gt.Summarize();

	return 0; 
}

int NewtonSolve(gsl_vector * xi,constants * modelConst, Grid* grid)
{
	double deltaT;
	int status;  // status of solver
	int power = -20;
	int iter = 0; 
	bool converge = false; 
	//set up solver
	Log(logINFO) <<"Setting up Solver";
	const gsl_multiroot_fsolver_type * Type = gsl_multiroot_fsolver_hybrids; //dnewton for newton solver;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(Type,xi->size);
	//for time marching, starting small and getting bigger works best. 
	//power here represents powers of 2 for deltaT.  
	double deltaEta = gsl_vector_get(grid->y, 0);
	do
	{
		//The first deltaT will run. The second is actual time marching. 
		//Comment out the next line and uncomment the following for time marching. 
		//This will lead to errors in the dissipation term however. 
		//deltaT = pow(2,-30);
		//deltaT = 1/modelConst->reyn + iter*pow(2,power); // start off 1/modelConst->reyn; 
		deltaT = pow(2,power); // start off 1/modelConst->reyn; 
		iter++;
		struct FParams p = {xi,deltaT,grid,modelConst};
		FParams * params = &p; 
		Log(logINFO) << "Setting up System";
		gsl_multiroot_function F = {&SysF,xi->size,params};
		gsl_multiroot_fsolver_set(s,&F,xi); 
		
		//only need one interation per deltaT since we don't care about temporal accuracy. 
		//We are just trying to get to the fully developed region of flow. 
		Log(logINFO) << "Iterating (deltaT = " << deltaT << ")";
		status = gsl_multiroot_fsolver_iterate(s);
		print_state(iter,s); 
		if(status)
			break;
		//status = gsl_multiroot_test_residual(s->f,0.001);
		Log(logINFO) << string(gsl_strerror(status));
		xi = s->x; 
		power +=2; 
		
		//check if we are in fully developed region
		for(unsigned int i=0; i<xi->size;i++)
		{
			if((fabs(gsl_vector_get(xi,i)-gsl_vector_get(params->XiN,i))<0.000001) && (iter > 10))
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

