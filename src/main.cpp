//--------------------------------------------------
// main.cpp Main program for running v2f code. 
//
// 12/3/2016 - (gyalla) Written for CSE380 final project. 
//--------------------------------------------------
#include<iostream>
#include<grvy.h>
#include<gsl/gsl_multiroots.h>
#include<math.h>
#include"systemSolve.h"
#include "Grid.h"
#include<string>
#include <sstream>

#define EP_MIN 1.0e-7
#define K_MIN  1.0e-7
#define V2_MIN 1.0e-12
#define T_MIN  1.0e-7
#define L_MIN  1.0e-5

using namespace GRVY;
using namespace std; 
//function declarations. 
int NewtonSolve(gsl_vector * xi,constants * modelConst, Grid* grid, int max_ts);
int print_state(int i,gsl_multiroot_fsolver * s);

std::string NumberToString ( int Number);
int main(int argc, char ** argv)
{
	// Parse inputs 
	Log(logINFO) << "Parsing inputs";
	bool   uniform_grid;
        int max_ts;
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename, outFile;
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	if(Grvy_Input_Parse(modelConst,filename,outFile, uniform_grid, max_ts))
	{
		Log(logERROR) << "Error parsing inputs";
		return 1; 
	}
	gt.EndTimer("Getting Inputs");

	// Make a new grid object
	Grid grid(uniform_grid, 1.0, 1.0/Const.reyn);
	Log(logINFO) << "---> Number of grid points = " << grid.getSize();

	// Solving for initial conditions 
	Log(logINFO) << "Solving initial conditions for U,k,ep,v2";
	gt.BeginTimer("Solving Initial Conditions");
	double I = grid.getSize();
	gsl_vector * xi = gsl_vector_calloc(5*(I));
	if(SolveIC(xi,modelConst,&grid,filename))
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

	SaveResults(xi,"../data/init.dat",&grid,modelConst);
	// Newton Solve. 
	Log(logINFO) << "Solving system...";
	gt.BeginTimer("Newton Solve + Time Marching");
	NewtonSolve(xi,modelConst,&grid,max_ts);
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

int NewtonSolve(gsl_vector * xi,constants * modelConst, Grid* grid, int max_ts)
{
	double deltaT;
	int status;  // status of solver
	int power = -10;
	int iter = 0; 
	int inner_iter = 0; 
	bool converge = false; 
	//set up solver
	Log(logINFO) <<"Setting up Solver";
	const gsl_multiroot_fsolver_type * Type = gsl_multiroot_fsolver_dnewton; //dnewton for newton solver;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(Type,xi->size);
	//for time marching, starting small and getting bigger works best. 
	//power here represents powers of 2 for deltaT.
	do
	{
		//deltaT = fmin(0.0001,pow(10,power));
		deltaT = 0.00001;
		if (iter > 100)
			deltaT = 0.0001;
		if (iter > 300)
			deltaT = 0.001;

		//deltaT = 1/modelConst->reyn + iter*pow(2,power); // start off 1/modelConst->reyn; 
		iter++;
		//deltaT = fmin(pow(10,3),pow(2,power)); // start off 1/modelConst->reyn; 
		struct FParams p = {xi,deltaT,grid,modelConst};
		FParams * params = &p; 
		gsl_multiroot_function F = {&SysF,xi->size,params};
		Log(logINFO) << "Setting up System";
		gsl_multiroot_fsolver_set(s,&F,xi); 
		//only need one iteration per deltaT since we don't care about temporal accuracy. 
		//We are just trying to get to the fully developed region of flow. 
		do
		{
			inner_iter++; 
			Log(logINFO) << "Iterating (deltaT = " << deltaT << ")";
			status = gsl_multiroot_fsolver_iterate(s);
			print_state(iter,s); 
			if(status)
				break;
			//status = gsl_multiroot_test_residual(s->f,0.001);
			Log(logINFO) << string(gsl_strerror(status));
		} while(inner_iter < 1);
		inner_iter =0;
		//xi = gsl_multiroot_fsolver_root(s); 
		//power +=1; 

		for (unsigned int i = 0; i < xi->size; i++)
		{
			gsl_vector_set(xi,i,gsl_vector_get(s->x,i));
			if(i%5==1)
				gsl_vector_set(xi,i,fmax(gsl_vector_get(xi,i),K_MIN));
			//if(i%5==2)
			//	gsl_vector_set(xi,i,fmax(gsl_vector_get(xi,i),EP_MIN));
			if(i%5==3)
				gsl_vector_set(xi,i,fmax(gsl_vector_get(xi,i),V2_MIN));
		}

		SaveResults(xi,"../data/solve.dat",grid,modelConst);

		if (iter%200 == 0)
			SaveResults(xi,"../data/test/solve" + NumberToString(iter)  + "_step001.dat",grid,modelConst);

		for(unsigned int i=0; i<xi->size;i++)
		{
			if((fabs(gsl_vector_get(xi,i)-gsl_vector_get(params->XiN,i))<0.000001) && iter > max_ts)
				converge = true;
		}
	}while(!converge);
	//gsl_multiroot_fsolver_free(s); 
	return 0; 
}

int print_state(int i,gsl_multiroot_fsolver * s)
{
	Log(logINFO) << "iteration = "<< i << ", U6 = " << gsl_vector_get(s->x,25);
	return 0; 
}

std::string NumberToString ( int Number)
{
	std::ostringstream ss;
        ss << Number;
        return ss.str();
}
