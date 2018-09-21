//--------------------------------------------------
// main.cpp Main program for running v2f code. 
//
// 12/3/2016 - (gyalla) Written for CSE380 final project. 
//--------------------------------------------------
#include<iostream>
#include<iomanip>
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

using namespace std; 
//function declarations. 
int NewtonSolve(gsl_vector * xi,constants * modelConst, Grid* grid, int max_ts);
int print_state(int i,gsl_multiroot_fsolver * s, string status, double deltaT, double maxres);
void Print_Program_Info();

std::string NumberToString ( int Number);
int main(int argc, char ** argv)
{
	// Parse inputs 
	Print_Program_Info();
	Log(logINFO) << "Parsing inputs";
	bool uniform_grid, restarting;
        int max_ts;
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename, outFile;
	if(Input_Parse(modelConst,filename,outFile, uniform_grid, max_ts,restarting,argc,argv))
	{
		Log(logERROR) << "Error parsing inputs";
		return 1; 
	}

	// Make a new grid object
	Grid grid(uniform_grid, 1.0, 1.0/Const.reyn);
	Log(logINFO) << "---> Number of grid points = " << grid.getSize();

	// Solving for initial conditions 
	Log(logINFO) << "Solving initial conditions for U,k,ep,v2";
	double I = grid.getSize();
	gsl_vector * xi = gsl_vector_calloc(5*(I));

	if(SolveIC(xi,modelConst,&grid,filename,restarting))
	{
		Log(logERROR) << "Error interpolating initial conditions.";
		return 1; 
	}

	if (!restarting)
	{
	        Log(logINFO) << "Solving initial conditions for f";
		if(Solve4f0(xi,modelConst,&grid))
		{
			Log(logERROR) << "Error initializing f";
			return 1; 
		}
	}
	

	SaveResults(xi,"../data/init.dat",&grid,modelConst);
	// Newton Solve. 
	Log(logINFO) << "Solving system...";
	NewtonSolve(xi,modelConst,&grid,max_ts);

	//writing data to output
	Log(logINFO) << "Writing results to " << outFile;
	SaveResults(xi,outFile,&grid,modelConst);

	return 0; 
}

int NewtonSolve(gsl_vector * xi,constants * modelConst, Grid* grid, int max_ts)
{
        double previous_residual=100; // The max residual at the previous step
        double max_residual = 100;    // The max residual at the current step
        double deltaT = 0.000001;
        int diverged_count = 0;       // Tracks the number of steps where
                                      // the residual is diverging.
        const double residual_switch = 0.5; // deltaT won't increase if residual is above this limit
        const double max_deltaT = 1000.0;      // maximum possible value of deltaT
	int status;  // status of solver
	int iter = 0; 
	int inner_iter = 0; 
	//set up solver
	Log(logINFO) <<"Setting up Solver";
	const gsl_multiroot_fsolver_type * Type = gsl_multiroot_fsolver_dnewton; //dnewton for newton solver;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(Type,xi->size);
	//for time marching, starting small and getting bigger works best. 
	//power here represents powers of 2 for deltaT.
	do
	{
		//deltaT = fmin(0.0001,pow(10,power));
		/*
                deltaT = 0.01;
		if (iter > 5000)
			deltaT = 0.1;
		if (iter > 7000)
			deltaT = 1;
		if (iter > 8000)
			deltaT = 10; 
		if (iter > 9000)
			deltaT = 100; 
		if (iter > 10000)
			deltaT = 1000; 
		if (iter > 10000)
			deltaT = 10;
			*/
	
                //deltaT = 1/modelConst->reyn + iter*pow(2,power); // start off 1/modelConst->reyn; 
		iter++;
		//deltaT = fmin(pow(10,3),pow(2,power)); // start off 1/modelConst->reyn; 
		struct FParams p = {xi,deltaT,grid,modelConst};
		FParams * params = &p; 
		gsl_multiroot_function F = {&SysF,xi->size,params};
		//Log(logINFO) << "Setting up System";
		gsl_multiroot_fsolver_set(s,&F,xi); 
		//only need one iteration per deltaT since we don't care about temporal accuracy. 
		//We are just trying to get to the fully developed region of flow. 
		do
		{
			inner_iter++; 
			status = gsl_multiroot_fsolver_iterate(s);
			if(status)
				break;
			//status = gsl_multiroot_test_residual(s->f,0.001);
			print_state(iter,s,string(gsl_strerror(status)),deltaT,max_residual); 
		} while(inner_iter < 1);
		inner_iter =0;
		//xi = gsl_multiroot_fsolver_root(s); 
		//power +=1; 

		for (unsigned int i = 0; i < xi->size; i++)
		{
			gsl_vector_set(xi,i,gsl_vector_get(s->x,i));
			if(i%5==1)
				gsl_vector_set(xi,i,fmax(gsl_vector_get(xi,i),K_MIN));
			if(i%5==3)
				gsl_vector_set(xi,i,fmax(gsl_vector_get(xi,i),V2_MIN));
		}

		if (iter%50 == 0)
			SaveResults(xi,"../data/test/solve" + NumberToString(iter)  + ".dat",grid,modelConst);
                
                // Change the time-step
                if (iter > 1) previous_residual = max_residual;
                max_residual = gsl_vector_max(s->f);
                if (max_residual/previous_residual > 2.0) {
                  diverged_count++;
                } else {
                  diverged_count = 0;
                }
                if (diverged_count > 3) deltaT /= 10;
                if (max_residual < residual_switch && deltaT < max_deltaT &&
                    max_residual/previous_residual > 0.2 &&
                    max_residual/previous_residual < 1.0) {
                  deltaT *= 2;
                }

		status = gsl_multiroot_test_residual (s->f, 1e-7);
	}while(status == GSL_CONTINUE && iter < max_ts);
	//gsl_multiroot_fsolver_free(s); 
	return 0; 
}

int print_state(int i,gsl_multiroot_fsolver * s,string status, double deltaT, double maxres)
{
	Log(logINFO) << setw(11)<< "Iteration: " << setw(7) << std::left <<  i << "\tdeltaT = " << setw(10) << std::left << setprecision(5) << deltaT 
		<< setw(14) << "\tMax Residual: " << setw(10) << std::left << setprecision(5) << maxres << "\t GSL SOLVER STATUS: " << status;
	//Log(logINFO) << "iteration = "<< i << ", max = " << gsl_vector_max(s->f);
	return 0; 
}

std::string NumberToString ( int Number)
{
	std::ostringstream ss;
        ss << Number;
        return ss.str();
}

void Print_Program_Info()
{
	cout<<"\n";
	cout<< "____   ____________          _______________ __________    \n"; 
	cout<< "\\   \\ /   /\\_____  \\         \\_   _____/    |   \\      \\   \n";
	cout<<" \\   Y   /  /  ____/   ______ |    __) |    |   /   |   \\  \n";
	cout<<"  \\     /  /       \\  /_____/ |     \\  |    |  /    |    \\ \n";
	cout<<"   \\___/   \\_______ \\         \\___  /  |______/\\____|__  / \n";
	cout<<"   		   \\/             \\/                   \\/  \n";
	cout<<"\n";
}
