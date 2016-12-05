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
using namespace GRVY;
using namespace std; 

//function declarations. 
int NewtonSolve(gsl_vector * xi,constants * modelConst, double deltaX);
int print_state(int i,gsl_multiroot_fsolver * s);

int main(int argc, char ** argv)
{
	// Parse inputs 
	log(verbose,0,"Parsing inputs.\n"); 
	double deltaEta;
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	string filename, outFile;
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	if(Grvy_Input_Parse(modelConst,filename,outFile,deltaEta))
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

	//writing data to output
	log(verbose,0,"Writing results to output file");
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
	int power = -2;
	int iter ; 
	if (verbose >=2)
	{
		for (int i=0; i<xi->size;i++)
			{
				cout << gsl_vector_get(xi,i) << endl; 
			}
	}
	//set up solver
	const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid; //dnewton for newton solver;
	gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc(T,xi->size);

	//for time marching, starting small and getting bigger works best. 
	//power here represents powers of 2 for deltaT.  
	do
	{
		iter = 0; 
		deltaT = pow(2,power); 
		struct FParams p = {xi,deltaT,deltaEta,modelConst}; 
		FParams * params = &p; 
		gsl_multiroot_function F = {&SysF,xi->size,params};
		gsl_multiroot_fsolver_set(s,&F,xi); 
		do
		{
			iter ++;
		//only need one interation per deltaT since we don't care about temporal accuracy. 
		//We are just trying to get to the fully developed region of flow. 
		status = gsl_multiroot_fsolver_iterate(s);
		print_state(power,s); 
		if(status)
			break;
		status = gsl_multiroot_test_residual(s->f,0.001);
		if (verbose >=2)
		{
			for (int i=0; i<s->x->size;i++)
			{
				cout << gsl_vector_get(s->x,i) << endl; 
			}
		}
		}while(status == GSL_CONTINUE && iter < 1);
		log(verbose,2," " + string(gsl_strerror(status)) +  "\n");
		xi = s->x; 
		power +=2; 
	}
	while(power <=4);

	gsl_multiroot_fsolver_free(s); 
	return 0; 
}

int print_state(int i,gsl_multiroot_fsolver * s)
{
	double ans = gsl_vector_get(s->x,s->x->size-5); 
	log(verbose,1," deltaT = 2^" + num2st(i) + ", Uend = " + num2st(gsl_vector_get(s->x,s->x->size-5)) + "\n"); 
	cout << setprecision(15) << "deltaT = 2^" << i << " Uend = " << ans << endl; 
	return 0; 
}

