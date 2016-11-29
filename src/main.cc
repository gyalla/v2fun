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
using namespace GRVY;
using namespace std; 


int main(int argc, char ** argv)
{
	double reyn, Cmu,C1,C2,sigmaEp,CL,Cep1,Cep2,Ceta;
	string filename; 
	GRVY_Timer_Class gt; 
	gt.BeginTimer("Getting Inputs");
	Grvy_Input_Parse(reyn,Cmu,C1,C2,sigmaEp,CL,Cep2,Cep1,Ceta,filename); 
	gt.EndTimer("Getting Inputs");

	double deltaEta = 1/reyn;
		gsl_vector * U = gsl_vector_calloc(1/deltaEta + 1); //mean velocity
		gsl_vector * k = gsl_vector_calloc(1/deltaEta + 1); //turbulent kinetic energy
		gsl_vector * ep = gsl_vector_calloc(1/deltaEta + 1); //energy dissipation
		gsl_vector * v2 = gsl_vector_calloc(1/deltaEta + 1); //turbulent velocity scale
		gsl_vector * f  = gsl_vector_calloc(1/deltaEta + 1); // redistribution term 
		gsl_vector * P  = gsl_vector_calloc(1/deltaEta + 1); 
		gsl_vector * T  = gsl_vector_calloc(1/deltaEta + 1); 
		gsl_vector * L  = gsl_vector_calloc(1/deltaEta + 1); 
		gsl_vector * vT = gsl_vector_calloc(1/deltaEta + 1); 

	SolveIC(U,k,ep,v2,deltaEta,filename);	

	ComputeT(k,ep,reyn,T); 
	ComputeEddyVisc(v2,T,Cmu,vT); 
	ComputeP(U,vT,deltaEta,P);
	ComputeL(k,ep,reyn,CL,Ceta,L); 
	Solve4f0(k,ep,v2,P,T,L,reyn,C2,C1,deltaEta,f); 

}

