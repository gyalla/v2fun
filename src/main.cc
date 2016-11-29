#include<cstdlib>
#include<iostream>
#include"options.h"
#include<grvy.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<fstream>
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

	SolveIC(U,k,ep,v2,deltaEta,filename);	
	for (int i =0;i<U->size;i++)
	{
		cout << gsl_vector_get(U,i) << endl;
	}
	return 0; 


}

