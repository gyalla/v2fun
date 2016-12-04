//--------------------------------------------------
// setup: Initializes terms for v2-f code. 
//
// 12/3/2016 (gry88) - Written for CSE380 Final Project. 
//--------------------------------------------------
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<sstream>
#include<grvy.h>
#include "setup.h"
#include "computeTerms.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<fstream>
#include<math.h>
using namespace std;
using namespace GRVY;

int verbose; 

int Grvy_Input_Parse(constants * modelConst,string & filename,double & deltaEta)
{
	// Use GRVY for inpute parsing as shown in grvy documentation. 
	
	GRVY_Input_Class iparse; // parsing object

	if(!iparse.Open("./input_file.txt"))
		exit(1);

	// read in each variable. 
	
	if(!iparse.Read_Var("verbose",&verbose))
		return 1; 
	if(!iparse.Read_Var("reyn", &(modelConst->reyn)))
		return 1; 
	log(verbose,1," ---> reyn = " + num2st(modelConst->reyn) + "\n");

	if(!iparse.Read_Var("Cmu", &(modelConst->Cmu)))
		return 1; 
	log(verbose,1," ---> Cmu = " + num2st(modelConst->Cmu) + "\n");

	if(!iparse.Read_Var("C1", &(modelConst->C1)))
		return 1; 
	log(verbose,1," ---> C1 = " + num2st(modelConst->C1) + "\n");

	if(!iparse.Read_Var("C2", &(modelConst->C2)))
		return 1; 
	log(verbose,1," ---> C2 = " + num2st(modelConst->C2) + "\n");

	if(!iparse.Read_Var("sigmaEp", &(modelConst->sigmaEp)))
		return 1; 
	log(verbose,1," ---> sigmaEp = " + num2st(modelConst->sigmaEp) + "\n");

	if(!iparse.Read_Var("CL", &(modelConst->CL)))
		return 1; 
	log(verbose,1," ---> CL = " + num2st(modelConst->CL) + "\n");

	if(!iparse.Read_Var("Cep1",&(modelConst->Cep1)))
		return 1; 
	log(verbose,1," ---> Cep1 = " + num2st(modelConst->Cep1) + "\n");

	if(!iparse.Read_Var("Cep2", &(modelConst->Cep2)))
		return 1; 
	log(verbose,1," ---> Cep2 = " + num2st(modelConst->Cep2) + "\n");

	if(!iparse.Read_Var("Ceta", &(modelConst->Ceta)))
		return 1; 
	log(verbose,1," ---> Ceta = " + num2st(modelConst->Ceta) + "\n");
	
	if(!iparse.Read_Var("filename",&filename))
		return 1; 
	log(verbose,1," ---> filename = " + filename + "\n");

	if(!iparse.Read_Var("deltaEta",&deltaEta))
		return 1; 
	log(verbose,1," ---> deltaEta = " + num2st(deltaEta) + "\n");

	iparse.Close();
	
	return 0;
}

void log(int verbose, int level, string message)
{
	// log level structure for output. 
	if (verbose  >= level)
		cout << message; 
}

string num2st(double number)
{
	string result; 
	ostringstream convert; 
	convert << number; 
	result = convert.str();
	return result; 
}
		
int LinInterp(gsl_vector * Vec,double pt1,double pt2, double U1,double U2,double gridPt,int i )
{
	double val = U1 + (gridPt-pt1)*( (U2-U1)/(pt2-pt1));
	if (!isfinite(val))
	{
		cerr << "Error: non-finite interpolation" << endl; 
		return 1; 
	}
	gsl_vector_set(Vec,i,val); 
	return 0;
}

int SolveIC(gsl_vector * xi, double deltaEta, string file)
{	
	ifstream inFile; 
	inFile.open(file.c_str()); 
	
	// terms for linear interpolation as defined in doc. 
	double pt1,tempU1,tempK1,tempEp1,tempV21;
	double pt2,tempU2,tempK2,tempEp2,tempV22;
	double gridPt,xiCounter;
		
	// read in first two files. 
	if (!(inFile >> pt1 >> tempU1 >> tempK1 >> tempEp1 >> tempV21))
		return 1; 

	if (!(inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
		return 1; 
	for(unsigned int i=0; i<(xi->size/float(5));i++)
	{
		xiCounter=5*i; //counter relative to xi vector. 
		gridPt=(i+1)*deltaEta; //which grid point we are working on. 
		log(verbose,2,"  Working on grid point " + num2st(gridPt) + "\n");
		//find two points in Moser data that surrounds our grid points. 
		while ((! ((pt1 <= gridPt) && (pt2 >= gridPt))) && (!inFile.eof()))
		{
			//transfer variables and read in next line of data. 
			pt1 = pt2; 
			tempU1 = tempU2;
			tempK1 = tempK2;
			tempEp1 = tempEp2; 
			tempV21 = tempV22; 
			if (! (inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
				cout << "error" <<endl; 
		}

		//interpolate each term U,k,ep,v2,f.
		log(verbose,3,"   interpolation points " + num2st(pt1) + "," + num2st(pt2) + "\n");
		log(verbose,3,"   interpolation values (U): " + num2st(tempU1) + "," + num2st(tempU2) + "\n");
		log(verbose,3,"   interpolation values (k): " + num2st(tempK1) + "," + num2st(tempK2) + "\n");
		log(verbose,3,"   interpolation values (ep): " + num2st(tempEp1) + "," + num2st(tempEp1) + "\n");
		log(verbose,3,"   interpolation values (v2): " + num2st(tempV21) + "," + num2st(tempV21) + "\n");
		if(LinInterp(xi,pt1,pt2,tempU1,tempU2,gridPt,xiCounter))
			return 1; 
		if(LinInterp(xi,pt1,pt2,tempK1,tempK2,gridPt,xiCounter+1))
			return 1; 
		if(LinInterp(xi,pt1,pt2,tempEp1,tempEp2,gridPt,xiCounter+2))
			return 1; 
		if(LinInterp(xi,pt1,pt2,tempV21,tempV22,gridPt,xiCounter+3))
			return 1; 
	}

	inFile.close();

	return 0; 
}

int Solve4f0(gsl_vector * xi, constants * modelConst,double deltaEta)
{
	//Solve for f_0 based on data from others terms. This function needs to be tested still
	//but I am not sure what correct values of f_0 should be!  

	int size = xi->size/float(5) + 1;  // size for f. 
	gsl_matrix * A = gsl_matrix_calloc(size,size);
	gsl_vector * b = gsl_vector_calloc(A->size1); 
	gsl_vector * f = gsl_vector_calloc(A->size1);

	double LOvrEta; //main terms in matrix. 
	double LHS1,LHS2; // LHS of f from finite difference. 
	unsigned int i,xiCounter; 

	// Set matrix so solve
	//
	//      |    1           0           0            |
	// A = 	| L^2/n^2 -(2*L^2/n^2 + 1) L^2/n^2        |
	// 	|    0         2*L^2/n^2 -(2*L^2/n^2 + 1) |
	//

	gsl_matrix_set(A,0,0,1); 
	gsl_vector_set(b,0,Computef0(xi,modelConst,deltaEta)); 

	for(i =1; i<A->size1-1;i++)
	{
		xiCounter = 5*(i-1);
		
		LOvrEta = pow(ComputeL(xi,modelConst,i),2)/pow(deltaEta,2);
		gsl_matrix_set(A,i,i-1,LOvrEta); 
		gsl_matrix_set(A,i,i,-( 2*LOvrEta + 1)); 
		gsl_matrix_set(A,i,i+1,LOvrEta); 

		LHS1 = (modelConst->C1/ComputeT(xi,modelConst,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1)) - 2/3); 
		LHS2 = (modelConst->C2*ComputeP(xi,modelConst,deltaEta,i))/gsl_vector_get(xi,xiCounter+1);

		if(!isfinite(LHS1-LHS2))
		{
			cerr << "Error: non-finite (" << LHS1-LHS2 << ")"<<endl; 
			return 1; 
		}

		gsl_vector_set(b,i,LHS1-LHS2); 
	}
	
	// set boundary term of matrix. 
	i = size-1;	
	xiCounter = 5*(i-1);
	LOvrEta = pow(ComputeL(xi,modelConst,i),2)/pow(deltaEta,2);
	gsl_matrix_set(A,i,i-1,2*LOvrEta);
	gsl_matrix_set(A,i,i,-( 2*LOvrEta + 1)); 

	LHS1 = (modelConst->C1/ComputeT(xi,modelConst,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1)) - 2/3); 
	if (!isfinite(LHS1))
	{
		cerr << "Error: non-finite (" << LHS1 << ")" << endl; 
		return 1; 
	}

	gsl_vector_set(b,i,LHS1); 

	// LU solve Af = b for initial values f. 
	log(verbose,1," Performing LU Solve for f_0\n");
	int s; 
	gsl_permutation * p = gsl_permutation_alloc(A->size1);
	gsl_linalg_LU_decomp(A,p,&s); 
	gsl_linalg_LU_solve(A,p,b,f); 

	// add f to xi. 
	log(verbose,1," Setting f_0 values\n");
	for(unsigned int i =1; i<f->size; i++)
	{
		log(verbose,2,"  f = " + num2st(gsl_vector_get(f,i)) + " at " + num2st(i*deltaEta) + "\n");
		xiCounter=5*(i-1)+4; 
		gsl_vector_set(xi,xiCounter,gsl_vector_get(f,i));
	}
	return 0; 
}
