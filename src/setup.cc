#include<cstdio>
#include<cstdlib>
#include<iostream>
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

int Grvy_Input_Parse(constants * modelConst,string & filename)
{
	GRVY_Input_Class iparse; // parsing object

	if(! iparse.Open("./input_file.txt"))
		exit(1);

	if(iparse.Read_Var("reyn", &(modelConst->reyn)))
		printf("---> %-11s = %f\n","reyn",(modelConst->reyn));
	
	if(iparse.Read_Var("Cmu", &(modelConst->Cmu)))
		printf("---> %-11s = %f\n","Cmu",(modelConst->Cmu));

	if(iparse.Read_Var("C1", &(modelConst->C1)))
		printf("---> %-11s = %f\n","C1",(modelConst->C1));

	if(iparse.Read_Var("C2", &(modelConst->C2)))
		printf("---> %-11s = %f\n","C2",(modelConst->C2));

	if(iparse.Read_Var("sigmaEp", &(modelConst->sigmaEp)))
		printf("---> %-11s = %f\n","sigmaEp",(modelConst->sigmaEp));

	if(iparse.Read_Var("CL", &(modelConst->CL)))
		printf("---> %-11s = %f\n","CL",(modelConst->CL));

	if(iparse.Read_Var("Cep1",&(modelConst->Cep1)))
		printf("---> %-11s = %f\n","Cep1",(modelConst->Cep1));

	if(iparse.Read_Var("Cep2", &(modelConst->Cep2)))
		printf("---> %-11s = %f\n","Cep2",(modelConst->Cep2));

	if(iparse.Read_Var("Ceta", &(modelConst->Ceta)))
		printf("---> %-11s = %f\n","Ceta",modelConst->Ceta);

	if(iparse.Read_Var("filename",&filename))
		cout << "---> filename = " << filename << endl; 

	iparse.Close();
	
	return 0;
}

void log(int verbose, int level, string message)
{
	if (verbose  >= level)
		cout << message; 
}

int LinInterp(gsl_vector * Vec,double pt1,double pt2, double U1,double U2,double gridPt,int i )
{
	double val = U1 + (gridPt-pt1)*( (U2-U1)/(pt2-pt1));
	gsl_vector_set(Vec,i,val); 
	return 0;
}

int SolveIC(gsl_vector * xi, double deltaEta, string file)
{	
	ifstream inFile; 
	inFile.open(file.c_str()); 
	
	double pt1,tempU1,tempK1,tempEp1,tempV21;
	double pt2,tempU2,tempK2,tempEp2,tempV22;
	double gridPt,xiCounter;
		
	if (!(inFile >> pt1 >> tempU1 >> tempK1 >> tempEp1 >> tempV21))
		return 1; 

	if (!(inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
		return 1; 
	for(unsigned int i=0; i<(xi->size/float(5));i++)
	{
		xiCounter=5*i;
		
		gridPt=(i+1)*deltaEta; 
		while ((! ((pt1 <= gridPt) && (pt2 >= gridPt))) && (!inFile.eof()))
		{
			pt1 = pt2; 
			tempU1 = tempU2;
			tempK1 = tempK2;
			tempEp1 = tempEp2; 
			tempV21 = tempV22; 
			if (! (inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
				cout << "error" <<endl; 
		}
		LinInterp(xi,pt1,pt2,tempU1,tempU2,gridPt,xiCounter); 	
		LinInterp(xi,pt1,pt2,tempK1,tempK2,gridPt,xiCounter+1); 	
		LinInterp(xi,pt1,pt2,tempEp1,tempEp2,gridPt,xiCounter+2); 	
		LinInterp(xi,pt1,pt2,tempV21,tempV22,gridPt,xiCounter+3); 	
	}

	inFile.close();

	return 0; 


}

int Solve4f0(gsl_vector * xi, constants * modelConst,double deltaEta)
{
	int size = xi->size/float(5) + 1; 
	gsl_matrix * A = gsl_matrix_calloc(size,size);
	gsl_vector * b = gsl_vector_calloc(A->size1); 
	gsl_vector * f = gsl_vector_calloc(A->size1);
	double LOvrEta; 
	double LHS1,LHS2; 
	unsigned int i,xiCounter; 

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
	
	i = size-1;	
	xiCounter = 5*(i-1);
	LOvrEta = pow(ComputeL(xi,modelConst,i),2)/pow(deltaEta,2);
	gsl_matrix_set(A,i,i-1,2*LOvrEta);
	gsl_matrix_set(A,i,i,-( 2*LOvrEta + 1)); 

	LHS1 = (modelConst->C1/ComputeT(xi,modelConst,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1)) - 2/3); 
	gsl_vector_set(b,i,LHS1); 

	int s; 
	gsl_permutation * p = gsl_permutation_alloc(A->size1);
	gsl_linalg_LU_decomp(A,p,&s); 
	gsl_linalg_LU_solve(A,p,b,f); 

	for(unsigned int i =1; i<f->size; i++)
	{
		xiCounter=5*(i-1)+4; 
		gsl_vector_set(xi,xiCounter,gsl_vector_get(f,i));
	}
	return 0; 

}
			
