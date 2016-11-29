#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<grvy.h>
#include "options.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<fstream>
#include<math.h>
using namespace std;
using namespace GRVY;

int Grvy_Input_Parse(double & reyn, double & Cmu, double & C1, double & C2, double & sigmaEp, double & CL, double & Cep2, double & Cep1, double & Ceta, string & filename)
{
	GRVY_Input_Class iparse; // parsing object

	if(! iparse.Open("./input_file.txt"))
		exit(1);

	if(iparse.Read_Var("reyn", &reyn))
		printf("---> %-11s = %f\n","reyn",reyn);
	
	if(iparse.Read_Var("Cmu", &Cmu))
		printf("---> %-11s = %f\n","Cmu",Cmu);

	if(iparse.Read_Var("C1", &C1))
		printf("---> %-11s = %f\n","C1",C1);

	if(iparse.Read_Var("C2", &C2))
		printf("---> %-11s = %f\n","C2",C2);

	if(iparse.Read_Var("sigmaEp", &sigmaEp))
		printf("---> %-11s = %f\n","sigmaEp",sigmaEp);

	if(iparse.Read_Var("CL", &CL))
		printf("---> %-11s = %f\n","CL",CL);

	if(iparse.Read_Var("Cep1", &Cep1))
		printf("---> %-11s = %f\n","Cep1",Cep1);

	if(iparse.Read_Var("Cep2", &Cep2))
		printf("---> %-11s = %f\n","Cep2",Cep2);

	if(iparse.Read_Var("Ceta", &Ceta))
		printf("---> %-11s = %f\n","Ceta",Ceta);

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

int LinInterp(gsl_vector * Vec,double pt1,double pt2, double tempU1,double tempU2,double gridPt,int i )
{
	double val = tempU1 + (gridPt-pt1)*( (tempU2-tempU1)/(pt2-pt1));
	gsl_vector_set(Vec,i,val); 
	return 0;
}

int SolveIC(gsl_vector * U, gsl_vector * k, gsl_vector * ep, gsl_vector * v2,double deltaEta,string file)
{	
	ifstream inFile; 
	inFile.open(file.c_str()); 
	
	double pt1,tempU1,tempK1,tempEp1,tempV21;
	double pt2,tempU2,tempK2,tempEp2,tempV22;
	double gridPt;
		
	if (!(inFile >> pt1 >> tempU1 >> tempK1 >> tempEp1 >> tempV21))
		return 1; 

	if (!(inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
		return 1; 
	for(int i=0; i<U->size;i++)
	{
		gridPt=i*deltaEta; 
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
		LinInterp(U,pt1,pt2,tempU1,tempU2,gridPt,i); 	
		LinInterp(k,pt1,pt2,tempK1,tempK2,gridPt,i); 	
		LinInterp(ep,pt1,pt2,tempEp1,tempEp2,gridPt,i); 	
		LinInterp(v2,pt1,pt2,tempV21,tempV22,gridPt,i); 	
	}

	inFile.close();

	return 0; 


}

int Solve4f0(gsl_vector * k, gsl_vector * ep, gsl_vector * v2, gsl_vector * P, gsl_vector *  T,gsl_vector * L, double reyn, double C2,double C1, double deltaEta, gsl_vector * f)
{
	gsl_matrix * A = gsl_matrix_calloc(k->size,k->size);
	gsl_vector * b = gsl_vector_calloc(A->size1); 
	double val; 
	double LOvrEta; 
	double LHS1,LHS2; 
	int i; 

	val = (20*gsl_vector_get(v2,1))/( pow(reyn,3)*gsl_vector_get(ep,0)*pow(deltaEta,4));  
	if (!isfinite(val))
	{
		cout <<"Error: non-finite (" << val << ")" << endl; 
		return 1; 
	}
	gsl_matrix_set(A,0,0,1); 
	gsl_vector_set(b,0,val); 

	for(i =1; i<A->size1-1;i++)
	{
		LOvrEta = pow(gsl_vector_get(L,i),2)/pow(deltaEta,2);
		gsl_matrix_set(A,i,i-1,LOvrEta); 
		gsl_matrix_set(A,i,i,-( 2*LOvrEta + 1)); 
		gsl_matrix_set(A,i,i+1,LOvrEta); 

		LHS1 = (C1/gsl_vector_get(T,i))*( (gsl_vector_get(v2,i)/gsl_vector_get(k,i)) - 2/3); 
		LHS2 = (C2*gsl_vector_get(P,i))/gsl_vector_get(k,i);

		if(!isfinite(LHS1-LHS2))
		{
			cerr << "Error: non-finite (" << LHS1-LHS2 << ")"<<endl; 
			return 1; 
		}

		gsl_vector_set(b,i,LHS1-LHS2); 
	}
	
	i = L->size-1;	
	LOvrEta = pow(gsl_vector_get(L,i),2)/pow(deltaEta,2);
	gsl_matrix_set(A,i,i-1,2*LOvrEta);
	gsl_matrix_set(A,i,i,-( 2*LOvrEta + 1)); 

	LHS1 = (C1/gsl_vector_get(T,i))*( (gsl_vector_get(v2,i)/gsl_vector_get(k,i)) - 2/3); 
	gsl_vector_set(b,i,LHS1); 

	int s; 
	gsl_permutation * p = gsl_permutation_alloc(A->size1);
	gsl_linalg_LU_decomp(A,p,&s); 
	gsl_linalg_LU_solve(A,p,b,f); 

	return 0; 

}
			
