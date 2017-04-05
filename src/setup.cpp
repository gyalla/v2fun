//--------------------------------------------------
// setup: Initializes terms for v2-f code. 
//
// 12/3/2016 (gry88) - Written for CSE380 Final Project. 
//--------------------------------------------------
#include<iomanip>
#include<grvy.h>
#include "setup.h"
#include "computeTerms.h"
#include<gsl/gsl_linalg.h>
#include<fstream>
#include<math.h>
using namespace std;
using namespace GRVY;
loglevel_e loglevel = logINFO;


int Grvy_Input_Parse(constants * modelConst,string & filename,string & outFile,double & deltaEta, bool &uniformGrid)
{
	// Use GRVY for inpute parsing as shown in grvy documentation. 
	int loglevelint;  
	
	GRVY_Input_Class iparse; // parsing object

	if(!iparse.Open("./input_file.txt"))
		exit(1);
	
	// read in each variable. 
	
	if(!iparse.Read_Var("loglevelint",&loglevelint))
		return 1; 

	loglevel=(loglevel_e)loglevelint; //tpye case int as loglevel

	if(!iparse.Read_Var("reyn", &(modelConst->reyn)))
		return 1; 
	
	Log(logINFO) << "---> reyn = " << modelConst->reyn;
	//log(verbose,1," ---> reyn = " + num2st(modelConst->reyn,1) + "\n");

	if(!iparse.Read_Var("Cmu", &(modelConst->Cmu)))
		return 1; 
	Log(logINFO) << "---> Cmu = " << modelConst->Cmu;

	if(!iparse.Read_Var("C1", &(modelConst->C1)))
		return 1; 
	Log(logINFO) <<"---> C1 = " << modelConst->C1;

	if(!iparse.Read_Var("C2", &(modelConst->C2)))
		return 1; 
	Log(logINFO) << "---> C2 = " << modelConst->C2;

	if(!iparse.Read_Var("sigmaEp", &(modelConst->sigmaEp)))
		return 1; 
	Log(logINFO) << "---> sigmaEp = " << modelConst->sigmaEp;

	if(!iparse.Read_Var("CL", &(modelConst->CL)))
		return 1; 
	Log(logINFO) << "---> CL = " << modelConst->CL;

	if(!iparse.Read_Var("Cep1",&(modelConst->Cep1)))
		return 1; 
	Log(logINFO) << "---> Cep1 = " << modelConst->CL;

	if(!iparse.Read_Var("Cep2", &(modelConst->Cep2)))
		return 1; 
	Log(logINFO) <<"---> Cep2 = " << modelConst->Cep2;

	if(!iparse.Read_Var("Ceta", &(modelConst->Ceta)))
		return 1; 
	Log(logINFO) << "---> Ceta = " << modelConst->Ceta;
	
	if(!iparse.Read_Var("data_filename",&filename))
		return 1; 
	Log(logINFO) << "---> data_filename = " << filename;

	if(!iparse.Read_Var("output_filename",&outFile))
		return 1; 
	Log(logINFO) << "---> output_filename = " << outFile;

	if(!iparse.Read_Var("deltaEta",&deltaEta))
		return 1; 
	Log(logINFO) << "---> deltaEta = " << deltaEta;

  if(!iparse.Read_Var("uniform-grid",&uniformGrid, true))
    return 1;
  Log(logINFO) << "---> deltaEta = " << uniformGrid;

	iparse.Close();
	
	return 0;
}

int LinInterp(gsl_vector * Vec,double pt1,double pt2, double U1,double U2,double gridPt,int i )
{
	double val = U1 + (gridPt-pt1)*( (U2-U1)/(pt2-pt1));
	if (!isfinite(val))
	{
		Log(logERROR) <<  "Error: non-finite interpolation";
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
		
	// read in first two line. 
	if (!(inFile >> pt1 >> tempU1 >> tempK1 >> tempEp1 >> tempV21))
		return 1; 

	if (!(inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
		return 1; 
	for(unsigned int i=0; i<(xi->size/float(5));i++)
	{
		xiCounter=5*i; //counter relative to xi vector. 
		gridPt=(i+1)*deltaEta; //which grid point we are working on. 
		Log(logDEBUG1) <<"Working on grid point: " << gridPt;
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
			{
				Log(logERROR) << "Error: Cannot get interpolation data."; 
			}
		}

		//interpolate each term U,k,ep,v2,f.
		Log(logDEBUG2) << "interpolation points: " << pt1 << "," <<  pt2;
		Log(logDEBUG2) << "interpolation values (U): " <<  tempU1 <<  "," << tempU2;
		Log(logDEBUG2) << "interpolation values (k): " <<tempK1<< "," << tempK2;
		Log(logDEBUG2) << "interpolation values (ep): " << tempEp1<< "," << tempEp1;
		Log(logDEBUG2) << "interpolation values (v2): " << tempV21 <<  "," << tempV21;
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

	gsl_vector * T = gsl_vector_calloc(A->size1);
	gsl_vector * vT = gsl_vector_calloc(T->size);
	for (unsigned int i = 1; i<vT->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(xi,modelConst,i));
		gsl_vector_set(vT,i,ComputeEddyVisc(xi,T,modelConst,i));
	}

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

		LHS1 = (modelConst->C1/gsl_vector_get(T,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1)) - 2/3); 
		LHS2 = (modelConst->C2*ComputeP(xi,vT,deltaEta,i))/gsl_vector_get(xi,xiCounter+1);

		if(!isfinite(LHS1-LHS2))
		{
			Log(logERROR) << "Error: non-finite b (" << LHS1-LHS2 << ")"; 
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
		Log(logERROR) << "Error: non-finite b (" << LHS1 << ")"; 
		return 1; 
	}

	gsl_vector_set(b,i,LHS1); 

	// LU solve Af = b for initial values f. 
	Log(logDEBUG) << "Performing LU Solve for f_0";
	int s; 
	gsl_permutation * p = gsl_permutation_alloc(A->size1);
	gsl_linalg_LU_decomp(A,p,&s); 
	gsl_linalg_LU_solve(A,p,b,f); 

	// add f to xi. 
	Log(logDEBUG) << "Setting f_0 values";
	for(unsigned int i =1; i<f->size; i++)
	{
		Log(logDEBUG2) <<"f = " << gsl_vector_get(f,i)<< " at " << i*deltaEta;
		xiCounter=5*(i-1)+4; 
		gsl_vector_set(xi,xiCounter,gsl_vector_get(f,i));
	}

	return 0; 
}

void SaveResults(gsl_vector * xi, string filename,double deltaEta,constants * modelConst)
{
	ofstream outFile; 
	outFile.open(filename.c_str()); 
	//output format: gridpoint U K EP V2 F 
	outFile << setprecision(15) << 0 <<" "<<0<<" "<< 0 <<" "<< ComputeEp0(xi,modelConst,deltaEta) <<" " << 0 << " " << Computef0(xi,modelConst,deltaEta) << endl; 

	for(unsigned int i = 0; i<xi->size;i+=5)
	{
		outFile << (i/float(5)+1)*deltaEta << " "; 	
		outFile << gsl_vector_get(xi,i) << " "; 
		outFile << gsl_vector_get(xi,i+1) << " "; 
		outFile << gsl_vector_get(xi,i+2) << " "; 
		outFile << gsl_vector_get(xi,i+3) << " "; 
		outFile << gsl_vector_get(xi,i+4) << " " << endl; 
	}

	outFile.close();

}
