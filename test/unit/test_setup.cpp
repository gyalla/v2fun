#include<iostream>
#include<grvy.h>
#include<fstream> 
#include<math.h>
#include"../../src/setup.h"
#include"test_computeTerms.h"

int test_interp()
{
	double deltaX = 0.25;
	Grid grid(true, 1.0, deltaX);
	gsl_vector * xi = gsl_vector_alloc(5*(1/deltaX));

	gsl_vector *truexi = gsl_vector_calloc(5*(1/deltaX)); 

	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst= & Const;  
	modelConst->reyn=1;
	modelConst->Cmu=3; 
	modelConst->C1=1; 
	modelConst->C2=2; 
	modelConst->Cep1=4; 
	modelConst->Cep2=5; 
	modelConst->Ceta=6; 
	modelConst->CL=7; 
	modelConst->sigmaEp=8;
	
	gsl_vector_set(truexi,0,0.5);
	gsl_vector_set(truexi,1,1);
	gsl_vector_set(truexi,3,-0.5);
	gsl_vector_set(truexi,5,1);
	gsl_vector_set(truexi,6,1);
	gsl_vector_set(truexi,8,-1);
	gsl_vector_set(truexi,10,1.5);
	gsl_vector_set(truexi,11,1);
	gsl_vector_set(truexi,13,-1.5);
	gsl_vector_set(truexi,15,2);
	gsl_vector_set(truexi,16,1);
	gsl_vector_set(truexi,18,-2);

	SolveIC(xi,modelConst,&grid,"./test_interp.data");

	if (gsl_vector_equal(xi,truexi)==0)
	{
		cout << "FAIL: Interpolation Test" << endl; 
		return 1; 
	}
	cout << "PASS: Interpolation Test" << endl; 
	return 0;
}

int test_Grvy_Input()
{ 	
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst= &Const; 
	string filename,outFile;
	bool uniformGrid;
	int max_ts; 
	if(Grvy_Input_Parse(modelConst,filename,outFile,uniformGrid,max_ts))
	{ 
		cout << "FAIL: Getting inputs" << endl; 
		return 1; 
	}

	if( (modelConst->reyn != 2000) || (modelConst->Cmu != 0.19) || (modelConst->C1 != 0.4)|| (modelConst->C2 != 0.3) || (modelConst->sigmaEp != 1.3) || (modelConst->CL != 0.3) || (modelConst->Cep2 != 1.9) || (modelConst->Cep1 != 1.55) || (modelConst->Ceta != 70) || (filename.compare("../data/Reyn_2000.dat")!=0) || (outFile.compare("../data/v2fResults_2000.dat")!=0) || (max_ts !=1))
	{
		cout << "FAIL: Getting inputs" << endl; 
		return 1; 
	}

	return 0; 
}

int test_Save_Results()
{
	gsl_vector * xi = gsl_vector_alloc(10); 
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst= &Const; 
	Setuptest(xi,modelConst);
	double deltaEta = 0.0625;
	Grid grid(true, 1, deltaEta);
	double temp; 

	SaveResults(xi,"test_output.txt",&grid,modelConst);
	ifstream inFile; 
	inFile.open("test_output.txt");
	inFile >> temp >> temp >> temp >> temp >> temp >> temp >> temp;  
	for(unsigned int i = 0; i<xi->size;i+=5)
	{
		inFile >> temp; 
		if (fabs(temp - (deltaEta*((i+5)/5)))>0.000001)
		{
			cout << "FAIL: Saving results" << endl; 
			return 1; 
		}
		for(unsigned int j=0;j<5;j++)
		{
			inFile >> temp; 
			if (fabs(temp - gsl_vector_get(xi,i+j))>0.000001)
			{
				cout <<"FAIL: Saving results" << endl; 
				return 1; 
			}
		}
	}
	cout << "PASS: Saving results" << endl; 
	inFile.close();
	return 0; 
}
