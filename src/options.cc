#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<grvy.h>
#include "options.h"
#include<gsl/gsl_vector.h>
#include<fstream>
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
	
	cout << "Opening File: "<<file<<endl;
	double pt1,tempU1,tempK1,tempEp1,tempV21;
	double pt2,tempU2,tempK2,tempEp2,tempV22;
	double gridPt;
		
	if (!(inFile >> pt1 >> tempU1 >> tempK1 >> tempEp1 >> tempV21))
		return 1; 

	if (!(inFile >> pt2 >> tempU2 >> tempK2 >> tempEp2 >> tempV22))
		return 1; 
	for(int i=0; i<U->size;i++)
	{
		cout << "Working on grid point" << i << "/"<< U->size << endl;
		gridPt=i*deltaEta; 
		cout << "Grid Point: " << gridPt << endl; 
		cout << "Point 1: " << pt1 << endl; 
		cout << "Point 2: " << pt2 << endl; 
		cout << "tempU1: " << tempU1<< endl; 
		cout << "tempU2: " << tempU2 << endl;

		while ((! ((pt1 <= gridPt) && (pt2 >= gridPt))) && (!inFile.eof()))
		{
			cout << "Here" << endl;
		cout << "Grid Point: " << gridPt << endl; 
		cout << "Point 1: " << pt1 << endl; 
		cout << "Point 2: " << pt2 << endl; 
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
