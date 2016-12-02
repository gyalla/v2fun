#include<iostream>
#include"../../src/setup.h"
#include"../../src/computeTerms.h"
#include"../../src/systemSolve.h"
#include<gsl/gsl_vector.h>

using namespace std; 

int test_interp();
int ComputeT_test();
int ComputeP_test();
int DeConstructXi_test();
int ReConstructXi_test();

int main()
{
	cout << "--------------------------------------------------" << endl;
	cout << "Running Unit Test Suite" << endl; 
	cout << "--------------------------------------------------" << endl << endl; 

	if (test_interp())
		cout << "PASS: Data Interpolation" << endl; 
	else
		cout << "FAIL: Data Interpolation"<<endl;

	/*
	if(ComputeT_test())
		cout << "PASS: Turbulent Time Scale Term" << endl; 
	else 
		cout << "FAIL: Turbulent Time Scale Term" << endl; 

	if(ComputeP_test())
		cout << "PASS: Production Rate Term" << endl; 
	else
		cout << "FAIL: Production Rate Term" << endl; 
	if(!DeConstructXi_test())
		cout << "PASS: Deconstructing Xi" << endl; 
	else
		cout << "FAIL: Deconstructing Xi" << endl;

	if(!ReConstructXi_test())
		cout << "PASS: Reconstructing Xi" << endl; 
	else
		cout << "FAIL: Reconstructing Xi" << endl; 
		*/
	return 0;

}


int test_interp()
{
	double deltaX = 0.25; 
	gsl_vector * xi = gsl_vector_alloc(5*(1/deltaX));


	gsl_vector *truexi = gsl_vector_calloc(5*(1/deltaX)); 
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

	SolveIC(xi,deltaX,"test_interp.data");

	if (gsl_vector_equal(xi,truexi)==0)
		return 0;
	return 1;
}
/*
int ComputeT_test()
{
	constants Const = {
		.reyn=2,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst = &Const; 
	gsl_vector * k = gsl_vector_alloc(2); 
	gsl_vector * ep = gsl_vector_alloc(2); 
	gsl_vector * T = gsl_vector_alloc(2); 
	gsl_vector * trueT = gsl_vector_alloc(2); 

	gsl_vector_set(k,0,6); 
	gsl_vector_set(k,1,2);
	gsl_vector_set(ep,0,0.5);
	gsl_vector_set(ep,1,2); 

	gsl_vector_set(trueT,0,12);
	gsl_vector_set(trueT,1,3);

	if(!ComputeT(k,ep,modelConst,T))
		return 0; 

	if(!gsl_vector_equal(T,trueT))
		return 0;
	return 1; 
}

int ComputeP_test()
{	
	
	double deltaEta = 0.25; 
	gsl_vector * U = gsl_vector_calloc(3); 
	gsl_vector * vt = gsl_vector_alloc(3); 
	gsl_vector * P = gsl_vector_calloc(3); 
	gsl_vector * trueP = gsl_vector_calloc(3); 

	gsl_vector_set(trueP,0,16);
	gsl_vector_set(trueP,1,24); 
	gsl_vector_set(trueP,2,0);

	gsl_vector_set(U,0,0);
	gsl_vector_set(U,1,4);
	gsl_vector_set(U,2,6); 

	gsl_vector_set(vt,0,1);
	gsl_vector_set(vt,1,2);
	gsl_vector_set(vt,2,3);

	if(!ComputeP(U,vt,deltaEta,P))
		return 0; 
	if(!gsl_vector_equal(P,trueP))
		return 0; 
	return 1; 
}

int DeConstructXi_test()
{
	gsl_vector * xi = gsl_vector_calloc(10);
	gsl_vector * U  = gsl_vector_calloc(3);
	gsl_vector * k  = gsl_vector_calloc(3);
	gsl_vector * ep = gsl_vector_calloc(3);
	gsl_vector * v2 = gsl_vector_calloc(3);
	gsl_vector * f  = gsl_vector_calloc (3);
	
	gsl_vector *trueU = gsl_vector_calloc(3); 
	gsl_vector *truek = gsl_vector_calloc(3);
	gsl_vector *trueep = gsl_vector_calloc(3);
	gsl_vector *truev2 = gsl_vector_calloc(3);
	gsl_vector *truef = gsl_vector_calloc(3);


	gsl_vector_set(xi,0,1);
	gsl_vector_set(xi,1,-1);
	gsl_vector_set(xi,2,1);
	gsl_vector_set(xi,3,0);
	gsl_vector_set(xi,4,5);
	gsl_vector_set(xi,5,2);
	gsl_vector_set(xi,6,-2);
	gsl_vector_set(xi,7,1);
	gsl_vector_set(xi,8,0);
	gsl_vector_set(xi,9,6);
	
	gsl_vector_set(trueU,1,1);
	gsl_vector_set(trueU,2,2);

	gsl_vector_set(truek,1,-1);
	gsl_vector_set(truek,2,-2);

	gsl_vector_set(trueep,0,0);
	gsl_vector_set(trueep,1,1);
	gsl_vector_set(trueep,2,1);
	
	gsl_vector_set(truef,0,0);
	gsl_vector_set(truef,1,5);
	gsl_vector_set(truef,2,6);

	if(DeconstructXi(xi,U,k,ep,v2,f))
		return 1;
	if(!gsl_vector_equal(U,trueU))
		return 1; 
	if(!gsl_vector_equal(k,truek))
		return 1;
	if(!gsl_vector_equal(ep,trueep))
		return 1; 
	if(!gsl_vector_equal(v2,truev2))
		return 1;
	if(!gsl_vector_equal(f,truef))
		return 1; 
	return 0; 
}

int ReConstructXi_test()
{

	gsl_vector * xi = gsl_vector_calloc(10);
	gsl_vector * U  = gsl_vector_calloc(3);
	gsl_vector * k  = gsl_vector_calloc(3);
	gsl_vector * ep = gsl_vector_calloc(3);
	gsl_vector * v2 = gsl_vector_calloc(3);
	gsl_vector * f  = gsl_vector_calloc (3);
	
	gsl_vector * truexi = gsl_vector_calloc(10);
	gsl_vector_set(truexi,0,1);
	gsl_vector_set(truexi,1,-1);
	gsl_vector_set(truexi,2,1);
	gsl_vector_set(truexi,3,0);
	gsl_vector_set(truexi,4,5);
	gsl_vector_set(truexi,5,2);
	gsl_vector_set(truexi,6,-2);
	gsl_vector_set(truexi,7,1);
	gsl_vector_set(truexi,8,0);
	gsl_vector_set(truexi,9,6);
	

	gsl_vector_set(U,1,1);
	gsl_vector_set(U,2,2);

	gsl_vector_set(ep,0,0);
	gsl_vector_set(ep,1,1);
	gsl_vector_set(ep,2,1);
	
	
	gsl_vector_set(k,1,-1);
	gsl_vector_set(k,2,-2);
	
	gsl_vector_set(f,0,0);
	gsl_vector_set(f,1,5);
	gsl_vector_set(f,2,6);
	
	if(ReconstructXi(xi,U,k,ep,v2,f))
		return 1; 
	if(!gsl_vector_equal(xi,truexi))
		return 1; 
	return 0; 

}
*/
