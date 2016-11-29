#include<iostream>
#include"../../src/setup.h"
#include"../../src/computeTerms.h"
#include<gsl/gsl_vector.h>


using namespace std; 

int test_interp();
int ComputeT_test();
int ComputeP_test();

int main()
{
	cout << "--------------------------------------------------" << endl;
	cout << "Running Unit Test Suite" << endl; 
	cout << "--------------------------------------------------" << endl << endl; 

	if (test_interp())
		cout << "PASS:  Data Interpolation" << endl; 
	else
		cout << "FAIL:  Data Interpolation"<<endl;

	if(ComputeT_test())
		cout << "PASS: Turbulent Time Scale Term" << endl; 
	else 
		cout << "FAIL: Turbulent Time Scale Term" << endl; 

	if(ComputeP_test())
		cout << "PASS: Production Rate Term" << endl; 
	else
		cout << "FAIL: Production Rate Term" << endl; 
	return 0;

}


int test_interp()
{
	double deltaX = 0.25; 
	gsl_vector * U = gsl_vector_alloc(1/deltaX + 1);
	gsl_vector * k = gsl_vector_alloc(1/deltaX + 1);
	gsl_vector * ep = gsl_vector_alloc(1/deltaX + 1);
	gsl_vector * v2 = gsl_vector_alloc(1/deltaX + 1);


	gsl_vector *trueU = gsl_vector_calloc(1/deltaX+1); 
	gsl_vector_set(trueU,1,0.5);
	gsl_vector_set(trueU,2,1);
	gsl_vector_set(trueU,3,1.5);
	gsl_vector_set(trueU,4,2); 

	gsl_vector * truek = gsl_vector_alloc(1/deltaX+1);
	gsl_vector_set_all(truek,1); 

	gsl_vector * trueEp = gsl_vector_calloc(1/deltaX+1);
	
	gsl_vector *truev2 = gsl_vector_calloc(1/deltaX+1); 
	gsl_vector_set(truev2,1,-0.5);
	gsl_vector_set(truev2,2,-1);
	gsl_vector_set(truev2,3,-1.5);
	gsl_vector_set(truev2,4,-2); 

	SolveIC(U,k,ep,v2,deltaX,"test_interp.data");

	if (gsl_vector_equal(U,trueU)==0)
		return 0;
	if (gsl_vector_equal(k,truek)==0)
		return 0;
	if (gsl_vector_equal(ep,trueEp)==0)
		return 0;
	if (gsl_vector_equal(v2,truev2)==0)
		return 0;
	return 1;
}

int ComputeT_test()
{
	int reyn = 2; 
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

	if(!ComputeT(k,ep,reyn,T))
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
