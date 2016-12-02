#include<iostream>
#include<gsl/gsl_vector.h>
#include"../../src/setup.h"
#include"test_computeTerms.h"
using namespace std; 

int test_interp();
int main()
{
	cout << "--------------------------------------------------" << endl;
	cout << "Running Unit Test Suite" << endl; 
	cout << "--------------------------------------------------" << endl << endl; 

	test_interp();
	ComputeT_test(); 
	ComputeL_test();
	ComputeEddyVisc_test();
	ComputeP_test();
	ComputeEp0_test();
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

	SolveIC(xi,deltaX,"../test/unit/test_interp.data");

	if (gsl_vector_equal(xi,truexi)==0)
	{
		cout << "FAIL: Interpolation Test" << endl; 
		return 1; 
	}
	cout << "PASS: Interpolation Test" << endl; 
	return 0;
}
