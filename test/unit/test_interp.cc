// g++ test_interp.cc  ../../src/options.cc -o test_interp -L$TACC_GRVY_LIB
// -lgrvy -I$TACC_GRVY_INC -L$TACC_GSL_LIB -lgsl -lgslcblas -I$TACC_GSL_INC 

#include<iostream>
#include"../../src/options.h"
#include<gsl/gsl_vector.h>

using namespace std; 

int main()
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

	SolveIC(U,k,ep,v2,deltaX,"temp.data");

	if (gsl_vector_equal(U,trueU)==0)
	{
		cout << "Test Failed" << endl; 
		return(1);
	}
	else if (gsl_vector_equal(k,truek)==0)
	{
		cout << "Test Failed" << endl;
		return(1);
	}

	else if (gsl_vector_equal(ep,trueEp)==0)
	{
		cout << "Test Failed" << endl;
		return(1);
	}
	else if (gsl_vector_equal(v2,truev2)==0)
	{
		cout << "Test Failed" << endl;
		return(1);
	}
	cout << "Interpolation Test Passed" << endl; 
	return 0;

}
