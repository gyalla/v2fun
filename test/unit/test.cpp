#include<iostream>
#include"../../src/setup.h"
#include"test_computeTerms.h"
#include"test_systemSolve.h"
#include"test_setup.h"
#include "test_finiteDiff.h"
using namespace std; 

int test_loglevel();
int main()
{
	cout << "--------------------------------------------------" << endl;
	cout << "Running Unit Test Suite" << endl; 
	cout << "--------------------------------------------------" << endl << endl; 

	test_interp();
	test_Grvy_Input();
	test_Save_Results();

	test_remapping();
	test_first_deriv(150);
	test_second_deriv(1200);

	ComputeT_test(); 
	ComputeL_test();
	ComputeEddyVisc_test();
	ComputeP_test();
	ComputeEp0_test();
	Computef0_test();

	SetUTerms_test(); 
	SetkTerms_test();
	SetEpTerms_test();
	Setv2Terms_test();
	SetFTerms_test();
	SysF_test();

	cout << "--------------------------------------------------" << endl << endl; 
	

}


