#ifndef TEST_COMPUTETERMS_H
#define TEST_COMPUTETERMS_H

#include<gsl/gsl_vector.h>


int ComputeT_test();
int ComputeL_test();
int ComputeEddyVisc_test();
int ComputeEp0_test();
int Setuptest(gsl_vector * xi,constants * modelConst);
int ComputeP_test();

#endif
