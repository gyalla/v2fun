#ifndef TEST_SYSTEMSOLVE_H
#define TEST_SYSTEMSOLVE_H

#include<gsl/gsl_vector.h>
#include"../../src/setup.h"
#include"../../src/systemSolve.h"
int SetUTerms_test();
int SetkTerms_test();
int SetEpTerms_test();
int Setv2Terms_test();
int SetFTerms_test();
int Setuptest_SS(gsl_vector *,struct FParams *);







#endif
