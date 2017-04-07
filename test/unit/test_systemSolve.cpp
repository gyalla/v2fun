#include<iostream>
#include<iomanip>
#include<math.h>
#include"../../src/systemSolve.h"
#include"../../src/computeTerms.h"
using namespace std; 

int Setuptest_SS(gsl_vector * xi, struct FParams * params)
{
	params->modelConst->reyn=2.0;
	params->modelConst->Cmu=3.0; 
	params->modelConst->C1=1.0; 
	params->modelConst->C2=2.0; 
	params->modelConst->Cep1=4.0; 
	params->modelConst->Cep2=5.0; 
	params->modelConst->Ceta=6.0; 
	params->modelConst->CL=7.0; 
	params->modelConst->sigmaEp=8.0;
	
	gsl_vector_set(xi,0,1.0);
	gsl_vector_set(xi,1,3.0);
	gsl_vector_set(xi,2,8.0);
	gsl_vector_set(xi,3,7.0);
	gsl_vector_set(xi,4,9.0);
	gsl_vector_set(xi,5,2.0);
	gsl_vector_set(xi,6,8.0);
	gsl_vector_set(xi,7,2.0);
	gsl_vector_set(xi,8,8.0);
	gsl_vector_set(xi,9,10.0);
	gsl_vector_set(xi,10,3.0);
	gsl_vector_set(xi,11,9.0);
	gsl_vector_set(xi,12,6.0);
	gsl_vector_set(xi,13,10.0);
	gsl_vector_set(xi,14,11.0);

	gsl_vector_set_all(params->XiN,1.0);

	return 0; 
}

int SetUTerms_test()
{
  Grid grid(true, 1.5, 0.5);
	gsl_vector * xi = gsl_vector_alloc(15); 
	gsl_vector * sysF = gsl_vector_alloc(15); 
	gsl_vector * trueF = gsl_vector_alloc(3);
	gsl_vector * xiN = gsl_vector_alloc(15); 
	gsl_vector * vT  = gsl_vector_calloc(4); 
	gsl_vector * T   = gsl_vector_calloc(4);
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst= &Const;  
	struct FParams p = {xiN,1.0,&grid,modelConst};
	FParams * params = &p; 
	Setuptest_SS(xi,params); 

	for (unsigned int i=1; i<vT->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(xi,modelConst,i));
		gsl_vector_set(vT,i,ComputeEddyVisc(xi,T,modelConst,i));
	}
	
	SetUTerms(xi,vT,params,sysF); 
	gsl_vector_set(trueF,0,193);
	gsl_vector_set(trueF,1,40.92304854); 
	gsl_vector_set(trueF,2,-420.69219382);

	for(unsigned int i =0; i<trueF->size;i++)
	{
		if(fabs(gsl_vector_get(sysF,5*i)-gsl_vector_get(trueF,i))>0.0000001)
		{
			cout <<"FAIL: Setting U terms in system" << endl; 
			return 1; 
		}
	}
	cout << "PASS: Setting U terms in system" << endl; 

	//uncomment to get values for tests
	//gsl_vector * P = gsl_vector_calloc(2);
	//gsl_vector * L = gsl_vector_calloc(4);

	//for (unsigned int i=1; i<L->size;i++)
	//{
	//	if (i!=3)
	//	cout << "P at " << i << " = " << ComputeP(xi,vT,0.5,i) << endl; 
	//	cout << "L at " << i << " = " << ComputeL(xi,modelConst,i) << endl; 
	//	cout << "T at " << i << " = " << gsl_vector_get(T,i) << endl; 
	//	cout << "vT at " << i << " = " << gsl_vector_get(vT,i) << endl; 
	//	 cout << "f0 = " << Computef0(xi,modelConst,0.5)<< endl; 
	//}

	return 0; 
}

int SetkTerms_test()
{
  Grid grid(true, 1.5, 0.5);
  gsl_vector * xi = gsl_vector_alloc(15);
  gsl_vector * sysF = gsl_vector_alloc(15);
  gsl_vector * trueF = gsl_vector_alloc(3);
  gsl_vector * xiN = gsl_vector_alloc(15);
  gsl_vector * vT  = gsl_vector_calloc(4);
  gsl_vector * T   = gsl_vector_calloc(4);
  struct constants Const = {
    .reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
  constants * modelConst= &Const;
  struct FParams p = {xiN,1.0,&grid,modelConst};
	FParams * params = &p; 
	Setuptest_SS(xi,params); 

	for (unsigned int i=1; i<vT->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(xi,params->modelConst,i));
		gsl_vector_set(vT,i,ComputeEddyVisc(xi,T,params->modelConst,i));
	}
	SetKTerms(xi,vT,params,sysF); 
	gsl_vector_set(trueF,0,1140);
	gsl_vector_set(trueF,1,-1046.23085463); 
	gsl_vector_set(trueF,2,-433.69219381);

	for(unsigned int i =0; i<trueF->size;i++)
	{
		if(fabs(gsl_vector_get(sysF,5*i+1)-gsl_vector_get(trueF,i))>0.0000001)
		{
			cout <<"FAIL: Setting k terms in system" << endl; 
			return 1; 
		}
	}
	cout << "PASS: Setting k terms in system" << endl; 
	return 0; 
}

int SetEpTerms_test()
{
  Grid grid(true, 1.5, 0.5);
  gsl_vector * xi = gsl_vector_alloc(15);
  gsl_vector * sysF = gsl_vector_alloc(15);
  gsl_vector * trueF = gsl_vector_alloc(3);
  gsl_vector * xiN = gsl_vector_alloc(15);
  gsl_vector * vT  = gsl_vector_calloc(4);
  gsl_vector * T   = gsl_vector_calloc(4);
  struct constants Const = {
    .reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
  constants * modelConst= &Const;
  struct FParams p = {xiN,1.0,&grid,modelConst};
	FParams * params = &p; 
	Setuptest_SS(xi,params); 

	for (unsigned int i=1; i<vT->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(xi,params->modelConst,i));
		gsl_vector_set(vT,i,ComputeEddyVisc(xi,T,params->modelConst,i));
	}

	SetEpTerms(xi,vT,T,params,sysF); 
	gsl_vector_set(trueF,0,146.83333333);
	gsl_vector_set(trueF,1,875.38461894); 
	gsl_vector_set(trueF,2,-246.166604983);

	for(unsigned int i =0; i<trueF->size;i++)
	{
		if(fabs(gsl_vector_get(sysF,5*i+2)-gsl_vector_get(trueF,i))>0.0000001)
		{
			cout <<"FAIL: Setting ep terms in system" << endl; 
			return 1; 
		}
	}
	cout << "PASS: Setting ep terms in system" << endl; 
	return 0; 
}	

int Setv2Terms_test()
{
  Grid grid(true, 1.5, 0.5);
  gsl_vector * xi = gsl_vector_alloc(15);
  gsl_vector * sysF = gsl_vector_alloc(15);
  gsl_vector * trueF = gsl_vector_alloc(3);
  gsl_vector * xiN = gsl_vector_alloc(15);
  gsl_vector * vT  = gsl_vector_calloc(4);
  gsl_vector * T   = gsl_vector_calloc(4);
  struct constants Const = {
    .reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
  constants * modelConst= &Const;
  struct FParams p = {xiN,1.0,&grid,modelConst};
	FParams * params = &p; 
	Setuptest_SS(xi,params); 

	for (unsigned int i=1; i<vT->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(xi,params->modelConst,i));
		gsl_vector_set(vT,i,ComputeEddyVisc(xi,T,params->modelConst,i));
	}

	
	SetV2Terms(xi,vT,params,sysF); 
	gsl_vector_set(trueF,0,2.33333333);
	gsl_vector_set(trueF,1,518.38457268119); 
	gsl_vector_set(trueF,2,-756.051054299727);

	for(unsigned int i =0; i<trueF->size;i++)
	{
		if(fabs(gsl_vector_get(sysF,5*i+3)-gsl_vector_get(trueF,i))>0.0000001)
		{
			cout <<"FAIL: Setting v2 terms in system" << endl; 
			return 1; 
		}
	}
	cout << "PASS: Setting v2 terms in system" << endl; 
	return 0; 
}	

int SetFTerms_test()
{
  Grid grid(true, 1.5, 0.5);
  gsl_vector * xi = gsl_vector_alloc(15);
  gsl_vector * sysF = gsl_vector_alloc(15);
  gsl_vector * trueF = gsl_vector_alloc(3);
  gsl_vector * xiN = gsl_vector_alloc(15);
  gsl_vector * vT  = gsl_vector_calloc(4);
  gsl_vector * T   = gsl_vector_calloc(4);
  struct constants Const = {
    .reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
  constants * modelConst= &Const;
  struct FParams p = {xiN,1.0,&grid,modelConst};
	FParams * params = &p; 
	Setuptest_SS(xi,params); 

	for (unsigned int i=1; i<T->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(xi,params->modelConst,i));
		gsl_vector_set(vT,i,ComputeEddyVisc(xi,T,params->modelConst,i));
	}

	SetFTerms(xi,vT,T,params,sysF); 
	gsl_vector_set(trueF,0,-27570.111111112);
	gsl_vector_set(trueF,1,76.916666666666); 
	gsl_vector_set(trueF,2,-7959.2566001196);

	for(unsigned int i =0; i<trueF->size;i++)
	{
		if(fabs(gsl_vector_get(sysF,5*i+4)-gsl_vector_get(trueF,i))>0.0000001)
		{
			cout <<"FAIL: Setting f terms in system" << endl; 
			return 1; 
		}
	}
	cout << "PASS: Setting f terms in system" << endl; 
	return 0; 
}	

int SysF_test()
{
  Grid grid(true, 1.5, 0.5);
	gsl_vector * xi = gsl_vector_alloc(15); 
	gsl_vector * F = gsl_vector_alloc(15); 
	gsl_vector * trueF = gsl_vector_alloc(15);
	gsl_vector * xiN = gsl_vector_alloc(15); 
	struct constants Const = {
		.reyn=0,.Cmu=0,.C1=0,.C2=0,.Cep1=0,.Cep2=0,.Ceta=0,.CL=0,.sigmaEp=0};
	constants * modelConst= &Const;  
	struct FParams p = {xiN,1.0,&grid,modelConst};
	FParams * params = &p; 
	Setuptest_SS(xi,params); 

	SysF(xi,params,F);

	//U TERMS
	gsl_vector_set(trueF,0,193);
	gsl_vector_set(trueF,5,40.92304854); 
	gsl_vector_set(trueF,10,-420.69219382);

	//K TERMS	
	gsl_vector_set(trueF,1,1140);
	gsl_vector_set(trueF,6,-1046.23085463); 
	gsl_vector_set(trueF,11,-433.69219381);

	//ep 
	gsl_vector_set(trueF,2,146.83333333);
	gsl_vector_set(trueF,7,875.38461894); 
	gsl_vector_set(trueF,12,-246.166604983);

	gsl_vector_set(trueF,3,2.33333333);
	gsl_vector_set(trueF,8,518.38457268119); 
	gsl_vector_set(trueF,13,-756.051054299727);

	gsl_vector_set(trueF,4,-27570.1111111112);
	gsl_vector_set(trueF,9,76.916666666666); 
	gsl_vector_set(trueF,14,-7959.2566001196);


	for(unsigned int i=0; i < trueF->size;i++)
	{
		if(fabs(gsl_vector_get(F,i)-gsl_vector_get(trueF,i))>0.0000001)
		{
			cout << "FAIL: Putting system together" << endl; 
			return 1; 
		}
	}
	cout << "PASS: Putting system together" << endl; 
	return 0; 
}









