#include<iostream>
#include<math.h>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>

#include"computeTerms.h"
#include"setup.h"
#include"systemSolve.h"
#include"finiteDiff.h"
using namespace std; 

int DeconstructXi(gsl_vector * xi, gsl_vector * U, gsl_vector * k,gsl_vector * ep,gsl_vector * v2, gsl_vector * f)
{
	int xiCounter;
	for(unsigned int i=1; i<U->size; i++)
	{
		xiCounter=5*(i-1); 
		gsl_vector_set(U,i,gsl_vector_get(xi,xiCounter));
		gsl_vector_set(k,i,gsl_vector_get(xi,xiCounter+1));
		gsl_vector_set(ep,i,gsl_vector_get(xi,xiCounter+2));
		gsl_vector_set(v2,i,gsl_vector_get(xi,xiCounter+3));
		gsl_vector_set(f,i,gsl_vector_get(xi,xiCounter+4));
	}

	return 0; 
}

int ReconstructXi(gsl_vector * xi,gsl_vector * U,gsl_vector *k,gsl_vector * ep,gsl_vector *v2,gsl_vector * f)
{
	int xiCounter; 
	for(unsigned int i=1; i<U->size;i++)
	{
		xiCounter=5*(i-1);
		gsl_vector_set(xi,xiCounter,gsl_vector_get(U,i));
		gsl_vector_set(xi,xiCounter+1,gsl_vector_get(k,i));
		gsl_vector_set(xi,xiCounter+2,gsl_vector_get(ep,i));
		gsl_vector_set(xi,xiCounter+3,gsl_vector_get(v2,i));
		gsl_vector_set(xi,xiCounter+4,gsl_vector_get(f,i));
	}
	return 0; 
}

int SysF(const gsl_vector * xi, void * p, gsl_vector * sysF)
{
	struct FParams * params = (struct FParams *)p; 
	int vecSize = ((xi->size)-2)/double(5)+1; 
	
	gsl_vector * tempxi;
	gsl_vector_memcpy(tempxi,xi);
	gsl_vector * U = gsl_vector_calloc(vecSize);
	gsl_vector * k = gsl_vector_calloc(vecSize);
	gsl_vector * ep = gsl_vector_calloc(vecSize);
	gsl_vector * v2 = gsl_vector_calloc(vecSize); 
	gsl_vector * f = gsl_vector_calloc(vecSize); 

	gsl_vector * P = gsl_vector_calloc(vecSize); 
	gsl_vector * T = gsl_vector_calloc(vecSize);
	gsl_vector * L = gsl_vector_calloc(vecSize); 
	gsl_vector * vT = gsl_vector_calloc(vecSize); 

	DeconstructXi(tempxi,U,k,ep,v2,f);
	ComputeT(k,ep,params->modelConst,T);
	ComputeL(k,ep,params->modelConst,L);
	ComputeEddyVisc(v2,T,params->modelConst,vT);
	ComputeP(U,vT,params->deltaEta,P); 


	//set middle terms
	SetUTerms(U,vT,params,sysF); 
	SetKTerms(k,P,ep,vT,params,sysF);
	SetEpTerms(ep,P,T,vT,params,sysF);
	SetV2Terms(v2,k,f,vT,ep,params,sysF);
	SetFTerms(f,L,P,k,T,v2,params,sysF); 

	return 0; 
}
int SetFTerms(gsl_vector *f, gsl_vector * L, gsl_vector*P, gsl_vector * k, gsl_vector * T, gsl_vector * v2,  FParams * params, gsl_vector * sysF)
{
	return 0; 
}


int SetV2Terms(gsl_vector * v2, gsl_vector * k, gsl_vector * f, gsl_vector * vT, gsl_vector * ep,FParams * params, gsl_vector * sysF)
{
	return 0; 
}


int SetEpTerms(gsl_vector * ep, gsl_vector * P, gsl_vector * T, gsl_vector * vT, FParams * params, gsl_vector * sysF)
{
	return 0; 
}


int SetUTerms(gsl_vector * U, gsl_vector * vT, FParams * params,gsl_vector * sysF)
{
	double firstTerm, secondTerm, thirdTerm;
	double xiCounter;
	for(unsigned int i=1; i<U->size-1; i++)
	{
		xiCounter = 5*(i-1);

		firstTerm = -(gsl_vector_get(U,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Diff2(U,params->deltaEta,i);
		thirdTerm = Diff1(U,params->deltaEta,i)*Diff1(U,params->deltaEta,i);

		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+1);
	}
	return 0; 
}

int SetKTerms(gsl_vector * k,gsl_vector * P, gsl_vector * ep, gsl_vector* vT,FParams * params,gsl_vector *sysF)
{
	return 0; 

}
	
int SetFirst2(gsl_vector * ep, gsl_vector * f,gsl_vector * k,  gsl_vector * v2, FParams * params,gsl_vector * sysF)
{

	double timeComp = (gsl_vector_get(ep,0)-gsl_vector_get(params->XiN,0))/params->deltaT;
	double spaceComp = (gsl_vector_get(ep,0) - ( (2*gsl_vector_get(k,1))/(params->modelConst->reyn*pow(params->deltaEta,2)))); 
	gsl_vector_set(sysF,0,timeComp+spaceComp); 

	timeComp = (gsl_vector_get(f,0)-gsl_vector_get(params->XiN,1))/params->deltaT;
	spaceComp = (gsl_vector_get(f,0) - ( (20*gsl_vector_get(v2,1))/( pow(params->modelConst->reyn,3)*gsl_vector_get(ep,0)*pow(params->deltaEta,4))));

	gsl_vector_set(sysF,1,timeComp+spaceComp); 

	return 0; 
}		
