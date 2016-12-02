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
/*	int xiCounter;
	for(unsigned int i=1; i<U->size; i++)
	{
		xiCounter=5*(i-1); 
		gsl_vector_set(U,i,gsl_vector_get(xi,xiCounter));
		gsl_vector_set(k,i,gsl_vector_get(xi,xiCounter+1));
		gsl_vector_set(ep,i,gsl_vector_get(xi,xiCounter+2));
		gsl_vector_set(v2,i,gsl_vector_get(xi,xiCounter+3));
		gsl_vector_set(f,i,gsl_vector_get(xi,xiCounter+4));
	}
*/
	return 0; 
}

int ReconstructXi(gsl_vector * xi,gsl_vector * U,gsl_vector *k,gsl_vector * ep,gsl_vector *v2,gsl_vector * f)
{
	/*
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
	*/
	return 0; 
}

int SysF(const gsl_vector * xi, void * p, gsl_vector * sysF)
{
	
	struct FParams * params = (struct FParams *)p; 
	int vecSize = ((xi->size))/double(5)+1; 
	
	gsl_vector * tempxi = gsl_vector_alloc(xi->size); 
	gsl_vector_memcpy(tempxi,xi);
	gsl_vector * vT = gsl_vector_calloc(vecSize);

	for (unsigned int i = 1; i<vT->size;i++)
	{
		gsl_vector_set(vT,i,ComputeEddyVisc(tempxi,params->modelConst,i));
	}
	SetUTerms(tempxi,vT,params,sysF); 
	SetKTerms(tempxi,vT,params,sysF);
	SetEpTerms(tempxi,vT,params,sysF);
	SetV2Terms(tempxi,vT,params,sysF);
	SetFTerms(tempxi,params,sysF); 

//	for(unsigned int i = 0; i<sysF->size;i++)
//	{
//		cout << gsl_vector_get(sysF,i) << endl; 
//	}


	gsl_vector_free(tempxi);

	return 0; 
}
int SetFTerms(gsl_vector * xi, FParams * params, gsl_vector * sysF)
{
	
	double firstTerm,secondTerm,thirdTerm,fourthTerm; 
	unsigned int i; 
	unsigned int size = xi->size/float(5) + 1; 
	double xiCounter; 

	for(i = 1; i<size-1; i++)
	{
		xiCounter=5*(i-1); 
		firstTerm = -(gsl_vector_get(xi,xiCounter+4)-gsl_vector_get(params->XiN,xiCounter+4))/params->deltaT;	
		secondTerm = pow(ComputeL(xi,params->modelConst,i),2)*Diff2(xi,params->deltaEta,Computef0(xi,params->modelConst,params->deltaEta),xiCounter+4); 
		thirdTerm = params->modelConst->C2*(ComputeP(xi,params->modelConst,params->deltaEta,i)/gsl_vector_get(xi,xiCounter+1)) - gsl_vector_get(xi,xiCounter+4);     
		fourthTerm = -(params->modelConst->C1)/(ComputeT(xi,params->modelConst,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1))-float(2)/3); 
		gsl_vector_set(sysF,xiCounter+4,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i=size-1;  
	xiCounter=5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+4)-gsl_vector_get(params->XiN,xiCounter+4))/params->deltaT;	
	secondTerm = pow(ComputeL(xi,params->modelConst,i),2)*BdryDiff2(xi,params->deltaEta,xiCounter+4); 
	thirdTerm = -gsl_vector_get(xi,xiCounter+4) -(params->modelConst->C1)/(ComputeT(xi,params->modelConst,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1))-float(2)/3); 
	gsl_vector_set(sysF,xiCounter+4,firstTerm+secondTerm+thirdTerm); 
	return 0; 
}


int SetV2Terms(gsl_vector * xi,gsl_vector * vT,FParams * params, gsl_vector * sysF)
{
	
	double firstTerm,secondTerm,thirdTerm,fourthTerm; 
	unsigned int i; 
	unsigned int size = xi->size/float(5)+1; 
	double xiCounter; 

	for(i = 1; i<size-1; i++)
	{
		xiCounter=5*(i-1); 
		firstTerm = -(gsl_vector_get(xi,xiCounter+3)-gsl_vector_get(params->XiN,xiCounter+3))/params->deltaT;	
		secondTerm = gsl_vector_get(xi,xiCounter+1)*gsl_vector_get(xi,xiCounter+4) - gsl_vector_get(xi,xiCounter+2)*( ( gsl_vector_get(vT,i)/gsl_vector_get(xi,xiCounter+1)));
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*Diff2(xi,params->deltaEta,0,xiCounter+3);
		fourthTerm = Diff1(xi,params->deltaEta,0,xiCounter+3)*Diff1vT(vT,params->deltaEta,i); 
		gsl_vector_set(sysF,xiCounter+3,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i=size-1; 
	xiCounter=5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+3)-gsl_vector_get(params->XiN,xiCounter+3))/params->deltaT;	
	secondTerm = gsl_vector_get(xi,xiCounter+1)*gsl_vector_get(xi,xiCounter+4) - gsl_vector_get(xi,xiCounter+2)*( (gsl_vector_get(vT,i)/gsl_vector_get(xi,xiCounter+1)));
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDiff2(xi,params->deltaEta,xiCounter+3);
	gsl_vector_set(sysF,xiCounter+3,firstTerm+secondTerm+thirdTerm); 

	return 0; 
}


int SetEpTerms(gsl_vector * xi, gsl_vector * vT,FParams * params, gsl_vector * sysF)
{
	double firstTerm,secondTerm,thirdTerm,fourthTerm; 
	unsigned int i; 
	unsigned int size = xi->size/float(5)+1; 
	double xiCounter; 

	for (i = 1; i<size-1;i++)
	{
		xiCounter=5*(i-1); 
		firstTerm = -(gsl_vector_get(xi,xiCounter+2)-gsl_vector_get(params->XiN,xiCounter+2))/params->deltaT;	
		secondTerm = (params->modelConst->Cep1*ComputeP(xi,params->modelConst,params->deltaEta,i) - params->modelConst->Cep2*gsl_vector_get(xi,xiCounter+2))/ComputeT(xi,params->modelConst,i); 
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*Diff2(xi,params->deltaEta,ComputeEp0(xi,params->modelConst,params->deltaEta),xiCounter+2);
		fourthTerm = (1/params->modelConst->sigmaEp)*Diff1(xi,params->deltaEta,ComputeEp0(xi,params->modelConst,params->deltaEta),xiCounter+2)*Diff1vT(vT,params->deltaEta,i); 
		gsl_vector_set(sysF,xiCounter+2,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i=size-1;  
	xiCounter=5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+2)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = - (params->modelConst->Cep2*gsl_vector_get(xi,xiCounter+2))/ComputeT(xi,params->modelConst,i); 
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*BdryDiff2(xi,params->deltaEta,xiCounter+2);
	gsl_vector_set(sysF,xiCounter+2,firstTerm+secondTerm+thirdTerm); 
	return 0; 
}


int SetKTerms(gsl_vector * xi, gsl_vector* vT,FParams * params,gsl_vector *sysF)
{
	double firstTerm, secondTerm,thirdTerm,fourthTerm; 
	unsigned int i;  
	unsigned int size=vT->size; 
	double xiCounter; 

	for(i=1; i<size-1;i++)
	{
		xiCounter=5*(i-1);
		firstTerm = -(gsl_vector_get(xi,xiCounter+1)-gsl_vector_get(params->XiN,xiCounter+1))/params->deltaT;	
		secondTerm = ComputeP(xi,params->modelConst,params->deltaEta,i)- gsl_vector_get(xi,xiCounter+2); 
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Diff2(xi,params->deltaEta,0,xiCounter+1); 
		fourthTerm = Diff1(xi,params->deltaEta,0,xiCounter+1)*Diff1vT(vT,params->deltaEta,i); 
		gsl_vector_set(sysF,xiCounter+1,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i = size-1;  
	xiCounter = 5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+1)-gsl_vector_get(params->XiN,xiCounter+1))/params->deltaT;	
	secondTerm = gsl_vector_get(xi,xiCounter+2); 
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDiff2(xi,params->deltaEta,xiCounter+1); 
	gsl_vector_set(sysF,xiCounter+1,firstTerm+secondTerm+thirdTerm); 
	
	return 0; 

}
int SetUTerms( gsl_vector * xi, gsl_vector * vT, FParams * params,gsl_vector * sysF)
{
	double firstTerm, secondTerm, thirdTerm;
	unsigned int i; 
	unsigned int size=vT->size; 
	double xiCounter;

	for(i=1; i<size-1; i++)
	{
		xiCounter = 5*(i-1);

		firstTerm = -(gsl_vector_get(xi,xiCounter)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Diff2(xi,params->deltaEta,0,xiCounter);
		thirdTerm = Diff1(xi,params->deltaEta,0,xiCounter)*Diff1vT(vT,params->deltaEta,i);

		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+1);
	}

	i =size-1; 
	xiCounter = 5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDiff2(xi,params->deltaEta,xiCounter);
	gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+1); 	
	return 0; 
	
}

