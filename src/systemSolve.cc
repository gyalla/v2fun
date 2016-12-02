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
	/*
	struct FParams * params = (struct FParams *)p; 
	int vecSize = ((xi->size))/double(5)+1; 
	
	gsl_vector * tempxi = gsl_vector_alloc(xi->size); 
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

	*for(int i = 0; i < k->size;i++)
	{
		cout << gsl_vector_get(U,i) << endl; 
	}
	for(int i = 0; i < k->size;i++)
	{
		cout << gsl_vector_get(k,i) << endl; 
	}
	for(int i = 0; i < k->size;i++)
	{
		cout << gsl_vector_get(ep,i) << endl; 
	}
	for(int i = 0; i < k->size;i++)
	{
		cout << gsl_vector_get(v2,i) << endl; 
	}
	for(int i = 0; i < k->size;i++)
	{
		cout << gsl_vector_get(f,i) << endl; 
	}*

	cout << "Computing Terms..." << endl;
	SetBdryEpf(ep,f,k,v2,params,sysF);
	ComputeT(k,ep,params->modelConst,T);
	ComputeL(k,ep,params->modelConst,L);
	ComputeEddyVisc(v2,T,params->modelConst,vT);
	ComputeP(U,vT,params->deltaEta,P); 

	//set boundary conditions here for ep and f 
	//set middle terms
	SetUTerms(U,vT,params,sysF); 
	SetKTerms(k,P,ep,vT,params,sysF);
	SetEpTerms(ep,k,P,T,vT,params,sysF);
	SetV2Terms(v2,k,f,vT,ep,params,sysF);
	SetFTerms(f,L,P,k,T,v2,params,sysF); 

	for(int i = 0; i<sysF->size;i++)
	{
		cout << gsl_vector_get(sysF,i) << endl; 
	}


	gsl_vector_free(tempxi);
	gsl_vector_free(U);
	gsl_vector_free(k);
	gsl_vector_free(ep);
	gsl_vector_free(v2);
	gsl_vector_free(f);
	gsl_vector_free(P);
	gsl_vector_free(T);
	gsl_vector_free(L);
	gsl_vector_free(vT);
	*/
	return 0; 
}
int SetFTerms(gsl_vector *f, gsl_vector * L, gsl_vector*P, gsl_vector * k, gsl_vector * T, gsl_vector * v2,  FParams * params, gsl_vector * sysF)
{
	/*
	double firstTerm,secondTerm,thirdTerm,fourthTerm; 
	unsigned int i; 
	double xiCounter; 

	for(i = 1; i<f->size-1; i++)
	{
		xiCounter=5*(i-1)+4; 
		firstTerm = -(gsl_vector_get(f,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = pow(gsl_vector_get(L,i),2)*Diff2(f,params->deltaEta,i); 
		thirdTerm = params->modelConst->C2*(gsl_vector_get(P,i)/gsl_vector_get(k,i)) - gsl_vector_get(f,i);     
		fourthTerm = -(params->modelConst->C1)/(gsl_vector_get(T,i))*( (gsl_vector_get(v2,i)/gsl_vector_get(k,i))-float(2)/3); 
		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i=f->size-1; 
	xiCounter=5*(i-1)+4; 
	firstTerm = -(gsl_vector_get(f,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = pow(gsl_vector_get(L,i),2)*BdryDiff2(f,params->deltaEta,i); 
	thirdTerm = -gsl_vector_get(f,i) -(params->modelConst->C1)/(gsl_vector_get(T,i))*( (gsl_vector_get(v2,i)/gsl_vector_get(k,i))-float(2)/3); 
	gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm); 

	return 0; 
	*/
}


int SetV2Terms(gsl_vector * v2, gsl_vector * k, gsl_vector * f, gsl_vector * vT, gsl_vector * ep,FParams * params, gsl_vector * sysF)
{
	/*
	double firstTerm,secondTerm,thirdTerm,fourthTerm; 
	unsigned int i; 
	double xiCounter; 

	for(i = 1; i<v2->size-1; i++)
	{
		xiCounter=5*(i-1)+3; 
		firstTerm = -(gsl_vector_get(v2,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = gsl_vector_get(k,i)*gsl_vector_get(f,i) - gsl_vector_get(ep,i)*( (gsl_vector_get(vT,i)/gsl_vector_get(k,i)));
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*Diff2(v2,params->deltaEta,i);
		fourthTerm = (1/params->modelConst->sigmaEp)*Diff1(v2,params->deltaEta,i)*Diff1(vT,params->deltaEta,i); 
		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i=v2->size-1; 
	xiCounter=5*(i-1)+3; 
	firstTerm = -(gsl_vector_get(v2,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = gsl_vector_get(k,i)*gsl_vector_get(f,i) - gsl_vector_get(ep,i)*( (gsl_vector_get(vT,i)/gsl_vector_get(k,i)));
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*BdryDiff2(v2,params->deltaEta,i);
	gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm); 

	return 0; 
	*/
}


int SetEpTerms(gsl_vector * ep, gsl_vector * k, gsl_vector * P, gsl_vector * T, gsl_vector * vT, FParams * params, gsl_vector * sysF)
{
	/*
	double firstTerm,secondTerm,thirdTerm,fourthTerm; 
	unsigned int i; 
	double xiCounter; 
	for (i = 1; i<ep->size-1;i++)
	{
		xiCounter=5*(i-1)+2; 
		firstTerm = -(gsl_vector_get(ep,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = (params->modelConst->Cep1*gsl_vector_get(P,i) - params->modelConst->Cep2*gsl_vector_get(ep,i))/gsl_vector_get(T,i); 
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*Diff2(ep,params->deltaEta,i);
		fourthTerm = (1/params->modelConst->sigmaEp)*Diff1(ep,params->deltaEta,i)*Diff1(vT,params->deltaEta,i); 
		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i=ep->size-1; 
	xiCounter=5*(i-1)+2; 
	firstTerm = -(gsl_vector_get(ep,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = - (params->modelConst->Cep2*gsl_vector_get(ep,i))/gsl_vector_get(T,i); 
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*BdryDiff2(ep,params->deltaEta,i);
	gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm); 
	return 0; 
	*/
}


int SetKTerms(gsl_vector * k,gsl_vector * P, gsl_vector * ep, gsl_vector* vT,FParams * params,gsl_vector *sysF)
{
	/*
	double firstTerm, secondTerm,thirdTerm,fourthTerm; 
	double xiCounter; 
	unsigned int i;  

	for(i=1; i<k->size-1;i++)
	{
		xiCounter=5*(i-1)+1;
		firstTerm = -(gsl_vector_get(k,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = gsl_vector_get(P,i) - gsl_vector_get(ep,i); 
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Diff2(k,params->deltaEta,i); 
		fourthTerm = Diff1(k,params->deltaEta,i)*Diff1(vT,params->deltaEta,i); 
		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+fourthTerm); 
	}

	i = k->size-1; 
	xiCounter = 5*(i-1)+1; 
	firstTerm = -(gsl_vector_get(k,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = gsl_vector_get(ep,i); 
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDiff2(k,params->deltaEta,i); 
	gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm); 
	
	return 0; 
	*/

}
int SetUTerms(gsl_vector * U, gsl_vector * vT, FParams * params,gsl_vector * sysF)
{
	/*double firstTerm, secondTerm, thirdTerm;
	double xiCounter;
	unsigned int i; 

	for(i=1; i<U->size-1; i++)
	{
		xiCounter = 5*(i-1);

		firstTerm = -(gsl_vector_get(U,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Diff2(U,params->deltaEta,i);
		thirdTerm = Diff1(U,params->deltaEta,i)*Diff1(U,params->deltaEta,i);

		gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+thirdTerm+1);
	}

	i = U->size-1; 
	xiCounter = 5*(i-1); 
	firstTerm = -(gsl_vector_get(U,i)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDiff2(U,params->deltaEta,i);
	gsl_vector_set(sysF,xiCounter,firstTerm+secondTerm+1); 	
	return 0; 
	*/
}

