//--------------------------------------------------
// systemSolve: Sets up system to be solved for gsl_multiroot solver. 
//
// 12/3/2016 - (gry88) Written for CSE380 final project. 
//--------------------------------------------------
#include<math.h>
#include<gsl/gsl_multiroots.h>
#include<omp.h>
#include<grvy.h>
#include"computeTerms.h"
#include"systemSolve.h"
#include"finiteDiff.h"
using namespace std; 
using namespace GRVY;

// The structure of this function is fixed by the definition of gsl_multiroot solvers. 
int SysF(const gsl_vector * xi, void * p, gsl_vector * sysF)
{
	Log(logDEBUG2) << "Setting up system F(xi)";
	struct FParams * params = (struct FParams *)p; //reference void pointer to parameter struct; 

	int vecSize = ((xi->size))/double(5)+1;  // size of single vector. I in doc. 
	
	//first paramter must be const. But since we need do access vector, we have to define
	//a temporary vector so work with. Perhaps there is a way around this. 
	gsl_vector * tempxi = gsl_vector_alloc(xi->size);  
	gsl_vector_memcpy(tempxi,xi);

	// Enough functions need eddy viscosity so we store results, same as T.
	// Note each of the Term vectors are full size (starting at i=0)
	gsl_vector * vT = gsl_vector_calloc(vecSize);
	gsl_vector * T  = gsl_vector_calloc(vecSize);
	for (unsigned int i = 1; i<vT->size;i++)
	{
		gsl_vector_set(T,i,ComputeT(tempxi,params->modelConst,i));
                // Set to 0 for laminar case
		gsl_vector_set(vT,i,ComputeEddyVisc(tempxi,T,params->modelConst,i));
	}

	GRVY_Timer_Class gt; 
	gt.BeginTimer("Setting up system");
	//Set each term based on functions below. 
	if(SetUTerms(tempxi,vT,params,sysF))
	{
		Log(logERROR) << "Error setting U terms in system";
		exit(1);
	}

	//cout << omp_get_thread_num() << " MADE IT" << endl;

	if(SetKTerms(tempxi,vT,params,sysF))
	{
		Log(logERROR) << "Error setting k terms in system";
		exit(1); 
	}

	if(SetEpTerms(tempxi,vT,T,params,sysF))
	{
		Log(logERROR) << "Error setting ep terms in system";
		exit(1); 
	}

	if(SetV2Terms(tempxi,vT,params,sysF))
	{
		Log(logERROR) << "Error setting v2 terms in system";
		exit(1); 
	}

	if(SetFTerms(tempxi,vT,T,params,sysF))
	{
		Log(logERROR) << "Error setting F terms in system";
		exit(1); 
	}

	gt.EndTimer("Setting up system");
//	for(unsigned int i = 0; i<sysF->size;i++)
//	{
//		cout << gsl_vector_get(sysF,i) << endl; 
//	}
	
	// Cleanup
	gsl_vector_free(vT);
	gsl_vector_free(T);
	gsl_vector_free(tempxi);

	return 0; 
}

int SetFTerms(gsl_vector * xi, gsl_vector * vT, gsl_vector * T, FParams * params, gsl_vector * sysF)
{
	Log(logDEBUG2) << "Setting f terms\n";
	unsigned int i; 
	unsigned int size = xi->size/float(5) + 1; //size of single vectors. Comes from old structure of code before restructure branch in git.  
	double xiCounter; 
	double f0 = Computef0(xi,params->modelConst,params->grid);

	#pragma omp parallel num_threads(4)
	{
	#pragma omp for private(xiCounter)
	for(i = 1; i<size-1; i++)
	{
		double firstTerm,secondTerm,thirdTerm,fourthTerm,val; 
		xiCounter=5*(i-1); 
		firstTerm = -(gsl_vector_get(xi,xiCounter+4)-gsl_vector_get(params->XiN,xiCounter+4))/params->deltaT;	
		secondTerm = pow(ComputeL(xi,params->modelConst,i),2)*Deriv2(xi,f0,xiCounter+4, params->grid);
		thirdTerm = params->modelConst->C2*(ComputeP(xi,vT,params->grid,i)/gsl_vector_get(xi,xiCounter+1)) - gsl_vector_get(xi,xiCounter+4);
		fourthTerm = -(params->modelConst->C1/gsl_vector_get(T,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1))-float(2)/3); 
		val = firstTerm + secondTerm + thirdTerm + fourthTerm;
		//if (!isfinite(val))
		//	return 1; 
		gsl_vector_set(sysF,xiCounter+4,val); 
		Log(logDEBUG3) << "f term = " << val<< " at " << i;
	}
	}

	double firstTerm,secondTerm,thirdTerm,fourthTerm,val; 
	//boundary terms. 
	i=size-1;  
	xiCounter=5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+4)-gsl_vector_get(params->XiN,xiCounter+4))/params->deltaT;	
	secondTerm = pow(ComputeL(xi,params->modelConst,i),2)*BdryDeriv2(xi,xiCounter+4,params->grid);
	thirdTerm = -gsl_vector_get(xi,xiCounter+4) -(params->modelConst->C1/gsl_vector_get(T,i))*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1))-float(2)/3); 
	val = firstTerm+secondTerm+thirdTerm;
	Log(logDEBUG3) << "f term = " << val << " at " << i;
	if(!isfinite(val))
		return 1; 
	gsl_vector_set(sysF,xiCounter+4,val);
	return 0; 
} 
int SetV2Terms(gsl_vector * xi,gsl_vector * vT,FParams * params, gsl_vector * sysF)
{
	
	Log(logDEBUG2) << "Setting V2 terms";
	unsigned int i; 
	unsigned int size = xi->size/float(5)+1; 
	double xiCounter; //counter for xi vector. 

	// i loops of size of U,k,ep,v2,f. 
	#pragma omp parallel num_threads(4) 
	{
	#pragma omp for private(xiCounter)
	for(i = 1; i<size-1; i++)
	{
		double firstTerm,secondTerm,thirdTerm,fourthTerm,val; 
		xiCounter=5*(i-1); //xiCounter is the counter for xi. 
		firstTerm = -(gsl_vector_get(xi,xiCounter+3)-gsl_vector_get(params->XiN,xiCounter+3))/params->deltaT;	
		secondTerm = gsl_vector_get(xi,xiCounter+1)*gsl_vector_get(xi,xiCounter+4) - gsl_vector_get(xi,xiCounter+2)*( ( gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1)));
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Deriv2(xi,0,xiCounter+3,params->grid);
		fourthTerm = Deriv1(xi,0,xiCounter+3,params->grid)*Deriv1vT(vT,i,params->grid);
		val = firstTerm + secondTerm + thirdTerm + fourthTerm;
		Log(logDEBUG3) << "V2 term = " << val<< " at " << i;
		//if (!isfinite(val))
		//	return 1; 
		gsl_vector_set(sysF,xiCounter+3,val); 
	}
	}

	double val; 
	double firstTerm,secondTerm,thirdTerm,fourthTerm;  //as defined in doc. 
	// compute boundary terms. 
	i=size-1; 
	xiCounter=5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+3)-gsl_vector_get(params->XiN,xiCounter+3))/params->deltaT;	
	secondTerm = gsl_vector_get(xi,xiCounter+1)*gsl_vector_get(xi,xiCounter+4) - gsl_vector_get(xi,xiCounter+2)*( (gsl_vector_get(xi,xiCounter+3)/gsl_vector_get(xi,xiCounter+1)));
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDeriv2(xi,xiCounter+3,params->grid);
	val = firstTerm + secondTerm + thirdTerm;
	Log(logDEBUG3) << "V2 term = " << val<< " at " << i;
	if (!isfinite(val))
		return 1; 
	gsl_vector_set(sysF,xiCounter+3,val); 

	return 0; 
}


int SetEpTerms(gsl_vector * xi, gsl_vector * vT,gsl_vector * T,FParams * params, gsl_vector * sysF)
{
	Log(logDEBUG2) << "Setting Ep terms";
	unsigned int i; 
	unsigned int size = vT->size;
	double xiCounter; 
	double ep0 = ComputeEp0(xi,params->modelConst,params->grid);

	//same loop as above. 
	#pragma omp parallel num_threads(4) 
	{
	#pragma omp for private(xiCounter)
	for (i = 1; i<size-1;i++)
	{
		double firstTerm,secondTerm,thirdTerm,fourthTerm,val; 
		xiCounter=5*(i-1); 
		firstTerm = -(gsl_vector_get(xi,xiCounter+2)-gsl_vector_get(params->XiN,xiCounter+2))/params->deltaT;	
		secondTerm = (params->modelConst->Cep1*ComputeP(xi,vT,params->grid,i) - params->modelConst->Cep2*gsl_vector_get(xi,xiCounter+2))/gsl_vector_get(T,i);
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*Deriv2(xi,ep0,xiCounter+2,params->grid);
		fourthTerm = (1/params->modelConst->sigmaEp)*Deriv1(xi,ep0,xiCounter+2,params->grid)*Deriv1vT(vT,i,params->grid);
		val = firstTerm + secondTerm + thirdTerm + fourthTerm;
		Log(logDEBUG3) << "Ep term = " << val << " at " << i;
		//if (!isfinite(val))
		//	return 1; 
		gsl_vector_set(sysF,xiCounter+2,val); 
	}
	}

	double val; 
	double firstTerm, secondTerm,thirdTerm,fourthTerm;  //as in doc. 
	i=size-1;  
	xiCounter=5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+2)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = - (params->modelConst->Cep2*gsl_vector_get(xi,xiCounter+2))/gsl_vector_get(T,i); 
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i)/params->modelConst->sigmaEp)*BdryDeriv2(xi,xiCounter+2,params->grid);
	val = firstTerm + secondTerm + thirdTerm;
	Log(logDEBUG3) << "Ep term = " << val << " at " << i;
	if (!isfinite(val))
		return 1; 
	gsl_vector_set(sysF,xiCounter+2,val);
	return 0; 
}


int SetKTerms(gsl_vector * xi, gsl_vector* vT,FParams * params,gsl_vector *sysF)
{
	Log(logDEBUG2) <<"Setting K terms";
	unsigned int i;  
	unsigned int size=vT->size; 
	double xiCounter; 
	//same loops as above. 
	#pragma omp parallel num_threads(4) 
	{
	#pragma omp for private(xiCounter)
	for(i=1; i<size-1;i++)
	{
		double firstTerm,secondTerm,thirdTerm,fourthTerm,val; 
		xiCounter=5*(i-1);
		firstTerm = -(gsl_vector_get(xi,xiCounter+1)-gsl_vector_get(params->XiN,xiCounter+1))/params->deltaT;	
		secondTerm = ComputeP(xi,vT,params->grid,i)-gsl_vector_get(xi,xiCounter+2);
		thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Deriv2(xi,0,xiCounter+1,params->grid);
		fourthTerm = Deriv1(xi,0,xiCounter+1,params->grid)*Deriv1vT(vT,i,params->grid);
		val = firstTerm + secondTerm + thirdTerm + fourthTerm; 
		Log(logDEBUG3) << "K term = " << val<< " at "<<i;
		//if(!isfinite(val))
		//	return 1; 
		gsl_vector_set(sysF,xiCounter+1,val); 
	}
	}

	double val; 
	double firstTerm, secondTerm, thirdTerm;
	i = size-1;  
	xiCounter = 5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter+1)-gsl_vector_get(params->XiN,xiCounter+1))/params->deltaT;	
	secondTerm = -gsl_vector_get(xi,xiCounter+2); 
	thirdTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDeriv2(xi,xiCounter+1,params->grid);
	val = firstTerm+secondTerm+thirdTerm; 
	Log(logDEBUG3) << "K term = " << val<< " at "<<i;
	if(!isfinite(val))
		return 1; 
	gsl_vector_set(sysF,xiCounter+1,val); 
	return 0; 

}

int SetUTerms( gsl_vector * xi, gsl_vector * vT, FParams * params,gsl_vector * sysF)
{
	Log(logDEBUG2) << "Setting U terms";
	//same structure as other functions. 
	unsigned int i; 
	unsigned int size=vT->size; 
	double xiCounter;

	#pragma omp parallel num_threads(4) 
	{
	#pragma omp for private(xiCounter)
	for(i=1; i<size-1; i++)
	{
		xiCounter = 5*(i-1);
		//cout << omp_get_thread_num() << endl; 

		double firstTerm, secondTerm, thirdTerm,val;
		firstTerm = -(gsl_vector_get(xi,xiCounter)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
		secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*Deriv2(xi,0,xiCounter, params->grid);
		thirdTerm = Deriv1(xi,0,xiCounter, params->grid)*Deriv1vT(vT,i, params->grid);
		val = firstTerm + secondTerm + thirdTerm+1;  
		Log(logDEBUG3) << "U term = " << val << " at " << i;
		//if(!isfinite(val))
		//	return 1; 
		gsl_vector_set(sysF,xiCounter,val);
	}
	}

	double val; 
	double firstTerm, secondTerm, thirdTerm;
	i =size-1; 
	xiCounter = 5*(i-1); 
	firstTerm = -(gsl_vector_get(xi,xiCounter)-gsl_vector_get(params->XiN,xiCounter))/params->deltaT;	
	secondTerm = (1/params->modelConst->reyn + gsl_vector_get(vT,i))*BdryDeriv2(xi,xiCounter, params->grid);
	val = firstTerm+secondTerm + 1; 
	Log(logDEBUG3) << "U term = " << val << " at " << i;
	//if(!isfinite(val))
	//	return 1; 
	gsl_vector_set(sysF,xiCounter,val); 	
	return 0; 
}


