//--------------------------------------------------
// computeTerms: Compute terms T,L,P,vT and f(0),ep(0).
// 
// 12/3/2016 - (gry88) Writen for final project CSE380.  
//-------------------------------------------------- 
#include<gsl/gsl_vector.h>
#include<math.h>
#include"setup.h"
#include"finiteDiff.h"
using namespace std;

double ComputeT(gsl_vector * xi, constants * modelConst,int i)
{
	double firstTerm,secondTerm; //1st and 2nd term as in documentation. 
	double xiCounter = 5*(i-1);  //counter relative to xi. 	
	
	log(verbose,3,"   Computing T\n");	
	firstTerm = (gsl_vector_get(xi,xiCounter+1)/gsl_vector_get(xi,xiCounter+2));
	log(verbose,3,"   Term1 = "+num2st(firstTerm,3) + " at " + num2st(i,3) + "\n");
	if (!isfinite(firstTerm))
	{
		cerr << "Error: T non-finite (" << firstTerm << ")" << endl; 
		cerr << "Note ep = " << gsl_vector_get(xi,xiCounter+2) << " at " << i  << endl; 
		exit(1);
	}

	secondTerm = 6*sqrt(1/(modelConst->reyn*gsl_vector_get(xi,xiCounter+2)));
	log(verbose,3,"   Term2 = " + num2st(secondTerm,3) + " at " + num2st(i,3) + "\n");
	if(!isfinite(secondTerm))
	{
		cerr << "Error: T non-finite (" << secondTerm << ")" << endl;
		cerr << "Note ep = " << gsl_vector_get(xi,xiCounter+2) << " at " << i  << endl; 
		exit(1);
	}

	if (firstTerm >= secondTerm)
			return firstTerm;
	else
			return secondTerm;
}


double ComputeL(gsl_vector * xi,constants * modelConst,int i)
{
	double firstTerm,secondTerm; //see doc.  
	double xiCounter = 5*(i-1); //counter relative to xi.  

	log(verbose,3,"   Computing L\n");
	firstTerm = pow(gsl_vector_get(xi,xiCounter+1),1.5)/gsl_vector_get(xi,xiCounter+2);
	log(verbose,3,"   Term1 = " + num2st(firstTerm,3)+ " at " + num2st(i,3) +"\n");
	if (!isfinite(firstTerm))
		cerr << "Error: L non-finite (" << firstTerm << ")" << endl;
		
	secondTerm = modelConst->Ceta*pow(1/(pow(modelConst->reyn,3)*gsl_vector_get(xi,xiCounter+2)),0.25);
	log(verbose,3,"   Term2 = " + num2st(firstTerm,3)+ " at " + num2st(i,3) + "\n");
	if (!isfinite(secondTerm))
		cerr << "Error: L non-finite (" << secondTerm << ")" << endl;

	if (firstTerm >= secondTerm)
		return modelConst->CL*firstTerm;
	else
		return modelConst->CL*secondTerm;

}

double ComputeEddyVisc(gsl_vector * xi, constants * modelConst,int i)
{
	double val; 
	double xiCounter = 5*(i-1); //counter relative to xi. -1 since U starts a 0. 

	log(verbose,3,"   Computing Eddy Viscosity\n");
	val = modelConst->Cmu*gsl_vector_get(xi,xiCounter+3)*ComputeT(xi,modelConst,i);
	log(verbose,3,"   vT = " + num2st(val,3) + " at " + num2st(i,3) + " at " + num2st(i,3) + "\n"); 
	if (!isfinite(val))
		cerr << "Error: vT non-finite (" << val << ")" << endl;
	return val; 
}

double ComputeP(gsl_vector * xi,constants * modelConst,double deltaEta,int i)
{
	double val; 

	log(verbose,3,"   Computing P\n");
	//We shouldn't need P at 0, but check just incase. 
	if (i==0)
	{
		cerr << "Error: i = 0" << endl; 
		return -1; 
	}
	val = ComputeEddyVisc(xi,modelConst,i)*pow(Diff1(xi,deltaEta,0,i),2);
	log(verbose,3,"   P = " + num2st(val,3) + " at " + num2st(i,3) + "\n");  
	if (!isfinite(val))
		cerr << "Error: P non-finite (" << val << ")" << endl;
	return val; 
}

double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) 
{
	log(verbose,3,"   Compute dissipation at wall boundary\n"); 
	double ep0 = ((2*gsl_vector_get(xi,1))/(modelConst->reyn*pow(deltaEta,2)));
	log(verbose,3,"   Ep(0) = " + num2st(ep0,3) + "\n");
	if (!isfinite(ep0) || ep0 < 0)
		cerr << "Error: unacceptable ep0 (" << ep0 << ")" << endl; 
	return ep0; 
}

double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta)
{
	log(verbose,3,"   Compute f at wall boundary\n"); 
	double f0  = -(( (20*gsl_vector_get(xi,3))/( pow(modelConst->reyn,3)*ComputeEp0(xi,modelConst,deltaEta)*pow(deltaEta,4))));
	log(verbose,3,"   f(0) = " + num2st(f0,3) + "\n");
	if(!isfinite(f0))
		cerr << "Error: f0 non-finite (" << f0 << ")" << endl; 
	return f0; 
}
