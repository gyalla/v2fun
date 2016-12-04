/**
 * \file 
 * \author Gopal Yalla
 * 
 * \brief Sets up system for gsl's multiroot solvers. 
 * 
 * This file defines the methods to set up the system \f$F(\xi)\f$ of finite difference aproximations
 * as defined in the project proposal.
 */
#ifndef SYSTEMSOLVE_H
#define SYSTEMSOLVE_H
#include<gsl/gsl_vector.h>
#include"setup.h"
using namespace std; 

/** 
 *\brief The struct needed for gsl multiroot solvers. Defined all the parameters
 * needed for the system. 
 */
struct FParams{
	gsl_vector * XiN; /**< vector of unknowns xi at previous time step. */
	double deltaT; /**< step size in time. */
	double deltaEta; /**< step size in wall normal direction. */
	constants * modelConst; /**< pointer to all of the model constants. */
}; 

/**
 * \brief Main function used by gsl_multiroot solver to set up system. 
 *
 * \param xi pointer to gsl_vector of unknowns at n+1 time step. 
 * \param p pointer to parameters for system. 
 * \param sysF gsl_vector defining multiroot function. 
 * \return Error code (0 = success). 
 */
int SysF(const gsl_vector * xi, void * p, gsl_vector * sysF);

/** 
 * \brief Sets terms in system related to mean veloctity, U. 
 * \param xi pointer to gsl_vector of unknowns \f$
 * U,k,\epsilon,\overline{v^2},f\f$. 
 * \param vT pointer to vector of eddy viscosity. 
 * \param params pointer to parameters of system. 
 * \param sysF gsl_vector defining multiroot function. 
 * \return Error code (0 = success). 
 */
int SetUTerms(gsl_vector * xi, gsl_vector * vT,FParams * params,gsl_vector * sysF);

/** 
 * \brief Sets terms in system related to kinetic energy, k. 
 * \param xi pointer to gsl_vector of unknowns \f$
 * U,k,\epsilon,\overline{v^2},f\f$. 
 * \param vT pointer to vector of eddy viscosity. 
 * \param params pointer to parameters of system. 
 * \param sysF gsl_vector defining multiroot function. 
 * \return Error code (0 = success). 
 */
int SetKTerms(gsl_vector * xi,gsl_vector * vT, FParams * params,gsl_vector *sysF);

/** 
 * \brief Sets terms in system related to dissipation, \f$\epsilon\f$. 
 * \param xi pointer to gsl_vector of unknowns \f$
 * U,k,\epsilon,\overline{v^2},f\f$. 
 * \param vT pointer to vector of eddy viscosity. 
 * \param params pointer to parameters of system. 
 * \param sysF gsl_vector defining multiroot function. 
 * \return Error code (0 = success). 
 */
int SetEpTerms(gsl_vector * xi, gsl_vector * vT, FParams * params, gsl_vector * sysF);

/** 
 * \brief Sets terms in system related to velocity scale, \f$\overline{v^2}\f$. 
 * \param xi pointer to gsl_vector of unknowns \f$
 * U,k,\epsilon,\overline{v^2},f\f$. 
 * \param vT pointer to vector of eddy viscosity. 
 * \param params pointer to parameters of system. 
 * \param sysF gsl_vector defining multiroot function. 
 * \return Error code (0 = success). 
 */
int SetV2Terms(gsl_vector * xi, gsl_vector * vT,  FParams * params, gsl_vector * sysF);

/** 
 * \brief Sets terms in system related to mean veloctity, U. 
 * \param xi pointer to gsl_vector of unknowns \f$
 * U,k,\epsilon,\overline{v^2},f\f$. 
 * \param params pointer to parameters of system. 
 * \param sysF gsl_vector defining multiroot function. 
 * \return Error code (0 = success). 
 */
int SetFTerms(gsl_vector * xi, FParams * params, gsl_vector * sysF);

#endif
