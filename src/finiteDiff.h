/**
 * \file 
 * \author Gopal Yalla
 *
 * \brief Defines finite difference approximations. 
 *
 * This file defines the method for finite difference approximations with respect to 
 * the structure of \f$\xi\f$. 
 *
 */
#ifndef FINITEDIFF_H
#define FINITEDIFF_H
#include<gsl/gsl_vector.h>
using namespace std;


/**
 * \brief Approximates first derivative of eddy viscosity vector using second order center
 * difference. 
 *
 *
 * \param x pointer to vector to compute finite difference approximation. 
 * \param deltaEta step size in wall normal direction. 
 * \param i point at which to compute first derivative around.  
 * \return centered difference approximation. 
 */
double Diff1vT(gsl_vector * x,double deltaEta,int i);

/**
 * \brief Approximates second derivative of gsl_vector using second order center
 * difference. 
 *
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction. 
 * \param bdry value at boundary.  
 * \param i point at which to compute second derivative around (relative to ordering of xi).  
 * \return centered difference approximation. 
 */
double Diff2(gsl_vector * x,double deltaEta,double bdry, int i);

/**
 * \brief Approximates first derivative of gsl_vector using second order center
 * difference. 
 *
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction. 
 * \param bdry value at boundary.  
 * \param i point at which to compute first derivative around (relative to ordering of xi).  
 * \return centered difference approximation. 
 */
double Diff1(gsl_vector *x, double deltaEta,double bdry,int i);

/**
 * \brief Approximates second derivative of gsl_vector at boundary using centered difference making use of ghost points and zero Nuemann boundary conditions.
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction. 
 * \param i point at which to compute second derivative around (relative to ordering of xi).  
 * \return centered difference approximation. 
 */
double BdryDiff2(gsl_vector * x ,double deltaEta,int i);

#endif
