/**
 * \file 
 * \author Gopal Yalla
 *
 * \brief Defines finite difference approximations. 
 *
 * This file defines the methods for finite difference approximations with respect to 
 * the structure of \f$\xi\f$ as defined in the model documentation.  
 *
 */
#ifndef FINITEDIFF_H
#define FINITEDIFF_H
#include<gsl/gsl_vector.h>
#include "Grid.h"
using namespace std;



/**
 * \brief Approximates first derivative of eddy viscosity vector using center
 * difference. 
 *
 *
 * \param x pointer to vector to compute finite difference approximation. 
 * \param deltaEta step size in wall normal direction. 
 * \param i point at which to compute first derivative around (not relative to \f$\xi\f$ in this case).  
 * \return centered difference approximation. 
 */
double Diff1vT(gsl_vector * x,double deltaEta,int i);

/**
 * \brief Approximates second derivative of gsl_vector using center
 * difference.
 *
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction.
 * \param bdry value at boundary.
 * \param i point at which to compute second derivative around (relative to ordering of \f$\xi\f$).
 * \return centered difference approximation.
 */
inline double Diff2(gsl_vector * x,double deltaEta,double bdry, int i);

/**
 * \brief Approximates first derivative of gsl_vector using center
 * difference.
 *
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction.
 * \param bdry value at boundary.
 * \param i point at which to compute first derivative around (relative to ordering of \f$\xi\f$).
 * \return centered difference approximation.
 */
inline double Diff1(gsl_vector *x, double deltaEta,double bdry,int i);

/**
 * \brief Approximates second derivative of gsl_vector using center
 * difference.
 *
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param bdry value at boundary.
 * \param i point at which to compute second derivative around (relative to ordering of \f$\xi\f$).
 * \return centered difference approximation.
 */
double Deriv2(gsl_vector * x, double bdry, int i, Grid* grid);

/**
 * \brief Approximates first derivative of gsl_vector using center
 * difference.
 *
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param bdry value at boundary.
 * \param i point at which to compute first derivative around (relative to ordering of \f$\xi\f$).
 * \return centered difference approximation.
 */
double Deriv1(gsl_vector *x, double bdry, int i, Grid* grid);

/**
 * \brief Approximates second derivative of gsl_vector at boundary using centered difference making use of ghost points and zero Nuemann boundary conditions.
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction. 
 * \param i point at which to compute second derivative around (relative to ordering of \f$\xi\f$).  
 * \return centered difference approximation. 
 */
double BdryDeriv2(gsl_vector * x , int i, Grid* grid);

/**
 * \brief Approximates second derivative of gsl_vector at boundary using centered difference making use of ghost points and zero Nuemann boundary conditions.
 *
 * \param x pointer to vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param deltaEta step size in wall normal direction.
 * \param i point at which to compute second derivative around (relative to ordering of \f$\xi\f$).
 * \return centered difference approximation.
 */
double BdryDiff2(gsl_vector * x ,double deltaEta,int i);

/**
 * \brief Approximates first derivative of eddy viscosity vector using center
 * difference.
 *
 *
 * \param x pointer to vector to compute finite difference approximation.
 * \param deltaEta step size in wall normal direction.
 * \param i point at which to compute first derivative around (not relative to \f$\xi\f$ in this case).
 * \return centered difference approximation.
 */
double Deriv1vT(gsl_vector * x, int i, Grid* grid);

#endif
