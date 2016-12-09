/**
 * \file 
 * \author Gopal Yalla
 *
 * \brief Computes terms \f$T,L,P,\nu_T,\f$ as well as \f$f(0),\epsilon(0)\f$. 
 *
 * This file defines the methods to compute various terms for the v2-f equations 
 * including turublent time scale, turbulent length scale, production rate, eddy viscosity, 
 * as well as the wall boundary terms for the dissipation and redistribution term. 
 */
#ifndef COMPUTETERMS_H
#define COMPUTETERMS_H

#include<gsl/gsl_vector.h>
#include"setup.h"
using namespace std;
/**
 * \brief Compute turbulent time scale, T.
 *
 * Compute T using \f[T = \max \left\{ \frac{k}{\epsilon}, 6 \sqrt{\frac{1}{Re_\tau \epsilon}} \right\} \f]
 * \param xi pointer to gsl_vector of unknowns \f$U,k,\epsilon,\overline{v^2},f\f$.
 * \param modelConst pointer to struct containing model constants. 
 * \param i position at which to compute T. 
 * \return T at i. 
 */
double ComputeT(gsl_vector * xi, constants * modelConst,int i);
/**
 * \brief Compute turbulent length scale, L. 
 *
 * Compute L using \f[ L = C_L \max \left\{ \frac{k^{3/2}}{\epsilon}, C_\eta \left(
	\frac{1}{Re_\tau^3 \epsilon}
	\right)^{1/4} \right\} \f]
 * \param xi pointer to gsl_vector of unknowns, \f$U,k,\epsilon,\overline{v^2},f\f$.
 * \param modelConst pointer to struct containing model constants.
 * \param i position at which to compute L.
 * \return L at i. 
 */
double ComputeL(gsl_vector * xi, constants * modelConst,int i);

/**
 * \brief Compute eddy viscosity. 
 *
 * Compute eddy viscosity using \f[ \nu_T = C_\mu \overline{v^2} T \f]
 * \param xi pointer to gsl_vector of unknowns \f$U,k,\epsilon,\overline{v^2},f\f$.
 * \param T pointer to gsl_vector of turbulent time scale
 * \param modelConst pointer to struct containing model constants. 
 * \param i position at which to compute \f$\nu_T\f$. 
 * \return \f$\nu_T\f$ at i. 
 */
double ComputeEddyVisc(gsl_vector * xi, gsl_vector * T, constants * modelConst,int i);

/**
 * \brief Comute production rate. 
 *
 * Compute production rate using \f[P=\nu_T \left( \frac{\partial U^+}{\partial \eta} \right)^2\f]
 * \param xi pointer to gsl_vector of unknowns \f$U,k,\epsilon,\overline{v^2},f\f$.
 * \param vT pointer to gsl_vector of eddy viscosity.
 * \param i position at which to compute P. 
 * \return P at i.  
 */
double ComputeP(gsl_vector * xi,gsl_vector * vT,double deltaEta,int i );

/**
 * \brief Compute redistribution term at wall boundary. 
 *
 * Compute f at boundary using \f[ f(0)= -\frac{20\overline{v^2}_1}{Re_\tau^2 \epsilon(0) \Delta \eta^4} \f]
 * \param xi pointer to gsl_vector of unknowns \f$U,k,\epsilon,\overline{v^2},f\f$. 
 * \param modelConst pointer to struct containing model constatns. 
 * \param deltaEta step size in wall normal direction. 
 * \return f at boundary. 
 */
double Computef0(gsl_vector * xi,constants * modelConst,double deltaEta);

/**
 * \brief Compute dissipation term at wall boundary. 
 *
 * Compute \f$\epsilon\f$ at boundary using \f[ \epsilon(0) = 2 \frac{k_1}{Re_\tau \Delta \eta^2}\f]
 * \param xi pointer to gsl_vector of unknowns \f$U,k,\epsilon,\overline{v^2},f\f$. 
 * \param modelConst pointer to struct containing model constants. 
 * \param deltaEta step size in wall normal direction. 
 * \return \f$\epsilon\f$ at boundary. 
 */
double ComputeEp0(gsl_vector * xi,constants * modelConst,double deltaEta) ;
#endif
