/**
 * \file 
 * \author Gopal Yalla
 *
 * \brief Initializes terms for v2-f simulation and provides some basic routines.  
 *
 * This file defines the methods used by v2fun, i.e., to get inputs, 
 * solve for the initial conditions, save output, and defines a struct of model constants. 
 */
#ifndef OPTIONS_H
#define OPTIONS_H

#include<gsl/gsl_vector.h>
#include<iostream>
#include"../include/loglevel.h"
#include "Grid.h"

using namespace std; 
/**
 * \brief Holds all of the model constants. 
 */
struct constants {
	double reyn; /**< Reynolds number. */
	double Cmu; /**< \f$C_\mu\f$ */
	double C1; /**< \f$C_1\f$ */
	double C2; /**< \f$C_2\f$ */
	double Cep1; /**< \f$C_{\epsilon 1}\f$ */
	double Cep2; /**< \f$C_{\epsilon 2}\f$ */
	double Ceta; /**< \f$C_{\eta}\f$ */
	double CL; /**< \f$C_L\f$ */
	double sigmaEp; /**< \f$\sigma_\epsilon\f$ */
};

/**
 * \brief Parse inputs. 
 *
 * Uses the Boost Library to parse inputs from file.
 * \param modelConst pointer to struct containing model constants. 
 * \param filename Defines data file to use for initial conditions. 
 * \param outFile Defines file to write results to.
 * \param uniformGrid If true, the grid will be uniform.
 * \param max_ts Defines maximum time steps taken before exiting..
 * \param restarting If true, the simulation is picking up where left off. Data file contains f.
 * \return Error code (0 = success).
 */
int Input_Parse(constants * modelConst,string &filename,string & outFile,
                     bool &uniformGrid, int &max_ts, bool &restarting,int ac,char ** av);

/**
 * \brief Solve for initial conditions.  
 * 
 * Solves for initial conditions of \f$ U,k,\epsilon,\overline{v^2}\f$ using linear 
 * interpolation from Moser channel flow code. 
 * \param xi pointer to vector of \f$ U,k,\epsilon,\overline{v^2},f\f$.
 * \param grid - A point to the grid of points
 * \param file filename to get data from. 
 * \param restarting If true, the simulation is picking up where left off. Data file contains f.
 * \return Error code (0 = success). 
 */


int SolveIC(gsl_vector* xi, constants * modelConst, Grid* grid, string file,bool restarting);

/**
 * \brief Linear interpolates inputs and places result in vector. 
 *
 * Does linear interpolation using, \f[ xi(i) = y_1 + (x - x_0)* \left( \frac{ (y_1 - y_0)}{x_1-x_0} \right) \f]
 * \param Vec xi in the formula given. 
 * \param pt1 \f$x_1\f$ in the formula given. 
 * \param pt2 \f$x_0\f$ in the formula given. 
 * \param U1  \f$y_1\f$ in the formula given. 
 * \param U2 \f$y_0\f$ in the formula given. 
 * \param gridPt \f$x\f$ in the formula given. 
 * \param i position in xi to place result. 
 * \param restarting If true, the simulation is picking up where left off. Data file contains f.
 * \param modelConstants pointer to structcontaining model constants. 
 * \return Error code (0 = success).
 */
int LinInterp(gsl_vector * Vec,double pt1,double pt2, double U1,double U2,double gridPt,int i, constants * modelConst, bool restarting );

/**
 * \brief Solve for initical conditions for f. 
 * 
 * Since Moser/MK data does not provide information for redistribution term, f, 
 * we solve for it using a finite difference formula. 
 * \param xi pointer to gsl_vector of unknowns \f$ U,k,\epsilon,\overline{v^2},f\f$. 
 * \param modelConst pointer to struct containing model constants.
 * \param grid - A point to the grid of points
 * \return Error code (0 = success). 
 */
int Solve4f0(gsl_vector * xi,constants * modelConst, Grid* grid);

/**
 * \brief Writes result to file
 *
 * \param xi pointer to gsl_vector of unknowns. 
 * \param filename name of output file.
 * \param grid - A point to the grid of points
 * \param modelConst pointer to struct of model constants. 
 */
void SaveResults(gsl_vector *xi,string filename, Grid* grid, constants * modelConst);

#endif

