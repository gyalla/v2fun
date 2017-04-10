/*
 * Grid.h
 *
 *  Created on: Apr 5, 2017
 *      Author: clarkp
 */

#ifndef SRC_GRID_H_
#define SRC_GRID_H_

#include<gsl/gsl_vector.h>

/**
 * \brief Used to define either a uniform or a non-uniform grid.
 */
class Grid {
 private:
   const double remap_param;  /// A parameter used to scale the mapping func.
   const double a;            /// Equal to remap_param * pi/2
   const double b;            /// Equal to the denominator in the mapping func.
   unsigned int size;
 public:
  const bool isUniform; /// True if the grid is uniform
  gsl_vector* chi;
  gsl_vector* y;

  /**
   * \brief Constructor
   * @param isUniform - True if the grid is to be uniform
   * @param delta - The channel half-width.
   * @param delta_v - The viscous length-scale.
   */
  Grid(bool isUniform, double delta, double delta_v);\

  /**
   * \brief Destructor
   */
  ~Grid();

  /**
   * \brief
  * @param chi - A coordinate on a uniform grid, with +/-1=wall and 0=center
  * @return The coordinate on a nonuniform grid
  */
  double remap(double chi) const;

  /**
   * \brief Gives the 1st derivative of a uniform grid w.r.t a nonuniform grid
   * @param \chi - A point on the uniformly spaced grid.
   * @return The first order derivative of xi w.r.t. y
   */
  double dChidY(double chi) const;

  /**
   * \brief Gives the 2nd derivative of a uniform grid w.r.t a nonuniform grid
   * @param chi - A point on the uniformly spaced grid.
   * @return The second order derivative of xi w.r.t. y
   */
  double d2ChidY2(double chi) const;

  unsigned int getSize() const { return size;};
};

#endif /* SRC_GRID_H_ */
