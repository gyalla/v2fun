/*
 * Grid.h
 *
 *  Created on: Apr 5, 2017
 *      Author: clarkp
 */

#ifndef SRC_GRID_H_
#define SRC_GRID_H_

#include<gsl/gsl_vector.h>
#include <cmath>

/**
 * \brief Used to define either a uniform or a non-uniform grid.
 */
class Grid {
 private:
   const double remap_param;  /// A parameter used to scale the mapping func.
   const double a;            /// Equal to remap_param * pi/2
   const double b;            /// Equal to the denominator in the mapping func.
   unsigned int size;
   mutable double cached_deriv1;
   mutable double cached_deriv2;
 public:
  const bool isUniform; /// True if the grid is uniform
  mutable bool isDeriv1Cached;
  mutable bool isDeriv2Cached;
  gsl_vector* chi;
  gsl_vector* y;

  /**
   * \brief Constructor
   * @param isUniform - True if the grid is to be uniform
   * @param delta - The channel half-width.
   * @param delta_v - The viscous length-scale.
   */
  Grid(bool isUniform, double delta, double delta_v);

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
  inline double dChidY(double chi) const;

  /**
   * \brief Gives the 2nd derivative of a uniform grid w.r.t a nonuniform grid
   * @param chi - A point on the uniformly spaced grid.
   * @return The second order derivative of xi w.r.t. y
   */
  inline double d2ChidY2(double chi) const;

  unsigned int getSize() const { return size;};
};

inline double Grid::dChidY(double chi) const {
  if (!isDeriv1Cached) {
    if (isUniform)
      cached_deriv1 = 1.0;
    else
      cached_deriv1 = b/(a*std::sqrt(1.0 - b*b*std::pow(remap(chi)-1,2)));
    isDeriv1Cached = true;
  }
  return cached_deriv1;
}

inline double Grid::d2ChidY2(double chi) const {
  if (!isDeriv2Cached) {
    if (isUniform)
      cached_deriv2 = 0.0;
    else
      cached_deriv2 = (b*b*(b*remap(chi)-b))/
                      (a*pow(1.0-pow(b-b*remap(chi),2),1.5));
    isDeriv2Cached = true;
  }
  return cached_deriv2;
}

#endif /* SRC_GRID_H_ */
