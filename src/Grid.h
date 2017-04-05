/*
 * Grid.h
 *
 *  Created on: Apr 5, 2017
 *      Author: clarkp
 */

#ifndef SRC_GRID_H_
#define SRC_GRID_H_

class Grid {
 private:
   const double remap_param;  /// A parameter used to scale the mapping func.
   const double a;            /// Equal to the denominator in the mapping func.
 public:
  const bool isUniform; /// True if the grid is uniform

  Grid(bool isUniform);
  /**
   * \brief
  * @param xi - A coordinate on a uniform grid, with +/-1=wall and 0=center
  * @return The coordinate on a nonuniform grid
  */
  double remap(double xi) const;

  /**
   * \brief Gives the 1st derivative of a uniform grid w.r.t a nonuniform grid
   * @param xi - A point on the uniformly spaced grid.
   * @return The first order derivative of xi w.r.t. y
   */
  double dXidY(double xi) const;

  /**
   * \brief Gives the 2nd derivative of a uniform grid w.r.t a nonuniform grid
   * @param xi - A point on the uniformly spaced grid.
   * @return The second order derivative of xi w.r.t. y
   */
  double d2XidY2(double xi) const;

  int getNumPoints(int Re_tau) const;
};

#endif /* SRC_GRID_H_ */
