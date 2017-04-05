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
   const double remap_param;
   const double a;
 public:
  const bool isUniform;

  Grid(bool isUniform);
  /**
  * @param xi - A coordinate on a uniform grid, with +/-1=wall and 0=center
  * @return The coordinate on a nonuniform grid
  */
  double remap(double xi) const;
  double dXidY(double xi) const;
  double d2XidY2(double xi) const;
};

#endif /* SRC_GRID_H_ */
