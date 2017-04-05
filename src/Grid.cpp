/*
 * Grid.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: clarkp
 */

#include "Grid.h"

#define _USE_MATH_DEFINES // Needed for M_PI
#include <cmath>

Grid::Grid(bool isUniform)
    : remap_param(0.97), a(std::sin(remap_param*M_PI/2)), isUniform(isUniform){

}

double Grid::remap(double xi) const {
  return sin(remap_param*xi*M_PI/2)/sin(remap_param*M_PI/2);
}

double Grid::dXidY(double xi) const {
 return 2.0*a/(remap_param*M_PI*sqrt(1.0 - pow(a*remap(xi),2)));
}

double Grid::d2XidY2(double xi) const {
 return (2.0*a*a*a*remap(xi))/
     (remap_param*M_PI*pow(1.0-pow(a*remap(xi),2),1.5));
}
