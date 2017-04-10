/*
 * Grid.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: clarkp
 */

#include "Grid.h"

#define _USE_MATH_DEFINES // Needed for M_PI
#include <cmath>

Grid::Grid(bool isUniform, double delta, double delta_v)
    : remap_param(0.97), a(remap_param*M_PI/2),
      b(std::sin(remap_param*M_PI/2)), isUniform(isUniform) {
  if (isUniform) {
    size = std::ceil(delta/delta_v);
  } else {
    size = std::ceil(a/(std::asin(b*delta_v/delta-b)+a));
  }

  chi = gsl_vector_alloc(size);
  y = gsl_vector_alloc(size);

  for (unsigned int i = 0; i<size; i++) {
    gsl_vector_set(chi, i, delta*(i+1)/size);
    if (isUniform) {
      gsl_vector_set(y, i, gsl_vector_get(chi, i));
    } else {
      gsl_vector_set(y, i, remap(gsl_vector_get(chi, i)));
    }
  }
}

Grid::~Grid() {
  gsl_vector_free(chi);
  gsl_vector_free(y);
}

double Grid::remap(double chi) const {
  if (isUniform)
    return chi;
  else
    return std::sin((chi-1)*a)/b+1;
}

double Grid::dChidY(double chi) const {
 if (isUniform)
   return 1.0;
 else
   return b/(a*std::sqrt(1.0 - b*b*std::pow(remap(chi)-1,2)));
}

double Grid::d2ChidY2(double chi) const {
 if (isUniform)
   return 0.0;
 else
   return (b*b*(b*remap(chi)-b))/
       (a*pow(1.0-pow(b-b*remap(chi),2),1.5));
}
