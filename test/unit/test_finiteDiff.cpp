/**
 * \file: test_finiteDiff.cpp
 * \brief: Tests whether the nonuniform finite differences are correct.
 */
#include "test_finiteDiff.h"
#include <iostream>
#include <limits>
#include <cmath>

#include "../../src/Grid.h"
#include "../../src/finiteDiff.h"

int test_remapping() {
  // ARRANGE
  Grid grid(false, 1.0, 0.1);

  double xi = 1.0;

  // ACT
  double y = grid.remap(xi);

  // ASSERT
  double tol = 5*std::numeric_limits<double>::epsilon();
  if (std::abs(y - 1.00) > tol) {
    std::cout << "FAIL: Calculated remapping incorrectly" << std::endl;
    return 1;
  }
  std::cout << "PASS: Calculated remapping correctly" << std::endl;
  return 0;
}

int test_first_deriv(int n) {
  // SETUP
  Grid grid(false, 1, 1.0/n);
  gsl_vector* phi = gsl_vector_alloc(5*grid.getSize());

  // ARRANGE
  // phi = y
  for (unsigned int i=0; i<grid.getSize(); i++) {
    gsl_vector_set(phi, i*5,  gsl_vector_get(grid.y, i));
  }

  // ACT
  double deriv[n-1];
  for (unsigned int index=0; index<grid.getSize()-1; index++)
    deriv[index] = Deriv1(phi, 0, index*5, &grid);

  // ASSERT
  double tol = 1e-3;
  for (unsigned int index = 1; index < grid.getSize()-1; index++) {
    if (std::abs(deriv[index] - 1.0) > tol) {
      std::cout << "FAIL: First derivative calculated incorrectly!" << std::endl;
      std::cout << "    Expected: " << 1.0 << "   Found: " << deriv[index] << std::endl;
      std::cout << "    At y = " << gsl_vector_get(grid.y, index) << std::endl;
      std::cout << "    Tolerance: " << tol << std::endl;
      return 1;
    }
  }
  std::cout << "PASS: First derivative on nonuniform grid" << std::endl;
  return 0;
}

int test_second_deriv(int n) {
  // SETUP
  Grid grid(false, 1.0, 1.0/n);
  gsl_vector* phi = gsl_vector_alloc(5*grid.getSize());

  // ARRANGE
  // phi = 0.5*y^2
  for (unsigned int i=0; i<grid.getSize(); i++) {
    gsl_vector_set(phi, i*5,  0.5*std::pow(gsl_vector_get(grid.y, i),2));
  }

  // ACT
  double deriv[grid.getSize()-1];
  for (unsigned int index=0; index<grid.getSize()-1; index++)
    deriv[index] = Deriv2(phi, 0, index*5, &grid);

  // ASSERT
  double tol = 1e-2;
  for (unsigned int index = 1; index < grid.getSize()-1; index++) {
    if (std::abs(deriv[index] - 1.0) > tol) {
      std::cout << "FAIL: Second derivative calculated incorrectly!" << std::endl;
      std::cout << "    Expected: " << 1.0 << "   Found: " << deriv[index] << std::endl;
      std::cout << "    At y = " << gsl_vector_get(grid.y, index) << std::endl;
      std::cout << "    Tolerance: " << tol << std::endl;
      return 1;
    }
  }
  std::cout << "PASS: Second derivative on nonuniform grid" << std::endl;
  return 0;
}
