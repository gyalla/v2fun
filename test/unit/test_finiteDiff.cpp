/**
 * \file: test_finiteDiff.cpp
 * \brief: Tests whether the nonuniform finite differences are correct.
 */
#include "test_finiteDiff.h"
#include <iostream>
#include <limits>
#include <cmath>

#include "../../src/finiteDiff.h"

int test_remapping() {
  // ARRANGE
  double xi = 0.1;

  // ACT
  double y = remap(xi);

  // ASSERT
  double tol = 5*std::numeric_limits<double>::epsilon();
  if (std::abs(y - 0.151947053906880373886) > tol) {
    std::cout << "FAIL: Calculated remapping incorrectly" << std::endl;
    return 1;
  }
  std::cout << "PASS: Calculated remapping correctly" << std::endl;
  return 0;
}

int test_first_deriv(int n) {
  // SETUP
  gsl_vector* phi = gsl_vector_alloc(5*n);

  // ARRANGE
  // phi = y
  for (int i=0; i<n; i++) {
    double xi = -1.0+1.0*i/(n-1);
    gsl_vector_set(phi, i*5,  remap(xi));
  }

  // ACT
  double deriv[n-1];
  for (int index=0; index<n-1; index++)
    deriv[index] = Deriv1(phi, 1.0/(n-1), 0, index*5);

  // ASSERT
  double tol = 1e-3;
  for (int index = 1; index < n-1; index++) {
    if (std::abs(deriv[index] - 1.0) > tol) {
      std::cout << "FAIL: First derivative calculated incorrectly!" << std::endl;
      std::cout << "    Expected: " << 1.0 << "   Found: " << deriv[index] << std::endl;
      std::cout << "    At y = " << remap(-1.0+index*1.0/(n-1)) << std::endl;
      std::cout << "    Tolerance: " << tol << std::endl;
      return 1;
    }
  }
  std::cout << "PASS: First derivative on nonuniform grid" << std::endl;
  return 0;
}

int test_second_deriv(int n) {
  // SETUP
  gsl_vector* phi = gsl_vector_alloc(5*n);

  // ARRANGE
  // phi = 0.5*y^2
  for (int i=0; i<n; i++) {
    double xi = -1.0+1.0*i/(n-1);
    gsl_vector_set(phi, i*5,  0.5*remap(xi)*remap(xi));
  }

  // ACT
  double deriv[n-1];
  for (int index=0; index<n-1; index++)
    deriv[index] = Deriv2(phi, 1.0/(n-1), 0, index*5);

  // ASSERT
  double tol = 1e-2;
  for (int index = 1; index < n-1; index++) {
    if (std::abs(deriv[index] - 1.0) > tol) {
      std::cout << "FAIL: First derivative calculated incorrectly!" << std::endl;
      std::cout << "    Expected: " << 1.0 << "   Found: " << deriv[index] << std::endl;
      std::cout << "    At y = " << remap(-1.0+index*1.0/(n-1)) << std::endl;
      std::cout << "    Tolerance: " << tol << std::endl;
      return 1;
    }
  }
  std::cout << "PASS: First derivative on nonuniform grid" << std::endl;
  return 0;
}
