
// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/numeric.h"
#include "openMVG/numeric/poly.h"

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

using namespace openMVG;

// Find the polynomial coefficients of x in the equation
//
//   (x - a)(x - b)(x - c) == 0
//
// by expanding to
//
//   x^3 - (c+b+a) * x^2 + (a*b+(b+a)*c) * x - a*b*c = 0.
//           = p               = q              = r
void CoeffsForCubicZeros(std::array<double,3> coeff, // a,b,c
                         std::array<double,3> & res) {
  const double & a = coeff[0];
  const double & b = coeff[1];
  const double & c = coeff[2];
  res[0] = -(c + b + a);
  res[1] = (a * b + (b + a) * c);
  res[2] = -a * b * c;
}

TEST(Poly, SolveCubicPolynomial) {
  std::array<double,3> coeff;
  std::array<double,3> poly_coeff;
  double found_roots[3];

  const double &p = poly_coeff[0];
  const double &q = poly_coeff[1];
  const double &r = poly_coeff[2];

  const double epsilon = 1e-10; 

  coeff = {1,2,3};
  CoeffsForCubicZeros(coeff, poly_coeff);
  CHECK_EQUAL(3, SolveCubicPolynomial(p, q, r, found_roots));
  EXPECT_NEAR(coeff[0], found_roots[0], epsilon);
  EXPECT_NEAR(coeff[1], found_roots[2], epsilon);
  EXPECT_NEAR(coeff[2], found_roots[1], epsilon);

  coeff = {0,1,3};
  CoeffsForCubicZeros(coeff, poly_coeff);
  CHECK_EQUAL(3, SolveCubicPolynomial(p, q, r, found_roots));
  EXPECT_NEAR(coeff[0], found_roots[0], epsilon);
  EXPECT_NEAR(coeff[1], found_roots[2], epsilon);
  EXPECT_NEAR(coeff[2], found_roots[1], epsilon);

  coeff = {-10,0,1};
  CoeffsForCubicZeros(coeff, poly_coeff);
  CHECK_EQUAL(3, SolveCubicPolynomial(p,q,r, found_roots));
  EXPECT_NEAR(coeff[0], found_roots[0], epsilon);
  EXPECT_NEAR(coeff[1], found_roots[2], epsilon);
  EXPECT_NEAR(coeff[2], found_roots[1], epsilon);

  coeff = {-8,1,3};
  CoeffsForCubicZeros(coeff, poly_coeff);
  CHECK_EQUAL(3, SolveCubicPolynomial(p, q, r, found_roots));
  EXPECT_NEAR(coeff[0], found_roots[0], epsilon);
  EXPECT_NEAR(coeff[1], found_roots[2], epsilon);
  EXPECT_NEAR(coeff[2], found_roots[1], epsilon);

  coeff = {28,28,105};
  CoeffsForCubicZeros(coeff, poly_coeff);
  CHECK_EQUAL(3, SolveCubicPolynomial(p, q, r, found_roots));
  EXPECT_NEAR(coeff[0], found_roots[0], epsilon);
  EXPECT_NEAR(coeff[1], found_roots[2], epsilon);
  EXPECT_NEAR(coeff[2], found_roots[1], epsilon);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
