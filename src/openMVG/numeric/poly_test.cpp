
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

TEST(Poly, QuarticPolishPreservesRepeatedRoots) {
  const std::array<std::array<double, 5>, 3> polynomials = {{
    {{1, -6, 13, -12, 4}},  // (x - 1)^2 (x - 2)^2
    {{1, -7, 17, -17, 6}},  // (x - 1)^2 (x - 2) (x - 3)
    {{1, -3, 2, 0, 0}}     // x^2 (x - 1) (x - 2)
  }};
  for (const auto & polynomial : polynomials) {
    std::array<double, 4> roots;
    solveQuarticPolynomial(polynomial, roots);
    const auto initial = roots;
    for (const auto root : initial)
      EXPECT_TRUE(std::isfinite(root));
    polishQuarticPolynomialRoots(polynomial, roots);
    for (size_t i = 0; i < roots.size(); ++i) {
      EXPECT_TRUE(std::isfinite(roots[i]));
      EXPECT_NEAR(initial[i], roots[i], 1e-12);
      const double x = roots[i];
      const double residual = polynomial[4] + x * (polynomial[3] +
        x * (polynomial[2] + x * (polynomial[1] + x * polynomial[0])));
      EXPECT_NEAR(0.0, residual, 1e-12);
    }
  }
}

TEST(Poly, QuarticPolishRefinesSimpleRoots) {
  const std::array<double, 5> polynomial = {{1, -10, 35, -50, 24}};
  std::array<double, 4> roots = {{0.99, 2.01, 2.99, 4.01}};
  const auto initial = roots;
  polishQuarticPolynomialRoots(polynomial, roots, 0);
  for (size_t i = 0; i < roots.size(); ++i)
    EXPECT_EQ(initial[i], roots[i]);
  polishQuarticPolynomialRoots(polynomial, roots, 4);
  for (size_t i = 0; i < roots.size(); ++i)
    EXPECT_NEAR(static_cast<double>(i + 1), roots[i], 1e-12);
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
