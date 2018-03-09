// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2017 Pierre MOULON.
// Copyright (c) 2017 Pascal Monasse

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_POLY_H
#define OPENMVG_NUMERIC_POLY_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

namespace openMVG
{

/**
* @brief Solve roots of a cubic polynomial
* \f$ x^3 + a x^2 + b x + c = 0 \f$
* @param a Coefficient of quadratic parameter
* @param b Coefficient of linear parameter
* @param c Coefficient of scalar parameter
* @param[out] x Found solution(s)
* @return Number of solution
*/
template<typename Real>
int SolveCubicPolynomial( Real a, Real b, Real c,
                          Real x[3] )
{
  const Real eps = std::numeric_limits<Real>::epsilon();
  a /= 3;
  Real p = (b - 3 * a * a) / 3;
  Real q = (2 * a * a * a - a * b + c) / 2;
  Real d = q * q + p * p * p;
  Real tolq = std::max(std::abs(2 * a * a * a),
                       std::max(std::abs(a * b), std::abs(c)));
  Real tolp = std::max(std::abs(b), std::abs(3 * a * a));
  int n = (d > eps * std::max(p * p * tolp, std::abs(q) * tolq)? 1 : 3);
  if (n == 1) // Single root: Cardano's formula
  {
    d = std::pow(std::abs(q) + std::sqrt(d), 1 / (Real)3);
    x[0] = d - p / d;
    if (q > 0)
      x[0] = -x[0];
  }
  else // Three roots: Viete's formula
  {
    if (3 * p >= -eps * tolp) // p=0 and d<=0 implies q=0: triple root 0
    {
      n = 1;
      x[0] = 0;
    }
    else
    {
      p = std::sqrt(-p);
      q /= p * p * p;
      d = Real((q <= -1)? M_PI : (q >= 1)? 0 : std::acos(q));
      for (int i = 0; i < 3; ++i)
        x[i] = Real(-2 * p * std::cos((d + 2 * M_PI * i) / 3));
    }
  }
  for (int i = 0; i < n; ++i)
    x[i] -= a;
  return n;
}


/**
* @brief Solve roots of a cubic polynomial
* @param coeffs Coefficients of the polynomial
* @param[out] solutions Solutions of the polynomial
* @return The number of solutions
*
* @note Input coefficients are in ascending order ( coeffs[N] * x^N )
* @note Assuming coeffs and solutions vectors have 4 values
*/
template<typename Real>
int SolveCubicPolynomial
(
  const Real *coeffs, Real *solutions
)
{
  if ( coeffs[0] == 0.0 )
  {
    return 0;
  }
  const Real a = coeffs[2] / coeffs[3];
  const Real b = coeffs[1] / coeffs[3];
  const Real c = coeffs[0] / coeffs[3];
  return SolveCubicPolynomial( a, b, c, solutions );
}

static std::complex<double> complex_cbrt
(
  const std::complex<double> & z
)
{
  return pow(z, 1. / 3.);
}

/**
* @brief Solve roots of a quartic polynomial.
* @param coeffs Coefficients of the polynomial
* @param[out] real_roots Found roots
*
* @note Input coefficients are in ascending order ( coeffs[N] * x^N )
*/
// Adapted from:
// https://github.com/sidneycadot/quartic/blob/master/solve-quartic.cc
static void solveQuarticPolynomial
(
  const std::array<double, 5> & coeffs,
  std::array<double, 4> & real_roots
)
{

  const double a = coeffs[0];
  const double b = coeffs[1] / a;
  const double c = coeffs[2] / a;
  const double d = coeffs[3] / a;
  const double e = coeffs[4] / a;

  const std::complex<double> Q1 = c * c - 3. * b * d + 12. * e;
  const std::complex<double> Q2 = 2. * c * c * c - 9. * b * c * d
                                  + 27. * d * d + 27. * b * b * e - 72. * c * e;
  const std::complex<double> Q3 = 8. * b * c - 16. * d - 2. * b * b * b;
  const std::complex<double> Q4 = 3. * b * b - 8. * c;

  const std::complex<double> Q5 = complex_cbrt(Q2 / 2.
                                  + sqrt(Q2 * Q2 / 4. - Q1 * Q1 * Q1));
  const std::complex<double> Q6 = (Q1 / Q5 + Q5) / 3.;
  const std::complex<double> Q7 = 2. * sqrt(Q4 / 12. + Q6);

  real_roots = {
    {(-b - Q7 - sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)).real() / 4.,
      (-b - Q7 + sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)).real() / 4.,
      (-b + Q7 - sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)).real() / 4.,
      (-b + Q7 + sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)).real() / 4.}};
}

/// Refine the quartic roots
static void polishQuarticPolynomialRoots
(
  const std::array<double, 5> & coeffs,
  std::array<double, 4> & roots,
  const int iterations = 2)
{
  for (int i = 0; i < iterations; ++i)
  {
    for (auto & root : roots)
    {
      const double error =
        coeffs[4] + root * (coeffs[3] +
                            root * (coeffs[2] +
                                    root * (coeffs[1] +
                                            root * coeffs[0])));

      const double derivative =
        coeffs[3] + root * (2 * coeffs[2] +
                            root * ((4 * coeffs[0] * root + 3 * coeffs[1])));

      root -= error / derivative;
    }
  }
}

}  // namespace openMVG
#endif  // OPENMVG_NUMERIC_POLY_H
