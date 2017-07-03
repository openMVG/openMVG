// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_POLY_H
#define OPENMVG_NUMERIC_POLY_H

#include <algorithm>
#include <cmath>

namespace openMVG
{

/**
* @brief Solve roots of a cubic polynomial
* \f$ x^3 + a x^2 + b x + c = 0 \f$
* @param a Coefficient of cubic parameter
* @param b Coefficient of linear parameter
* @param c Coefficient of scalar parameter
* @param[out] x0 potential solution
* @param[out] x1 potential solution
* @param[out] x2 potential solution
* @return Number of solution
* @note If number of solution is less than 3, only the first output parameters are computed
* @note The GSL cubic solver was used as a reference for this routine.
*/
template<typename Real>
int SolveCubicPolynomial( Real a, Real b, Real c,
                          Real *x0, Real *x1, Real *x2 )
{
  const Real q = a * a - 3 * b;
  const Real r = 2 * a * a * a - 9 * a * b + 27 * c;

  const Real Q = q / 9;
  const Real R = r / 54;

  const Real Q3 = Q * Q * Q;
  const Real R2 = R * R;

  const Real CR2 = 729 * r * r;
  const Real CQ3 = 2916 * q * q * q;

  if ( R == 0 && Q == 0 )
  {
    // Tripple root in one place.
    *x0 = *x1 = *x2 = -a / 3;
    return 3;

  }
  else if ( CR2 == CQ3 )
  {
    // This test is actually R2 == Q3, written in a form suitable for exact
    // computation with integers.
    //
    // Due to finite precision some double roots may be missed, and considered
    // to be a pair of complex roots z = x +/- epsilon i close to the real
    // axis.
    const Real sqrtQ = sqrt ( Q );
    if ( R > 0 )
    {
      *x0 = -2 * sqrtQ - a / 3;
      *x1 =      sqrtQ - a / 3;
      *x2 =      sqrtQ - a / 3;
    }
    else
    {
      *x0 =     -sqrtQ - a / 3;
      *x1 =     -sqrtQ - a / 3;
      *x2 =  2 * sqrtQ - a / 3;
    }
    return 3;

  }
  else if ( CR2 < CQ3 )
  {
    // This case is equivalent to R2 < Q3.
    const Real sqrtQ = sqrt ( Q );
    const Real sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
    const Real theta = acos ( R / sqrtQ3 );
    const Real norm = -2 * sqrtQ;
    *x0 = norm * cos ( theta / 3 ) - a / 3;
    *x1 = norm * cos ( ( theta + 2.0 * M_PI ) / 3 ) - a / 3;
    *x2 = norm * cos ( ( theta - 2.0 * M_PI ) / 3 ) - a / 3;

    // Put the roots in ascending order.
    if ( *x0 > *x1 )
    {
      std::swap( *x0, *x1 );
    }
    if ( *x1 > *x2 )
    {
      std::swap( *x1, *x2 );
      if ( *x0 > *x1 )
      {
        std::swap( *x0, *x1 );
      }
    }
    return 3;
  }
  const Real sgnR = ( R >= 0 ? 1 : -1 );
  const Real A = -sgnR * pow ( std::abs ( R ) + sqrt ( R2 - Q3 ), 1.0 / 3.0 );
  const Real B = Q / A;
  *x0 = A + B - a / 3;
  return 1;
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
int SolveCubicPolynomial( const Real *coeffs, Real *solutions )
{
  if ( coeffs[0] == 0.0 )
  {
    // TODO(keir): This is a quadratic not a cubic. Implement a quadratic
    // solver!
    return 0;
  }
  const Real a = coeffs[2] / coeffs[3];
  const Real b = coeffs[1] / coeffs[3];
  const Real c = coeffs[0] / coeffs[3];
  return SolveCubicPolynomial( a, b, c,
                               solutions + 0,
                               solutions + 1,
                               solutions + 2 );
}
}  // namespace openMVG
#endif  // OPENMVG_NUMERIC_POLY_H
