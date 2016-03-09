
// Copyright (c) 2010 libmv authors.
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


// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
//
// Five point relative pose computation using Groebner basis.
// We follow the algorithm of [1] and apply some optimization hints of [2].
//
// [1] H. Stewénius, C. Engels and D. Nistér,  "Recent Developments on Direct
//     Relative Orientation",  ISPRS 2006
//
// [2] D. Nistér,  "An Efficient Solution to the Five-Point Relative Pose",
//     PAMI 2004


#ifndef OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_FIVE_POINT_H_
#define OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_FIVE_POINT_H_

#include "openMVG/numeric/numeric.h"
#include <vector>

namespace openMVG
{
using namespace std;

/**
 * @brief Computes the relative pose of two calibrated cameras from 5 correspondences.
 *
 * \param x1 Points in the first image.  One per column.
 * \param x2 Corresponding points in the second image. One per column.
 * \param E  A list of at most 10 candidate essential matrix solutions.
 */
void FivePointsRelativePose( const Mat2X &x1, const Mat2X &x2,
                             vector<Mat3> *E );

/**
* @brief Compute the nullspace of the linear constraints given by the matches.
* @param x1 Match position in first camera
* @param x2 Match position in second camera
* @return Nullspace (homography) that maps x1 points to x2 points
*/
Mat FivePointsNullspaceBasis( const Mat2X &x1, const Mat2X &x2 );

/**
* @brief Multiply two polynomials of degree 1.
* @param a Polynomial a1 + a2 x + a3 y + a4 z
* @param b Polynomial b1 + b2 x + b3 y + b4 z
* @return Product of a and b :
* res = a1 b1 +
        (a1b2 + b1a2) x +
        (a1b3 + b1a3) y +
        (a1b4 + b1a4) z +
        (a2b3 + b2a3) xy +
        (a2b4 + b2a4) xz +
        (a3b4 + b3a4) yz +
        a2b2 x^2 +
        a3b3 y^2 +
        a4b4 z^2
* @note Ordering is defined as follow :
* [xxx xxy xyy yyy xxz xyz yyz xzz yzz zzz xx xy yy xz yz zz x y z 1]
*/
Vec o1( const Vec &a, const Vec &b );

/**
* @brief Multiply two polynomials of degree 2
* @param a Polynomial a1 + a2 x + a3 y + a4 z + a5 x^2 + a6 y^2 + a7 z^2
* @param b Polynomial b1 + b2 x + b3 y + b4 z + b5 x^2 + b6 y^2 + b7 z^2
* Product of a and b
* @note Ordering is defined as follow :
* [xxx xxy xyy yyy xxz xyz yyz xzz yzz zzz xx xy yy xz yz zz x y z 1]
*/
Vec o2( const Vec &a, const Vec &b );

/**
* Builds the polynomial constraint matrix M.
* @param E_basis Basis essential matrix
* @return polynomial constraint associated to the essential matrix
*/
Mat FivePointsPolynomialConstraints( const Mat &E_basis );

// In the following code, polynomials are expressed as vectors containing
// their coeficients in the basis of monomials:
//
//  [xxx xxy xyy yyy xxz xyz yyz xzz yzz zzz xx xy yy xz yz zz x y z 1]
//
// Note that there is an error in Stewenius' paper.  In equation (9) they
// propose to use the basis:
//
//  [xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz xx xy xz yy yz zz x y z 1]
//
// But this is not the basis used in the rest of the paper, neither in
// the code they provide.  I (pau) have spend 4 hours debugging and
// reverse engineering their code to find the problem. :(
enum
{
  coef_xxx,
  coef_xxy,
  coef_xyy,
  coef_yyy,
  coef_xxz,
  coef_xyz,
  coef_yyz,
  coef_xzz,
  coef_yzz,
  coef_zzz,
  coef_xx,
  coef_xy,
  coef_yy,
  coef_xz,
  coef_yz,
  coef_zz,
  coef_x,
  coef_y,
  coef_z,
  coef_1
};

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_FIVE_POINT_H_

