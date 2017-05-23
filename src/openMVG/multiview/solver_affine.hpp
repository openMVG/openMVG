
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

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SOLVER_MULTIVIEW_AFFINE_HPP
#define OPENMVG_SOLVER_MULTIVIEW_AFFINE_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

/** 2D Affine transformation estimation
 *
 * This function can be used in order to estimate the affine transformation
 * between two sets of points with known 2D correspondences.
 *
 * \param[in] x1 The first 2xN matrix of euclidean points
 * \param[in] x2 The second 2xN matrix of euclidean points
 * \param[out] M The 3x3 affine transformation matrix (6 dof)
 *          with the following parametrization
 *              |a b tx|
 *          M = |c d ty|
 *              |0 0 1 |
 *          such that    x2 = M * x1
 * \param[in] expected_precision The expected precision in order for instance
 *        to accept almost affine matrices.
 *
 * \return true if the transformation estimation has succeeded
 *
 * \note Need at least 3 non aligned points
 * \note Points coordinates must be normalized (euclidean)
 */
bool Affine2DFromCorrespondencesLinear
(
  const Mat &x1,
  const Mat &x2,
  Mat3 *M,
  double expected_precision = EigenDoubleTraits::dummy_precision()
);

/** 3D Affine transformation estimation
 *
 * This function can be used in order to estimate the affine transformation
 * between two sets of points with known 3D correspondences.
 *
 * \param[in] x1 The first 3xN matrix of euclidean points
 * \param[in] x2 The second 3xN matrix of euclidean points
 * \param[out] M The 3x3 affine transformation matrix (12 dof)
 *          with the following parametrization
 *              |a b c tx|
 *          M = |d e f ty|
 *              |g h i tz|
 *              |0 0 0 1 |
 *          such that    x2 = M * x1
 * \param[in] expected_precision The expected precision in order for instance
 *        to accept almost affine matrices.
 *
 * \return true if the transformation estimation has succeeded
 *
 * \note Need at least 4 non coplanar points
 * \note Points coordinates must be normalized (euclidean)
 */
bool Affine3DFromCorrespondencesLinear
(
  const Mat &x1,
  const Mat &x2,
  Mat4 *M,
  double expected_precision = EigenDoubleTraits::dummy_precision());
} // namespace openMVG

#endif  // OPENMVG_SOLVER_MULTIVIEW_AFFINE_HPP
