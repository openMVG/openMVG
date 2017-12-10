
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

#ifndef OPENMVG_MULTIVIEW_ESSENTIAL_HPP
#define OPENMVG_MULTIVIEW_ESSENTIAL_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/geometry/pose3.hpp"

namespace openMVG
{

/**
* @brief Compute the relative camera motion between two cameras.
* Given the motion parameters of two cameras, computes the motion parameters
* of the second one assuming the first one to be at the origin.
* If T1 and T2 are the camera motions, the computed relative motion is
*
*      T = T2 T1^{-1}
* @param R1 Rotation of first camera
* @param t1 Translation of first camera
* @param R2 Rotation of second camera
* @param t2 Translation of second camera
* @param[out] R Relative rotation between camera 1 and camera 2
* @param[out] t Relative translation between camera 1 and camera 2
*/
void RelativeCameraMotion( const Mat3 &R1,
                           const Vec3 &t1,
                           const Mat3 &R2,
                           const Vec3 &t2,
                           Mat3 *R,
                           Vec3 *t );

/**
* @brief Given F, Left/Right K matrix it compute the Essential matrix
* @param F Fundamental matrix
* @param K1 Intrinsic matrix of first camera
* @param K2 Intrinsic matrix of second camera
* @param[out] E Essential matrix
*/
void EssentialFromFundamental( const Mat3 &F,
                               const Mat3 &K1,
                               const Mat3 &K2,
                               Mat3 *E );

/**
* @brief Compute E as E = [t12]x R12.
* @param R1 First camera rotation matrix
* @param t1 First camera translation vector
* @param R2 Second camera rotation matrix
* @param t2 Second camera translation vector
* @param[out] Essential matrix
*/
void EssentialFromRt( const Mat3 &R1,
                      const Vec3 &t1,
                      const Mat3 &R2,
                      const Vec3 &t2,
                      Mat3 *E );

/**
* @brief Given E, Left/Right K matrix it compute the Fundamental matrix
* @param E Essential matrix
* @param K1 Intrinsic matrix of first camera
* @param K2 Intrinsic matrix of second camera
* @param[out] F Fundamental matrix
*/
void FundamentalFromEssential( const Mat3 &E,
                               const Mat3 &K1,
                               const Mat3 &K2,
                               Mat3 *F );

/**
 * @brief Given an essential matrix computes the 2 rotations and the 2 translations
 * that can generate four possible motions.
 * @param E Essential matrix
 * @param[out] relative_poses The 4 possible relative poses
 * @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
 * @see HZ 9.7 page 259 (Result 9.19)
 */
void MotionFromEssential(
  const Mat3 &E,
  std::vector<geometry::Pose3> * relative_poses );


} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_ESSENTIAL_HPP
