
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

#ifndef OPENMVG_MULTIVIEW_PROJECTION_H_
#define OPENMVG_MULTIVIEW_PROJECTION_H_

#include "openMVG/numeric/numeric.h"

/// Collection of function related to the classic Projection matrix used
///  in computer vision. P = K[R|t] with [t]=[-RC] Cf HZ
namespace openMVG {

/// Compute P = K[R|t]
void P_From_KRt(const Mat3 &K, const Mat3 &R, const Vec3 &t,Mat34 *P);

/// Decompose using the RQ decomposition HZ A4.1.1 pag.579.
void KRt_From_P(const Mat34 &P, Mat3 *Kp, Mat3 *Rp, Vec3 *tp);

/// Compute a fundamental matrix from projection matrices
Mat3 F_from_P(const Mat34 & P1, const Mat34 & P2);

// Compute the depth of the X point. R*X[2]+t[2].
double Depth(const Mat3 &R, const Vec3 &t, const Vec3 &X);

// Compute P*[X|1.0]. Transformed from homogeneous to euclidean coordinates.
Vec2 Project(const Mat34 &P, const Vec3 &X);

// Compute P*[X|1.0] for the X list of point (3D point).
void Project(const Mat34 &P, const Mat3X &X, Mat2X *x);

// Compute P*[X|1.0] for the X list of point (4D point).
void Project(const Mat34 &P, const Mat4X &X, Mat2X *x);

// Return P*[X|1.0] for the X list of point (3D point).
Mat2X Project(const Mat34 &P, const Mat3X &X);

// Return P*[X|1.0] for the X list of point (4D point).
Mat2X Project(const Mat34 &P, const Mat4X &X);

// Change homogeneous coordinates to euclidean.
void HomogeneousToEuclidean(const Vec4 &H, Vec3 *X);

// Change euclidean coordinates to homogeneous.
void EuclideanToHomogeneous(const Mat &X, Mat *H);

// Change euclidean coordinates to homogeneous.
Vec3 EuclideanToHomogeneous(const Vec2 &x);

void HomogeneousToEuclidean(const Mat &H, Mat *X);

Mat3X EuclideanToHomogeneous(const Mat2X &x);

void EuclideanToHomogeneous(const Mat2X &x, Mat3X *h);

void HomogeneousToEuclidean(const Mat3X &h, Mat2X *e);

/// Project x point in camera coordinates
void EuclideanToNormalizedCamera(const Mat2X &x, const Mat3 &K, Mat2X *n);

/// Project x point in camera coordinates
void HomogeneousToNormalizedCamera(const Mat3X &x, const Mat3 &K, Mat2X *n);

/// Estimates the root mean square error (2D)
double RootMeanSquareError(const Mat2X &x_image,
  const Mat4X &X_world,
  const Mat34 &P);

/// Estimates the root mean square error (2D)
double RootMeanSquareError(const Mat2X &x_image,
  const Mat3X &X_world,
  const Mat3 &K,
  const Mat3 &R,
  const Vec3 &t);

} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_PROJECTION_H_
