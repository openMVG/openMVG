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

#ifndef OPENMVG_MULTIVIEW_PROJECTION_HPP
#define OPENMVG_MULTIVIEW_PROJECTION_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

/// Collection of function related to the classic Projection matrix used
///  in computer vision. P = K[R|t] with [t]=[-RC] Cf HZ
namespace openMVG
{

/**
* @brief Compute P = K[R|t]
* @param K Intrinsic matrix
* @param R Rotation matrix
* @param t Translation vector
* @param[out] P Projection matrix
*/
void P_From_KRt( const Mat3 &K, const Mat3 &R, const Vec3 &t, Mat34 *P );

/**
* Decompose using the RQ decomposition
* @param P Projection matrix
* @param[out] Kp Intrinsic matrix
* @param[out] Rp Rotation matrix
* @param[out] tp Translation vector
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
* @see  HZ A4.1.1 pag.579.
*/
void KRt_From_P( const Mat34 &P, Mat3 *Kp, Mat3 *Rp, Vec3 *tp );

/**
* @brief Compute a fundamental matrix from projection matrices
* @param P1 Projection matrix of first camera
* @param P2 Projection matrix of second camera
* @return Fundamental matrix between the two camera
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
* @see  HZ Equation (17.3), page 412.
*/
Mat3 F_from_P( const Mat34 & P1, const Mat34 & P2 );

/**
* @brief Compute the depth of the X point. (R*X)[2]+t[2]
* @param R Rotation matrix
* @param t Translation vector
* @param X 3d points
* @return Depth of X wtr to the camera center
*/
double Depth( const Mat3 &R, const Vec3 &t, const Vec3 &X );

/**
* @brief Compute P*[X|1.0]. Transformed from homogeneous to euclidean coordinates
* @param P Camera projection matrix
* @param X Input 3d point
* @return Projected point
*/
Vec2 Project( const Mat34 &P, const Vec3 &X );

/**
* @brief Compute P*[X|1.0] for the X list of point (3D point)
* @param P Camera projection matrix
* @param X Input 3d points
* @param[out] x Projected points
*/
void Project( const Mat34 &P, const Mat3X &X, Mat2X *x );

/**
* @brief Compute P*X for the X list of point (4D point)
* @param P Camera projection matrix
* @param X Input 4d points
* @param[out] x Projected points
*/
void Project( const Mat34 &P, const Mat4X &X, Mat2X *x );

/**
* @brief Return P*[X|1.0] for the X list of point (3D point)
* @param P Camera projection matrix
* @param X Input 3d points
* @return Projected points
*/
Mat2X Project( const Mat34 &P, const Mat3X &X );

/**
* @brief Return P*X for the X list of point (4D point)
* @param P Camera projection matrix
* @param X Input 4d points
* @return Projected points
*/
Mat2X Project( const Mat34 &P, const Mat4X &X );

/**
* @brief Estimates the root mean square error (2D)
* @param x_image Points in image frame
* @param X_world Points in world frame
* @param P Projection matrix
* @return RMS of projection error
*/
double RootMeanSquareError( const Mat2X &x_image,
                            const Mat4X &X_world,
                            const Mat34 &P );

/**
* @brief Estimates the root mean square error (2D)
* @param x_image Points in image frame
* @param X_world Points in world frame
* @param K Intrinsic matrix
* @param R Rotation matrix
* @param t translation vector
* @note K[R|t] defines a projection
*/
double RootMeanSquareError( const Mat2X &x_image,
                            const Mat3X &X_world,
                            const Mat3 &K,
                            const Mat3 &R,
                            const Vec3 &t );

} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_PROJECTION_HPP
