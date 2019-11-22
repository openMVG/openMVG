// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_HPP
#define OPENMVG_MULTIVIEW_TRIANGULATION_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/multiview/triangulation_method.hpp"

namespace openMVG
{

/**
* @brief Generic triangulation routine (Aggregate all the 2view triangulation solvers in one place).
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @return true if the point pass the adequacy or cheirality test (depending of the solver), false otherwise
* (invalid 3d point or method that does not exist)
**/
bool Triangulate2View
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &bearing0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &bearing1,
  Vec3 &X,
  ETriangulationMethod etri_method = ETriangulationMethod::DEFAULT
);

/**
* @brief Linear DLT triangulation
* @param P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 bearing vector of the landmark observation in the first camera
* @param x2 bearing vector of the landmark observation in the second camera
* @param[out] X_homogeneous Homogeneous triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT
(
  const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec4 *X_homogeneous
);

/**
* @brief Linear DLT triangulation
* @param P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 bearing vector of the landmark observation in the first camera
* @param x2 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT
(
  const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec3 *X_euclidean
);

/**
* @brief Linear DLT triangulation
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @return true if the point pass the cheirality test, false otherwise
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
bool TriangulateDLT
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
);

/**
* @brief Optimal L1 Angular triangulation
* @brief Minimize the L1 norm of angular errors
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @return true if the point pass the cheirality test, false otherwise
* @ref S.H. Lee, J. Civera - Closed-Form Optimal Triangulation Based on Angular Errors - ICCV 2019 - https://arxiv.org/pdf/1903.09115.pdf
*/
bool TriangulateL1Angular
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
);


/**
* @brief Optimal LInfinity Angular triangulation
* @brief Minimize the LInfinity norm of angular errors
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @return true if the point pass the cheirality test, false otherwise
* @ref S.H. Lee, J. Civera - Closed-Form Optimal Triangulation Based on Angular Errors - ICCV 2019 - https://arxiv.org/pdf/1903.09115.pdf
*/
bool TriangulateLInfinityAngular
(
  const Mat3 &R0,
  const Vec3 &t0,
  const Vec3 &x0,
  const Mat3 &R1,
  const Vec3 &t1,
  const Vec3 &x1,
  Vec3 *X_euclidean
);

/**
* @brief Inverse Depth Weighted Midpoint method
* @brief and its ad hoc adequacy test (a replacement for cheirality tests)
* @brief should be better than DLT for low and high parallax angles
* @param R0 First Camera rotation matrix
* @param t0 First Camera translation vector
* @param x0 bearing vector of the landmark observation in the first camera
* @param R1 Second Camera rotation matrix
* @param t1 Second Camera translation vector
* @param x1 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @return true if the point pass the adequacy test, false otherwise
* @ref S.H. Lee, J. Civera - Triangulation: Why Optimize? - BMVC 2019 - https://arxiv.org/pdf/1907.11917.pdf
*/
bool TriangulateIDWMidpoint(
  const Mat3 & R0,
  const Vec3 & t0,
  const Vec3 & x0,
  const Mat3 & R1,
  const Vec3 & t1,
  const Vec3 & x1,
  Vec3 *X_euclidean
);

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_HPP
