// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_HPP
#define OPENMVG_MULTIVIEW_TRIANGULATION_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{

/**
* @brief Linear DLT triangulation
* @brief P1 First camera projection matrix
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
* @brief P1 First camera projection matrix
* @param P2 Second camera projection matrix
* @param x1 bearing vector of the landmark observation in the first camera
* @param x2 bearing vector of the landmark observation in the second camera
* @param[out] X_euclidean Euclidean triangulated point
* @see HZ 12.2 pag.312
* @ref Multiple View Geometry - Richard Hartley, Andrew Zisserman - second edition
*/
void TriangulateDLT
( const Mat34 &P1,
  const Vec3 &x1,
  const Mat34 &P2,
  const Vec3 &x2,
  Vec3 *X_euclidean
);

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_HPP
