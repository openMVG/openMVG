// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_SIMILARITY3_HPP
#define OPENMVG_GEOMETRY_SIMILARITY3_HPP

#include "openMVG/geometry/pose3.hpp"

namespace openMVG
{
namespace geometry
{

/**
* @brief Define a 3D Similarity transform encoded as a 3D pose plus a scale
*/
struct Similarity3
{
  /// Pose
  Pose3 pose_;

  /// Scale
  double scale_;

  /**
  * @brief Default constructor
  * @note This define identity transformation centered at origin of a cartesian frame
  */
  Similarity3();

  /**
  * @brief Constructor
  * @param pose a 3d pose
  * @param scale a scale factor
  */
  Similarity3( const Pose3 & pose, const double scale );

  /**
  * @brief Apply transformation to a point
  * @param point Input point
  * @return transformed point
  */
  Mat3X operator () ( const Mat3X & point ) const;

  /**
  * @brief Concatenation of pose
  * @param pose Pose to be concatenated with the current one
  * @return Concatenation of poses
  */
  Pose3 operator () ( const Pose3 & pose ) const;

  /**
  * @brief Get inverse of the similarity
  * @return Inverse of the similarity
  */
  Similarity3 inverse() const;

};

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_SIMILARITY3_HPP
