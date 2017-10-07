// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_MOTION_FROM_ESSENTIAL_HPP
#define OPENMVG_MULTIVIEW_MOTION_FROM_ESSENTIAL_HPP

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{

/**
* @brief Estimate the best possible relative pose from E.
*  Four relative poses can be build from the Ematrix decomposition.
*  We keep the one with most of the point in front of the camera.
*
* @param[in] x1 bearing vectors corresponding to image observation in image 1
* @param[in] x2 bearing vectors corresponding to image observation in image 2
* @param[in] E essential matrix
* @param[in] bearing_vector_index_to_use selection of x1, x2 columns that are used
* @param[out] relative_pose the estimated relative pose
* @param[out] vec_selected_points return the index of bearing_vector_index_to_use that are
*    in front of the cameras
* @param[out] vec_points return the 3D point corresponding to
*    vec_selected_points indexes
* @param[in] positive_depth_solution_ratio Pourcentage ratio threshold used
*    to discard if there is two good solution that have many points in front
*    of the cameras
*/

bool RelativePoseFromEssential
(
  const Mat3X & x1,
  const Mat3X & x2,
  const Mat3 & E,
  const std::vector<uint32_t> & bearing_vector_index_to_use,
  geometry::Pose3 * relative_pose = nullptr,
  std::vector<uint32_t> * vec_selected_points = nullptr,
  std::vector<Vec3> * vec_points = nullptr,
  const double positive_depth_solution_ratio = 0.7
);

} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_MOTION_FROM_ESSENTIAL_HPP
