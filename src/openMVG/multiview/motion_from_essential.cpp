// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/motion_from_essential.hpp"

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/extract_columns.hpp"

using namespace openMVG::cameras;
using namespace openMVG::geometry;

namespace openMVG
{

bool RelativePoseFromEssential
(
  const Mat3X & x1,
  const Mat3X & x2,
  const Mat3 & E,
  const std::vector<uint32_t> & bearing_vector_index_to_use,
  Pose3 * relative_pose,
  std::vector<uint32_t> * vec_selected_points,
  std::vector<Vec3> * vec_points,
  const double positive_depth_solution_ratio
)
{
  // Recover plausible relative poses from E.
  std::vector<Pose3> relative_poses;
  MotionFromEssential(E, &relative_poses);

  // Accumulator to find the best solution
  std::vector<uint32_t> cheirality_accumulator(relative_poses.size(), 0);

  // Find which solution is the best:
  // - count how many triangulated observations are in front of the cameras
  std::vector<std::vector<uint32_t>> vec_newInliers(relative_poses.size());
  std::vector<std::vector<Vec3>> vec_3D(relative_poses.size());

  const Pose3 pose1(Mat3::Identity(), Vec3::Zero());
  const Mat34 P1 = pose1.asMatrix();

  for (size_t i = 0; i < relative_poses.size(); ++i)
  {
    const Pose3 pose2 = relative_poses[i];
    const Mat34 P2 = pose2.asMatrix();
    Vec3 X;

    for (const uint32_t inlier_idx : bearing_vector_index_to_use)
    {
      const auto
        f1 = x1.col(inlier_idx),
        f2 = x2.col(inlier_idx);
      TriangulateDLT(P1, f1, P2, f2, &X);
      // Test if X is visible by the two cameras
      if (CheiralityTest(f1, pose1, f2, pose2, X))
      {
        ++cheirality_accumulator[i];
        vec_newInliers[i].push_back(inlier_idx);
        vec_3D[i].push_back(X);
      }
    }
  }

  // Check if there is a valid solution:
  const auto iter = std::max_element(cheirality_accumulator.cbegin(), cheirality_accumulator.cend());
  if (*iter == 0)
  {
    // There is no right solution with points in front of the cameras
    return false;
  }

  // Export the best solution data
  const size_t index = std::distance(cheirality_accumulator.cbegin(), iter);
  if (relative_pose)
  {
    (*relative_pose) = relative_poses[index];
  }
  if (vec_selected_points)
    (*vec_selected_points) = vec_newInliers[index];
  if (vec_points)
    (*vec_points) = vec_3D[index];

  // Test if the best solution is good by using the ratio of the two best solution score
  std::sort(cheirality_accumulator.begin(), cheirality_accumulator.end());
  const double ratio = cheirality_accumulator.rbegin()[1]
    / static_cast<double>(cheirality_accumulator.rbegin()[0]);
  return (ratio < positive_depth_solution_ratio);
}

} // namespace openMVG
