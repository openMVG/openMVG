// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_RELATIVE_POSE_ENGINE_HPP
#define OPENMVG_SFM_RELATIVE_POSE_ENGINE_HPP

#include "openMVG/types.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/triangulation_method.hpp"

namespace openMVG {
namespace sfm {

struct SfM_Data;
struct Matches_Provider;
struct Features_Provider;

/// An engine to compute relative pose
class Relative_Pose_Engine
{
public:
  using Relative_Pair_Poses = Hash_Map<Pair, geometry::Pose3>;

  Relative_Pose_Engine () = default;

  // Try to compute all the possible relative pose.
  bool Process(
    const SfM_Data & sfm_data_,
    const Matches_Provider * matches_provider_,
    const Features_Provider * features_provider_
  );

  // Try to compute the depicted relative pose pairs
  bool Process(
    const Pair_Set & relative_pose_pairs,
    const SfM_Data & sfm_data_,
    const Matches_Provider * matches_provider_,
    const Features_Provider * features_provider_
  );

  // Relative poses accessor
  const Relative_Pair_Poses& Get_Relative_Poses() const;

  /// Configure the 2view triangulation method used by the RelativePose computing engine
  void SetTriangulationMethod(const ETriangulationMethod method)
  {
    triangulation_method_ = method;
  }

private:
  Relative_Pair_Poses relative_poses_;

  ETriangulationMethod triangulation_method_ = ETriangulationMethod::DEFAULT;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_RELATIVE_POSE_ENGINE_HPP
