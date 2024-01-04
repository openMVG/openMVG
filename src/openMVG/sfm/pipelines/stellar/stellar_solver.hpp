// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_STELLAR_STELLAR_SOLVER_HPP
#define OPENMVG_SFM_STELLAR_STELLAR_SOLVER_HPP

#include "openMVG/sfm/pipelines/stellar/relative_scale.hpp"
#include "openMVG/sfm/sfm_data.hpp"

namespace openMVG {
namespace sfm {

struct Matches_Provider;
struct Features_Provider;

using StellarPod = Pair_Set;

struct Stellar_Solver
{
  Stellar_Solver
  (
    const StellarPod & stellar_pod,
    const Hash_Map<Pair, Pose3> & relative_poses,
    const SfM_Data & sfm_data,
    const Matches_Provider * matches_provider,
    const Features_Provider * features_provider,
    const bool use_all_matches = false,
    const bool use_threading = false
  );

  bool Solve(Poses & poses);

private:

  std::vector<Pair_Set> ListEdge2Uplets();

  // Solve the relative scale for some triplet of poses (2-uplet of edges)
  bool Solve2UpletsRelativeScales
  (
    const std::vector<Pair_Set> & edge_two_uplets,
    std::vector<Relative_Scale> & relative_scales
  );

  // Solve the local to global coordinates of the poses.
  // Re-conciliate the translation relative scale to the same coordinate system.
  bool SolveStellarPoses
  (
    const IndexT & central_node_id,
    const std::vector<Relative_Scale> & relative_scales,
    Hash_Map<IndexT, geometry::Pose3> & triplet_poses
  );

  /// Optimize a stellar pod by using track triangulation & BA (Bundle Adjustment)
  /// Default mode use only the tracks linked to the triplet pair
  /// The b_use_all_matches make usage of all the tracks linked to the used view index
  ///   so a larger number of tracks can be used. It enhances the results, but make the BA slower.
  bool Optimize
  (
    SfM_Data & stellar_pod_reconstruction,
    const Hash_Map<IndexT, geometry::Pose3> & triplet_poses,
    const std::vector<Relative_Scale> & relative_scales
  );

private:
  StellarPod stellar_pod_;
  // Store if we use all the matches relating to the considered view id,
  // or only the one linked to the input relative pairs.
  bool use_all_matches_;
  bool use_threading_;

  // General property from the scene:
  const SfM_Data & sfm_data_;
  const Matches_Provider * matches_provider_;
  const Features_Provider * features_provider_;
  // Relative motion cache:
  const Hash_Map<Pair, Pose3> & relative_poses_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_STELLAR_STELLAR_DEFINITIONS_HPP
