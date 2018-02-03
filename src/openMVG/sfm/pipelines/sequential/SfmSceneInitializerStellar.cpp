
// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerStellar.hpp"

#include "openMVG/sfm/pipelines/relative_pose_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/stellar/stellar_solver.hpp"
#include "openMVG/sfm/sfm_data.hpp"

namespace openMVG {

using namespace matching;

namespace sfm {

SfMSceneInitializerStellar::SfMSceneInitializerStellar(
  SfM_Data & sfm_data,
  const Features_Provider * features_provider,
  const Matches_Provider * matches_provider)
  :SfMSceneInitializer(sfm_data, features_provider, matches_provider)
{
  sfm_data_.poses.clear();
}

bool SfMSceneInitializerStellar::Process()
{
  if (sfm_data_.GetIntrinsics().empty())
    return false;

  // Compute a relative pose for each edge of the pose pair graph
  const Relative_Pose_Engine::Relative_Pair_Poses relative_poses = [&]
  {
    Relative_Pose_Engine relative_pose_engine;
    if (!relative_pose_engine.Process(sfm_data_,
        matches_provider_,
        features_provider_))
      return Relative_Pose_Engine::Relative_Pair_Poses();
    else
      return relative_pose_engine.Get_Relative_Poses();
  }();

  // List all possible stellar configurations
  using StellarPods = Hash_Map<IndexT, Pair_Set>;
  StellarPods stellar_pods;
  for (const auto & it : relative_poses)
  {
    const Pair & pairIt = it.first;
    stellar_pods[pairIt.first].insert(pairIt);
    stellar_pods[pairIt.second].insert(pairIt);
  }

  // Keep the possible stellar configuration that has the most of average matches
  std::vector<float> matches_count_per_stellar_pod(stellar_pods.size(), 0);
  IndexT id = 0;
  for (const auto & stellar_pod : stellar_pods)
  {
    const Pair_Set & pairs = stellar_pod.second;
    for (const auto & pair : pairs)
    {
      const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(pair);
      matches_count_per_stellar_pod[id] += matches.size();
    }
    matches_count_per_stellar_pod[id] /= pairs.size();
    ++id;
  }

  for (const auto & val : matches_count_per_stellar_pod)
  {
    std::cout << val << std::endl;
  }

  const auto stellar_pod_max_matches_iterator =
    std::max_element(matches_count_per_stellar_pod.cbegin(),
                     matches_count_per_stellar_pod.cend());

  std::cout << *stellar_pod_max_matches_iterator << std::endl;

  if (stellar_pod_max_matches_iterator == matches_count_per_stellar_pod.cend())
    return false;

  auto stellar_pod_it = std::next(stellar_pods.cbegin(), std::distance(matches_count_per_stellar_pod.cbegin(), stellar_pod_max_matches_iterator));

  std::cout << "The choosen stellar pod has " << stellar_pod_it->second.size() << " pairs." << std::endl;

  // Configure the stellar pod optimization to use all the matches relating to the considered view id.
  // The found camera poses will be better thanks to the larger point support.
  const bool use_all_matches = true;
  const bool use_threading = true;
  Stellar_Solver stellar_pod_solver(
    stellar_pod_it->second,
    relative_poses,
    sfm_data_,
    matches_provider_,
    features_provider_,
    use_all_matches,
    use_threading);

  if (stellar_pod_solver.Solve(sfm_data_.poses))
  {
    return true;
  }

  return false;
}

} // namespace sfm
} // namespace openMVG
