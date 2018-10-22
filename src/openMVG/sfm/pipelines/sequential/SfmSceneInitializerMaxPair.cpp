
// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerMaxPair.hpp"

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/stl/stl.hpp"

namespace openMVG {

using namespace cameras;
using namespace geometry;
using namespace matching;

namespace sfm {

SfMSceneInitializerMaxPair::SfMSceneInitializerMaxPair(
  SfM_Data & sfm_data,
  const Features_Provider * features_provider,
  const Matches_Provider * matches_provider)
  :SfMSceneInitializer(sfm_data, features_provider, matches_provider)
{
  sfm_data_.poses.clear();
}

bool SfMSceneInitializerMaxPair::Process()
{
  if (sfm_data_.GetIntrinsics().empty())
    return false;

  //
  // Sort the PairWiseMatches by matches count.
  // Keep the first pair that provides a valid relative pose.
  //
  std::vector<IndexT> matches_count_per_pair;
  matches_count_per_pair.reserve(matches_provider_->pairWise_matches_.size());
  std::transform(matches_provider_->pairWise_matches_.cbegin(), matches_provider_->pairWise_matches_.cend(),
    std::back_inserter(matches_count_per_pair),
    [](const std::pair<Pair, IndMatches> & match_pair) -> IndexT { return match_pair.second.size(); });

  // sort the Pairs in descending order according their matches count
  using namespace stl::indexed_sort;
  std::vector<sort_index_packet_descend<IndexT, IndexT>> packet_vec(matches_count_per_pair.size());
  sort_index_helper(packet_vec, &matches_count_per_pair[0], std::min((IndexT)10, (IndexT)matches_count_per_pair.size()));

  // Iterate through the pairs and try to find the relative pose
  for (size_t i = 0; i < matches_count_per_pair.size(); ++i)
  {
    const IndexT index = packet_vec[i].index;
    openMVG::matching::PairWiseMatches::const_iterator iter = matches_provider_->pairWise_matches_.cbegin();
    std::advance(iter, index);

    std::cout << "(" << iter->first.first << "," << iter->first.second <<"): "
      << iter->second.size() << " matches" << std::endl;

    const IndexT
      I = iter->first.first,
      J = iter->first.second;

    const View
      * view_I = sfm_data_.views[I].get(),
      * view_J = sfm_data_.views[J].get();

    // Check that the pair has valid intrinsic
    if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
        sfm_data_.GetIntrinsics().count(view_J->id_intrinsic) == 0)
      continue;

    const IntrinsicBase
      * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get(),
      * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

    // Compute for each feature the un-distorted camera coordinates
    const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(iter->first);
    size_t number_matches = matches.size();
    Mat2X x1(2, number_matches), x2(2, number_matches);
    number_matches = 0;
    for (const auto & match : matches)
    {
      x1.col(number_matches) = cam_I->get_ud_pixel(features_provider_->feats_per_view.at(I)[match.i_].coords().cast<double>());
      x2.col(number_matches) = cam_J->get_ud_pixel(features_provider_->feats_per_view.at(J)[match.j_].coords().cast<double>());
      ++number_matches;
    }

    RelativePose_Info relativePose_info;
    relativePose_info.initial_residual_tolerance = Square(2.5);
    if (!robustRelativePose(cam_I, cam_J,
                            x1, x2, relativePose_info,
                            {cam_I->w(), cam_I->h()},
                            {cam_J->w(), cam_J->h()},
                            2048))
    {
      continue;
    }
    sfm_data_.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    sfm_data_.poses[view_J->id_pose] = relativePose_info.relativePose;
    return true;
  }
  return false;
}

} // namespace sfm
} // namespace openMVG
