
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/matching_image_collection/Matcher.hpp"

namespace openMVG {
namespace matching_image_collection {

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
/// Spurious correspondences are discarded by using the
///  a threshold over the distance ratio of the 2 nearest neighbours.
///
class Matcher_Regions_AllInMemory : public Matcher
{
  public:
  Matcher_Regions_AllInMemory
  (
    float dist_ratio,
    matching::EMatcherType eMatcherType
  );

  /// Find corresponding points between some pair of view Ids
  void Match
  (
    const sfm::SfM_Data & sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair_Set & pairs,
    matching::PairWiseMatches & map_PutativesMatches // the pairwise photometric corresponding points
  )const;

  private:
  // Distance ratio used to discard spurious correspondence
  float _f_dist_ratio;
  // Matcher Type
  matching::EMatcherType _eMatcherType;
};

} // namespace openMVG
} // namespace matching_image_collection
