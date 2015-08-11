
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/features/regions.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

namespace openMVG {
namespace matching_image_collection {

enum EMatcherType
{
  BRUTE_FORCE_L2,
  ANN_L2,
  CASCADE_HASHING_L2,
  BRUTE_FORCE_HAMMING
};

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
/// Spurious correspondences are discarded by using the
///  a threshold over the distance ratio of the 2 neighbours points.
///
class Matcher_Regions_AllInMemory : public Matcher
{
  public:
  Matcher_Regions_AllInMemory(float distRatio, EMatcherType eMatcherType);

  /// Find corresponding points between some pair of view Ids
  void Match(
    const sfm::SfM_Data & sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair_Set & pairs,
    matching::PairWiseMatches & map_PutativesMatches)const; // the pairwise photometric corresponding points

  private:
  std::map<IndexT, std::unique_ptr<features::Regions> > regions_perImage;
  // Distance ratio used to discard spurious correspondence
  float fDistRatio;
  // Matcher Type
  EMatcherType _eMatcherType;
};

} // namespace openMVG
} // namespace matching_image_collection
