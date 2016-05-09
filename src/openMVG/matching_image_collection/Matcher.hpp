
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

#include <string>
#include <vector>

namespace openMVG {
namespace matching_image_collection {

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
class Matcher
{
  public:
  Matcher() = default ;

  virtual ~Matcher() = default ;

  /// Find corresponding points between some pair of view Ids
  virtual void Match(
    const sfm::SfM_Data & sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair_Set & pairs, // list of pair to consider for matching
    matching::PairWiseMatches & map_putatives_matches // the output pairwise photometric corresponding points
    )const = 0;
};

} // namespace matching_image_collection
} // namespace openMVG 
