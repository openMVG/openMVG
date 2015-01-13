
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/matching/indMatch.hpp"
//#include "openMVG/matching_image_collection/Pair_Builder.hpp"

#include <string>
#include <vector>

namespace openMVG {

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
class Matcher
{
  public:
  Matcher() {};

  virtual ~Matcher() {};

  // Visidyn shit
  //virtual void Match(
    

  /// Build point indexes correspondences lists between images ids
  virtual void Match(
    const std::vector<std::string> & vec_filenames,
    //const PairsT & pairs, // list of pair to consider for matching
    matching::PairWiseMatches & map_putatives_matches // the output pairwise photometric corresponding points
    )const = 0;
};

}; // namespace openMVG
