
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/features/regions.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

namespace openMVG {

enum EMatcherType
{
  BRUTE_FORCE_L2,
  ANN_L2,
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

  /// Load all features and descriptors in memory
  bool loadData(
    std::unique_ptr<features::Regions>& region_type, // interface to load computed regions
    const std::vector<std::string> & vec_fileNames, // input filenames
    const std::string & sMatchDir); // where the data are saved

  void Match(
    const std::vector<std::string> & vec_fileNames, // input filenames,
    const Pair_Set & pairs,
    matching::PairWiseMatches & map_PutativesMatches)const; // the pairwise photometric corresponding points

  private:
  std::map<IndexT, std::unique_ptr<features::Regions> > regions_perImage;
  // Distance ratio used to discard spurious correspondence
  float fDistRatio;
  // Matcher Type
  EMatcherType _eMatcherType;
};

}; // namespace openMVG
