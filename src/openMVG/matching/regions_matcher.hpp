
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "indMatch.hpp"
#include "indMatchDecoratorXY.hpp"
#include "matching_filters.hpp"

#include "openMVG/numeric/numeric.h"
#include "openMVG/features/regions.hpp"

#include <vector>

namespace openMVG {
namespace matching {


/**
 * Match one View Regions with other Regions.
 */
template < class ArrayMatcherT >
class RegionsMatcher
{
private:
  ArrayMatcherT _matcher;
  const features::Regions& _origRegions;

public:
  typedef typename ArrayMatcherT::ScalarT Scalar;
  typedef typename ArrayMatcherT::DistanceType DistanceType;

  /**
   * @brief Create the matcher with a reference image "I".
   */
  RegionsMatcher(const features::Regions& regionsI)
    : _origRegions(regionsI)
  {
    if (regionsI.RegionCount() == 0)
      return;

    const Scalar * tab = reinterpret_cast<const Scalar *>(regionsI.DescriptorRawData());
    _matcher.Build(tab, regionsI.RegionCount(), regionsI.DescriptorLength());
  }

  /**
   * @brief Match a Region group "J" (aka. an image) with the reference.
   */
  bool MatchRatioTest(
    std::vector<IndMatch>& vec_FilteredMatches,
    const features::Regions& regionsJ,
    const float fDistRatio)
  {
      const Scalar * tab = reinterpret_cast<const Scalar *>(regionsJ.DescriptorRawData());

      const size_t NNN__ = 2;
      std::vector<int> vec_nIndice10;
      std::vector<DistanceType> vec_fDistance10;

      // Search the 2 closest features neighbours for each feature
      if (!_matcher.SearchNeighbours(tab, regionsJ.RegionCount(), &vec_nIndice10, &vec_fDistance10, NNN__))
        return false;
      
      std::vector<int> vec_NNRatioIndexes;
      // Filter the matches using the ratio test:
      //   The probability that a match is correct can be determined by taking
      //   the ratio of distance from the closest neighbor to the distance
      //   of the second closest.
      NNdistanceRatio(
        vec_fDistance10.begin(), // distance start
        vec_fDistance10.end(),  // distance end
        NNN__, // Number of neighbor in iterator sequence (minimum required 2)
        vec_NNRatioIndexes, // output (indices that respect Lowe Ratio)
        fDistRatio);

      for (size_t k=0; k < vec_NNRatioIndexes.size(); ++k)
      {
        const size_t index = vec_NNRatioIndexes[k];
        vec_FilteredMatches.push_back(
          IndMatch(vec_nIndice10[index*NNN__], index) );
      }

      // Remove duplicates
      IndMatch::getDeduplicated(vec_FilteredMatches);

      // Remove matches that have the same (X,Y) coordinates
      IndMatchDecorator<float> matchDeduplicator(vec_FilteredMatches, _origRegions.GetRegionsPositions(), regionsJ.GetRegionsPositions());
      matchDeduplicator.getDeduplicated(vec_FilteredMatches);
  }
};

}  // namespace matching
}  // namespace openMVG
