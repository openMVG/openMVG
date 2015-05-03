
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/matching_filters.hpp"
using namespace openMVG;
using namespace openMVG::matching;

/// Return the corresponding nearest neighbor in descsR for descsL descriptor
/// Filtered with Nearest neighbor distance ratio in order to remove too
///  similar looking matches
template< typename DESCRIPTOR_TYPE, typename MatcherT>
void getPutativesMatches(
  const std::vector<DESCRIPTOR_TYPE > & descsL,
  const std::vector<DESCRIPTOR_TYPE > & descsR,
  const float fNearestNeighborDistanceRatio,
  IndMatches & vec_PutativeMatches)
{
  typedef typename DESCRIPTOR_TYPE::bin_type DescTbin;

  const int NNN__ = 2; // look for the two nearest distance for each Left descriptor
  std::vector<int> vec_nIndice;
  std::vector<typename MatcherT::MetricT::ResultType> vec_Distance;

  MatcherT matcher;
  const DescTbin* descsLReinterpret =
    reinterpret_cast<const DescTbin *>(&descsL[0]);

  const DescTbin* descsRReinterpret =
    reinterpret_cast<const DescTbin *>(&descsR[0]);

  matcher.Build(descsLReinterpret, descsL.size(), DESCRIPTOR_TYPE::static_size);
  matcher.SearchNeighbours(descsRReinterpret, descsR.size(), &vec_nIndice, &vec_Distance, NNN__);

  // Filter the correspondences with the Ratio of distance
  std::vector<int> vec_loweRatioIndexes;
  NNdistanceRatio( vec_Distance.begin(), vec_Distance.end(),
    NNN__, // Number of neighbor in iterator sequence (minimum required 2)
    vec_loweRatioIndexes, // output (index that respect Lowe Ratio)
    fNearestNeighborDistanceRatio);

  vec_PutativeMatches.reserve(vec_loweRatioIndexes.size());
  // Get back feature index : (Left, right) index.
  for (size_t k=0; k < vec_loweRatioIndexes.size(); ++k) {
    vec_PutativeMatches.push_back(
    IndMatch(vec_nIndice[vec_loweRatioIndexes[k]*NNN__],
             vec_loweRatioIndexes[k]) );
  }
}
