
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
  std::vector<IndMatch> & vec_PutativeMatches)
{
    //-- Copy to contiguous array for matching
    //-- In order to be generic as possible, the array matcher use raw pointer.
    typedef typename DESCRIPTOR_TYPE::bin_type descTbin;
    const size_t DSIZE = DESCRIPTOR_TYPE::static_size;
    descTbin * tab0 = new descTbin[descsL.size()*DSIZE];
    descTbin * tab1 = new descTbin[descsR.size()*DSIZE];
    for(size_t i=0; i < descsL.size(); ++i)
      memcpy(&tab0[i*DSIZE],descsL[i].getData(),DSIZE*sizeof(typename DESCRIPTOR_TYPE::bin_type));

    for(size_t i=0; i < descsR.size(); ++i)
      memcpy(&tab1[i*DSIZE],descsR[i].getData(),DSIZE*sizeof(typename DESCRIPTOR_TYPE::bin_type));

    const int NNN__ = 2; // look for the two nearest distance for each Left descriptor
    std::vector<int> vec_nIndice;
    std::vector<typename MatcherT::MetricT::ResultType> vec_Distance;

    MatcherT matcher;
    matcher.Build(tab0, descsL.size(), DESCRIPTOR_TYPE::static_size);
    matcher.SearchNeighbours(tab1, descsR.size(), &vec_nIndice, &vec_Distance, NNN__);

    // Free temporary used memory
    delete [] tab0;
    delete [] tab1;

    // Filter the correspondences with the Ratio of distance
    std::vector<int> vec_loweRatioIndexes;
    NNdistanceRatio( vec_Distance.begin(), vec_Distance.end(),
      NNN__, // Number of neighbor in iterator sequence (minimum required 2)
      vec_loweRatioIndexes, // output (index that respect Lowe Ratio)
      fNearestNeighborDistanceRatio);

    // Get back feature index : (Left, right) index.
    for (size_t k=0; k < vec_loweRatioIndexes.size()-1
        && vec_loweRatioIndexes.size()>0; ++k) {
      vec_PutativeMatches.push_back(
      IndMatch(vec_nIndice[vec_loweRatioIndexes[k]*NNN__],
               vec_loweRatioIndexes[k]) );
    }
  }

