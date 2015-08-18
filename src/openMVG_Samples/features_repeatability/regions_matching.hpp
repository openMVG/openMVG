
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/features/features.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"

namespace openMVG{

/// Template Regions nearest neighbor matching with distance ratio filter
template <typename MatcherT>
void Template_Matcher(
  const features::Regions * regionsI,
  const features::Regions * regionsJ,
  IndMatches & matches, // photometric corresponding points
  // Distance ratio used to discard spurious correspondence
  const float fDistRatio)
{
  const size_t regions_countI = regionsI->RegionCount();
  if (regions_countI == 0)
  {
    return;
  }

  const std::vector<features::PointFeature> pointFeaturesI = regionsI->GetRegionsPositions();
  const typename MatcherT::ScalarT * tabI =
    reinterpret_cast<const typename MatcherT::ScalarT *>(regionsI->DescriptorRawData());

  MatcherT matcher10;
  if (!matcher10.Build(tabI, regions_countI, regionsI->DescriptorLength()))
  {
    return;
  }

  const size_t regions_countJ = regionsJ->RegionCount();
  if(regions_countJ == 0)
  {
    return;
  }

  const typename MatcherT::ScalarT * tabJ =
    reinterpret_cast<const typename MatcherT::ScalarT *>(regionsJ->DescriptorRawData());

  const size_t NNN__ = 2;
  matching::IndMatches vec_nIndice10;
  std::vector<typename MatcherT::DistanceType> vec_fDistance10;

  //Find left->right
  if (matcher10.SearchNeighbours(tabJ, regions_countJ, &vec_nIndice10, &vec_fDistance10, NNN__))
  {
    IndMatches vec_FilteredMatches;
    std::vector<int> vec_NNRatioIndexes;
    NNdistanceRatio(
      vec_fDistance10.begin(), // distance start
      vec_fDistance10.end(),  // distance end
      NNN__, // Number of neighbor in iterator sequence (minimum required 2)
      vec_NNRatioIndexes, // output (indices that respect Lowe Ratio)
      fDistRatio);

    for (size_t k=0; k < vec_NNRatioIndexes.size(); ++k)
    {
      const size_t index = vec_NNRatioIndexes[k];
      vec_FilteredMatches.emplace_back(
        IndMatch(vec_nIndice10[index*NNN__]._j, vec_nIndice10[index*NNN__]._i));
    }

    // Remove duplicates
    IndMatch::getDeduplicated(vec_FilteredMatches);

    // Remove matches that have the same (X,Y) coordinates
    const std::vector<features::PointFeature> pointFeaturesJ = regionsJ->GetRegionsPositions();
    IndMatchDecorator<float> matchDeduplicator(vec_FilteredMatches, pointFeaturesI, pointFeaturesJ);
    matchDeduplicator.getDeduplicated(vec_FilteredMatches);

    if (!vec_FilteredMatches.empty())
    {
      matches = std::move(vec_FilteredMatches);
    }
  }
}

/// Generic function to match regions (use here a BRUTE_FORCE matcher)
void Regions_Matcher(
  const features::Regions * regions_I,
  const features::Regions * regions_J,
  IndMatches & matches, // photometric corresponding points
  // Distance ratio used to discard spurious correspondence
  const float fDistRatio)
{
  if (regions_I->Type_id() != regions_J->Type_id())
    return;

  if (regions_I->IsScalar())
  {
    if(regions_I->Type_id() == typeid(unsigned char).name())
    {
      // Build on the fly unsigned char based Matcher
      typedef L2_Vectorized<unsigned char> MetricT;
      typedef ArrayMatcherBruteForce<unsigned char, MetricT> MatcherT;
      /// Match the distRatio to the used metric
      Template_Matcher<MatcherT>(regions_I, regions_J, matches, Square(fDistRatio));
    }
    else
    if(regions_I->Type_id() == typeid(float).name())
    {
      typedef L2_Vectorized<float> MetricT;
      typedef ArrayMatcherBruteForce<float, MetricT> MatcherT;
      // Build on the fly float based Matcher
      Template_Matcher<MatcherT>(regions_I, regions_J, matches, Square(fDistRatio));
    }
    else
    if(regions_I->Type_id() == typeid(double).name())
    {
      typedef L2_Vectorized<double> MetricT;
      typedef ArrayMatcherBruteForce<double, MetricT> MatcherT;
      /// Match the distRatio to the used metric
      Template_Matcher<MatcherT>(regions_I, regions_J, matches, Square(fDistRatio));
    }
  }
  else
  if (regions_I->IsBinary() && regions_I->Type_id() == typeid(unsigned char).name())
  {
    typedef Hamming<unsigned char> Metric;
    typedef ArrayMatcherBruteForce<unsigned char, Metric> MatcherT;
    Template_Matcher<MatcherT>(regions_I, regions_J, matches, Square(fDistRatio));
  }
}

} // namespace openMVG
