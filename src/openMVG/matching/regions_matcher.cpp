// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/metric.hpp"
#include "openMVG/matching/metric_hamming.hpp"

namespace openMVG {
namespace matching {

void Match
(
  const matching::EMatcherType & matcher_type,
  const features::Regions & database_regions,
  const features::Regions & query_regions,
  matching::IndMatches & matches
)
{
  const std::unique_ptr<RegionsMatcher> matcher =
    RegionMatcherFactory(matcher_type, database_regions);
  if (matcher)
  {
    matcher->Match(query_regions, matches);
  }
}

void DistanceRatioMatch
(
  float f_dist_ratio,
  const matching::EMatcherType & matcher_type,
  const features::Regions & database_regions,
  const features::Regions & query_regions,
  matching::IndMatches & matches
)
{
  const std::unique_ptr<RegionsMatcher> matcher =
    RegionMatcherFactory(matcher_type, database_regions);
  if (matcher)
  {
    matcher->MatchDistanceRatio(f_dist_ratio, query_regions, matches);
  }
}

std::unique_ptr<RegionsMatcher> RegionMatcherFactory
(
  matching::EMatcherType eMatcherType,
  const features::Regions & regions
)
{
  // Handle invalid request
  if (regions.IsScalar() && eMatcherType == BRUTE_FORCE_HAMMING)
    return {};
  if (regions.IsBinary() && eMatcherType != BRUTE_FORCE_HAMMING)
    return {};

  std::unique_ptr<RegionsMatcher> region_matcher;
  // Switch regions type ID, matcher & Metric: initialize the Matcher interface
  if (regions.IsScalar())
  {
    if (regions.Type_id() == typeid(unsigned char).name())
    {
      // Build on the fly unsigned char based Matcher
      switch (eMatcherType)
      {
        case BRUTE_FORCE_L2:
        {
          using MetricT = L2<unsigned char>;
          using MatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        case ANN_L2:
        {
          using MetricT = flann::L2<unsigned char>;
          using MatcherT = ArrayMatcher_Kdtree_Flann<unsigned char, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        case CASCADE_HASHING_L2:
        {
          using MetricT = L2<unsigned char>;
          using MatcherT = ArrayMatcherCascadeHashing<unsigned char, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        default:
          std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
    else if (regions.Type_id() == typeid(float).name())
    {
      // Build on the fly float based Matcher
      switch (eMatcherType)
      {
        case BRUTE_FORCE_L2:
        {
          using MetricT = L2<float>;
          using MatcherT = ArrayMatcherBruteForce<float, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        case ANN_L2:
        {
          using MetricT = flann::L2<float>;
          using MatcherT = ArrayMatcher_Kdtree_Flann<float, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        case CASCADE_HASHING_L2:
        {
          using MetricT = L2<float>;
          using MatcherT = ArrayMatcherCascadeHashing<float, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        default:
          std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
    else if (regions.Type_id() == typeid(double).name())
    {
      // Build on the fly double based Matcher
      switch (eMatcherType)
      {
        case BRUTE_FORCE_L2:
        {
          using MetricT = L2<double>;
          using MatcherT = ArrayMatcherBruteForce<double, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        case ANN_L2:
        {
          using MetricT = flann::L2<double>;
          using MatcherT = ArrayMatcher_Kdtree_Flann<double, MetricT>;
          region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, true));
        }
        break;
        case CASCADE_HASHING_L2:
        {
          std::cerr << "Not implemented" << std::endl;
        }
        break;
        default:
          std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
  }
  else if (regions.IsBinary() && regions.Type_id() == typeid(unsigned char).name())
  {
    switch (eMatcherType)
    {
      case BRUTE_FORCE_HAMMING:
      {
        using MetricT = Hamming<unsigned char>;
        using MatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
        region_matcher.reset(new matching::RegionsMatcherT<MatcherT>(regions, false));
      }
      break;
      default:
        std::cerr << "Using unknown matcher type" << std::endl;
    }
  }
  else
  {
    std::cerr << "Please consider add this region type_id to RegionMatcherFactory(...)\n"
      << "typeid: " << regions.Type_id() << std::endl;
  }
  return region_matcher;
}

}  // namespace matching
}  // namespace openMVG
