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

void DistanceRatioMatch
(
  float f_dist_ratio,
  matching::EMatcherType eMatcherType,
  const features::Regions & regions_I, // database
  const features::Regions & regions_J, // query
  matching::IndMatches & matches // photometric corresponding points
)
{
  Matcher_Regions_Database matcher(eMatcherType, regions_I);
  matcher.Match(f_dist_ratio, regions_J, matches);
}

bool Matcher_Regions_Database::Match
(
  float dist_ratio, // Distance ratio used to discard spurious correspondence
  const features::Regions & query_regions,
  matching::IndMatches & matches // photometric corresponding points
)const
{
  if (query_regions.RegionCount() == 0 || ! matching_interface_)
  {
    return false;
  }

  matching_interface_->Match(dist_ratio, query_regions, matches);
  return true;
}

Matcher_Regions_Database::Matcher_Regions_Database():
  eMatcherType_(BRUTE_FORCE_L2),
  matching_interface_(nullptr)
{}

Matcher_Regions_Database::Matcher_Regions_Database
(
  matching::EMatcherType eMatcherType,
  const features::Regions & database_regions // database
):
  eMatcherType_(eMatcherType)
{
  // Handle invalid request
  if (database_regions.IsScalar() && eMatcherType == BRUTE_FORCE_HAMMING)
    return;
  if (database_regions.IsBinary() && eMatcherType != BRUTE_FORCE_HAMMING)
    return;

  // Switch regions type ID, matcher & Metric: initialize the Matcher interface
  if (database_regions.IsScalar())
  {
    if (database_regions.Type_id() == typeid(unsigned char).name())
    {
      // Build on the fly unsigned char based Matcher
      switch (eMatcherType_)
      {
        case BRUTE_FORCE_L2:
        {
          using MetricT = L2<unsigned char>;
          using MatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        case ANN_L2:
        {
          using MetricT = flann::L2<unsigned char>;
          using MatcherT = ArrayMatcher_Kdtree_Flann<unsigned char, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        case CASCADE_HASHING_L2:
        {
          using MetricT = L2<unsigned char>;
          using MatcherT = ArrayMatcherCascadeHashing<unsigned char, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        default:
          std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
    else if (database_regions.Type_id() == typeid(float).name())
    {
      // Build on the fly float based Matcher
      switch (eMatcherType)
      {
        case BRUTE_FORCE_L2:
        {
          using MetricT = L2<float>;
          using MatcherT = ArrayMatcherBruteForce<float, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        case ANN_L2:
        {
          using MetricT = flann::L2<float>;
          using MatcherT = ArrayMatcher_Kdtree_Flann<float, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        case CASCADE_HASHING_L2:
        {
          using MetricT = L2<float>;
          using MatcherT = ArrayMatcherCascadeHashing<float, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        default:
          std::cerr << "Using unknown matcher type" << std::endl;
      }
    }
    else if (database_regions.Type_id() == typeid(double).name())
    {
      // Build on the fly double based Matcher
      switch (eMatcherType)
      {
        case BRUTE_FORCE_L2:
        {
          using MetricT = L2<double>;
          using MatcherT = ArrayMatcherBruteForce<double, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
        }
        break;
        case ANN_L2:
        {
          using MetricT = flann::L2<double>;
          using MatcherT = ArrayMatcher_Kdtree_Flann<double, MetricT>;
          matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, true));
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
  else if (database_regions.IsBinary() && database_regions.Type_id() == typeid(unsigned char).name())
  {
    switch (eMatcherType)
    {
      case BRUTE_FORCE_HAMMING:
      {
        using MetricT = Hamming<unsigned char>;
        using MatcherT = ArrayMatcherBruteForce<unsigned char, MetricT>;
        matching_interface_.reset(new matching::RegionsMatcherT<MatcherT>(database_regions, false));
      }
      break;
      default:
          std::cerr << "Using unknown matcher type" << std::endl;
    }
  }
  else
  {
    std::cerr << "Please consider add this region type_id to Matcher_Regions_Database::Match(...)\n"
      << "typeid: " << database_regions.Type_id() << std::endl;
  }
}

}  // namespace matching
}  // namespace openMVG
