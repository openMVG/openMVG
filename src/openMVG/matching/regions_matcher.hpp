// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_REGION_MATCHER_HPP
#define OPENMVG_MATCHING_REGION_MATCHER_HPP

#include <vector>

#include "openMVG/features/regions.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace matching {

/**
 * @brief Match two Regions according a chosen MatcherType.
 * For each query region the nearest entry in the database is kept.
 * @param[in] matcher_type The Matcher type.
 * @param[in] database_regions The database regions.
 * @param[in] query_regions The query regions.
 * @param[out] matches The computed correspondences indices.
 */
void Match
(
  const matching::EMatcherType & matcher_type,
  const features::Regions & database_regions,
  const features::Regions & query_regions,
  matching::IndMatches & matches
);

/**
 * @brief Match two Regions according a chosen MatcherType. A distance ratio
 * between the two nearest neighbor is used to discard some false positive.
 * @param[in] dist_ratio The distance ratio value. Range [0;1].
 *   0.8 or 0.6 are commonly used values.
 * @param[in] matcher_type The Matcher type.
 * @param[in] database_regions The database regions.
 * @param[in] query_regions The query regions.
 * @param[out] matches The computed correspondences indices.
 */
void DistanceRatioMatch
(
  float dist_ratio,
  const matching::EMatcherType & matcher_type,
  const features::Regions & database_regions,
  const features::Regions & query_regions,
  matching::IndMatches & matches
);

class RegionsMatcher
{
  public:
  virtual ~RegionsMatcher() = default;

  /**
   * @brief Match some regions to the database
   * Look for each query to the nearest neighbor.
   */
  virtual bool Match
  (
    const features::Regions & query_regions,
    matching::IndMatches & vec_putative_matches
  ) = 0;

  /**
   * @brief Match some regions to the database
   * Look for each query to the 2 nearest neighbor and keep the match if it pass
   * the distance ratio test: (first_distance < second_distance * f_dist_ratio).
   */
  virtual bool MatchDistanceRatio
  (
    const float dist_ratio,
    const features::Regions & query_regions,
    matching::IndMatches & vec_putative_matches
  ) = 0;
};

/**
 * @brief Create a region matcher according a matcher type and the regions type.
 * @param[in] matcher_type The Matcher type.
 * @param[in] regions The database regions.
 * @return The created RegionsMatcher or an empty smart pointer if the a matcher
 * for the region type asked matcher type cannot be created.
 */
std::unique_ptr<RegionsMatcher> RegionMatcherFactory
(
  matching::EMatcherType matcher_type,
  const features::Regions & regions
);

/**
 * Match two Regions with one stored as a "database" according a Template ArrayMatcher.
 * Template is required in order to make the allocation of the distance array in the good data type.
 */
template < class ArrayMatcherT >
class RegionsMatcherT : public RegionsMatcher
{
private:
  ArrayMatcherT matcher_;
  const features::Regions * regions_;
  const bool b_squared_metric_; // Store if the metric is squared or not
public:
  using Scalar = typename ArrayMatcherT::ScalarT;
  using DistanceType = typename ArrayMatcherT::DistanceType;

  RegionsMatcherT(): regions_(nullptr), b_squared_metric_(false) {}

  /**
   * @brief Init the matcher with some reference regions.
   */
  RegionsMatcherT
  (
    const features::Regions & regions,
    bool b_squared_metric = false
  ):
    regions_(&regions),
    b_squared_metric_(b_squared_metric)
  {
    if (regions_->RegionCount() == 0)
      return;

    const Scalar * tab = reinterpret_cast<const Scalar *>(regions_->DescriptorRawData());
    matcher_.Build(tab, regions_->RegionCount(), regions_->DescriptorLength());
  }

  bool Match
  (
    const features::Regions & query_regions,
    matching::IndMatches & matches
  ) override
  {
    if (!regions_)
      return false;

    const Scalar * queries = reinterpret_cast<const Scalar *>(query_regions.DescriptorRawData());

    // Search the closest neighbour for each query descriptor
    std::vector<DistanceType> distances;
    matches.clear();
    if (!matcher_.SearchNeighbours(queries, query_regions.RegionCount(), &matches, &distances, 1))
      return false;

    for ( auto & match : matches )
    {
      std::swap(match.i_, match.j_);
    }

    return (!matches.empty());
  }

  /**
   * @brief Match some regions to the database of internal regions.
   */
  bool MatchDistanceRatio
  (
    const float distance_ratio,
    const features::Regions & query_regions,
    matching::IndMatches & matches
  ) override
  {
    if (!regions_)
      return false;

    const Scalar * queries = reinterpret_cast<const Scalar *>(query_regions.DescriptorRawData());

    const size_t number_neighbor = 2;
    matching::IndMatches nn_matches;
    std::vector<DistanceType> nn_distances;

    // Search the 2 closest neighbours for each query descriptor
    if (!matcher_.SearchNeighbours(queries,
                                   query_regions.RegionCount(),
                                   &nn_matches,
                                   &nn_distances,
                                   number_neighbor))
      return false;

    std::vector<int> nn_ratio_indexes;
    // Filter the matches using a distance ratio test:
    //   The probability that a match is correct is determined by taking
    //   the ratio of distance from the closest neighbor to the distance
    //   of the second closest.
    matching::NNdistanceRatio(
      nn_distances.cbegin(), // distance start
      nn_distances.cend(),   // distance end
      number_neighbor,       // Number of neighbor in iterator sequence (minimum required 2)
      nn_ratio_indexes,      // output (indices that respect the distance Ratio)
      b_squared_metric_ ? Square(distance_ratio) : distance_ratio);

    matches.clear();
    matches.reserve(nn_ratio_indexes.size());
    for (const auto & index : nn_ratio_indexes)
    {
      matches.emplace_back(nn_matches[index * number_neighbor].j_,
                           nn_matches[index * number_neighbor].i_);
    }

    return (!matches.empty());
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_REGION_MATCHER_HPP
