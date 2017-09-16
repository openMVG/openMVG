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
 * Match two Regions according a chosen MatcherType.
 */
void DistanceRatioMatch
(
  float dist_ratio,   // Distance ratio
  matching::EMatcherType eMatcherType, // Matcher
  const features::Regions & regions_I, // database
  const features::Regions & regions_J, // query
  matching::IndMatches & matches // corresponding points
);

class RegionsMatcher
{
  public:
  virtual ~RegionsMatcher() = default;

  /**
   * @brief Initialize the retrieval database
   */
  virtual void Init_database
  (
    const features::Regions& regions
  ) = 0;

  /**
   * @brief Match some regions to the database
   */
  virtual bool Match
  (
    const float f_dist_ratio,
    const features::Regions& query_regions,
    matching::IndMatches & vec_putative_matches
  ) =0;
};

/**
 * Regions matching in query to database mode
 */
class Matcher_Regions_Database
{
  public:

  Matcher_Regions_Database();

  /**
   * @brief Initialize the retrieval database
   */
  Matcher_Regions_Database
  (
    matching::EMatcherType eMatcherType,
    const features::Regions & database_regions // database
  );

  /// Find corresponding points between the query regions and the database one
  bool Match
  (
    float dist_ratio, // Distance ratio used to discard spurious correspondence
    const features::Regions & query_regions,
    matching::IndMatches & matches // photometric corresponding points
  )const;

  private:
  // Matcher Type
  matching::EMatcherType eMatcherType_;
  // The matching interface
  std::unique_ptr<RegionsMatcher> matching_interface_;
};

/**
 * Match two Regions with one stored as a "database" according a Template ArrayMatcher.
 */
template < class ArrayMatcherT >
class RegionsMatcherT : public RegionsMatcher
{
private:
  ArrayMatcherT matcher_;
  const features::Regions* regions_;
  bool b_squared_metric_; // Store if the metric is squared or not
public:
  using Scalar = typename ArrayMatcherT::ScalarT;
  using DistanceType = typename ArrayMatcherT::DistanceType;

  RegionsMatcherT() :regions_(nullptr), b_squared_metric_(false) {}

  /**
   * @brief Init the matcher with some reference regions.
   */
  RegionsMatcherT(const features::Regions& regions, bool b_squared_metric = false)
    : regions_(&regions), b_squared_metric_(b_squared_metric)
  {
    if (regions_->RegionCount() == 0)
      return;

    const Scalar * tab = reinterpret_cast<const Scalar *>(regions_->DescriptorRawData());
    matcher_.Build(tab, regions_->RegionCount(), regions_->DescriptorLength());
  }

  void Init_database
  (
    const features::Regions& regions
  ) override
  {
    regions_ = &regions;
    if (regions_->RegionCount() == 0)
      return;

    const Scalar * tab = reinterpret_cast<const Scalar *>(regions_->DescriptorRawData());
    matcher_.Build(tab, regions_->RegionCount(), regions_->DescriptorLength());
  }

  /**
   * @brief Match some regions to the database of internal regions.
   */
  bool Match
  (
    const float f_dist_ratio,
    const features::Regions& queryregions_,
    matching::IndMatches & vec_putative_matches
  ) override
  {
    if (regions_ == nullptr)
      return false;

    const Scalar * queries = reinterpret_cast<const Scalar *>(queryregions_.DescriptorRawData());

    const size_t NNN__ = 2;
    matching::IndMatches vec_Indice;
    std::vector<DistanceType> vec_Distance;

    // Search the 2 closest features neighbours for each query descriptor
    if (!matcher_.SearchNeighbours(queries, queryregions_.RegionCount(), &vec_Indice, &vec_Distance, NNN__))
      return false;

    std::vector<int> vec_nn_ratio_idx;
    // Filter the matches using a distance ratio test:
    //   The probability that a match is correct is determined by taking
    //   the ratio of distance from the closest neighbor to the distance
    //   of the second closest.
    matching::NNdistanceRatio(
      vec_Distance.begin(), // distance start
      vec_Distance.end(),   // distance end
      NNN__, // Number of neighbor in iterator sequence (minimum required 2)
      vec_nn_ratio_idx, // output (indices that respect the distance Ratio)
      b_squared_metric_ ? Square(f_dist_ratio) : f_dist_ratio);

    vec_putative_matches.reserve(vec_nn_ratio_idx.size());
    for ( const auto & index : vec_nn_ratio_idx )
    {
      vec_putative_matches.emplace_back(vec_Indice[index*NNN__].j_, vec_Indice[index*NNN__].i_);
    }

    // Remove duplicates
    matching::IndMatch::getDeduplicated(vec_putative_matches);

    // Remove matches that have the same (X,Y) coordinates
    matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
      regions_->GetRegionsPositions(), queryregions_.GetRegionsPositions());
    matchDeduplicator.getDeduplicated(vec_putative_matches);

    return (!vec_putative_matches.empty());
  }
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_REGION_MATCHER_HPP
