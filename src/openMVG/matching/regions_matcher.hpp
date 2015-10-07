
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"

#include "openMVG/numeric/numeric.h"
#include "openMVG/features/regions.hpp"

#include <vector>

namespace openMVG {
namespace matching {

/**
 * @brief Match two Regions according a chosen MatcherType using the ratio test
 * to assure the robustness of the matches.
 * It returns the matching features.
 * 
 * @param[in] dist_ratio The threshold for the ratio test.
 * @param[in] eMatcherType The type of matcher to use.
 * @param[in] regions_I The first Region to match
 * @param[in] regions_J The second Region to match
 * @param[out] matches It contains the indices of the matching features
 */
void DistanceRatioMatch
(
  float dist_ratio,   // Distance ratio
  matching::EMatcherType eMatcherType, // Matcher
  const features::Regions & regions_I, // database
  const features::Regions & regions_J, // query
  matching::IndMatches & matches // corresponding points
);


/**
 * @brief Interface for a matcher of Regions that keep on of the Regions stored
 * as a "database". The other region is then a query that is matched with respect
 * to the database.
 */
class RegionsMatcher
{
  public:
    
   /**
    * @brief The destructor.
    */ 
   ~RegionsMatcher() {}


    /**
     * @brief Initialize the matcher by setting one region as a "database".
     * 
     * @param[in] regions The Regions to be used as reference when matching another Region.
     */
    virtual void Init_database
    (
      const features::Regions& regions
    ) = 0;

    /**
     * @brief Match a Regions to the internal database using the test ratio to improve
     * the robustness of the match.
     * 
     * @param[in] f_dist_ratio The threshold for the ratio test.
     * @param[in] query_regions The Regions to match.
     * @param[out] vec_putative_matches It contains the indices of the matching features
     * of the database and the query Regions.
     * @return True if everything went well.
     */
    virtual bool Match
    (
      const float f_dist_ratio,
      const features::Regions& query_regions,
      matching::IndMatches & vec_putative_matches
    ) =0;
};

/**
 * @brief A simple class that allows to match Regions with respect a given Regions
 * used as "database".
 */
class Matcher_Regions_Database
{
  public:

    /**
     * @brief Empty constructor, by default it initializes the matcher to BRUTE_FORCE_L2
     * and the database to an empty database.
     */
    Matcher_Regions_Database();

    /**
     * @brief Initialize the internal database
     * 
     * @param[in] eMatcherType The type of matcher to use to match the Regions.
     * @param[in] database_regions The Regions that will be used as database to 
     * match other Regions (query).
     */
    Matcher_Regions_Database
    (
      matching::EMatcherType eMatcherType,
      const features::Regions & database_regions // database
    );

    /**
     * @brief Find corresponding points between the query Regions and the database one
     * 
     * @param[in] dist_ratio The threshold for the ratio test.
     * @param[in] query_regions The Regions to match.
     * @param[out] matches It contains the indices of the matching features
     * of the database and the query Regions.
     * @return True if everything went well.
     */
    bool Match
    (
    float dist_ratio, // Distance ratio used to discard spurious correspondence
      const features::Regions & query_regions,
    matching::IndMatches & matches // photometric corresponding points
    )const;

  private:
  // Matcher Type
  matching::EMatcherType _eMatcherType;
  // The matching interface
  std::unique_ptr<RegionsMatcher> _matching_interface;
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
  typedef typename ArrayMatcherT::ScalarT Scalar;
  typedef typename ArrayMatcherT::DistanceType DistanceType;

  /**
   * @brief Empty constructor, by default it initializes the database to an empty database.
   */
  RegionsMatcherT() : regions_(nullptr) {}

  /**
   * @brief Initialize the matcher with a Regions that will be used as database
   * 
   * @param regions The Regions to be used as database.
   * @param b_squared_metric Whether to use a squared metric for the ratio test 
   * when matching two Regions.
   */
  RegionsMatcherT(const features::Regions& regions, bool b_squared_metric = false)
    : regions_(&regions), b_squared_metric_(b_squared_metric)
  {
    if (regions_->RegionCount() == 0)
      return;

    const Scalar * tab = reinterpret_cast<const Scalar *>(regions_->DescriptorRawData());
    matcher_.Build(tab, regions_->RegionCount(), regions_->DescriptorLength());
  }

  /**
   * @brief Initialize the matcher with a Regions that will be used as database
   * 
   * @param regions The Regions to be used as database.
   */
  void Init_database
  (
    const features::Regions& regions
  )
  {
    regions_ = &regions;
    if (regions_->RegionCount() == 0)
      return;

    const Scalar * tab = reinterpret_cast<const Scalar *>(regions_->DescriptorRawData());
    matcher_.Build(tab, regions_->RegionCount(), regions_->DescriptorLength());
  }

  /**
   * @brief Match a Regions to the internal database using the test ratio to improve
   * the robustness of the match.
   * 
   * @param f_dist_ratio The threshold for the ratio test.
   * @param query_regions The Regions to match.
   * @param vec_putative_matches It contains the indices of the matching features
     * of the database and the query Regions.
   * @return True if everything went well.
   */
  bool Match(
    const float f_dist_ratio,
    const features::Regions& queryregions_,
    matching::IndMatches & vec_putative_matches)
  {
    if (regions_ == nullptr)
      return false;

    const Scalar * queries = reinterpret_cast<const Scalar *>(queryregions_.DescriptorRawData());

    const size_t NNN__ = 2;
    matching::IndMatches vec_nIndice;
    std::vector<DistanceType> vec_fDistance;

    // Search the 2 closest features neighbours for each query descriptor
    if (!matcher_.SearchNeighbours(queries, queryregions_.RegionCount(), &vec_nIndice, &vec_fDistance, NNN__))
      return false;

    std::vector<int> vec_nn_ratio_idx;
    // Filter the matches using a distance ratio test:
    //   The probability that a match is correct is determined by taking
    //   the ratio of distance from the closest neighbor to the distance
    //   of the second closest.
    matching::NNdistanceRatio(
      vec_fDistance.begin(), // distance start
      vec_fDistance.end(),   // distance end
      NNN__, // Number of neighbor in iterator sequence (minimum required 2)
      vec_nn_ratio_idx, // output (indices that respect the distance Ratio)
      b_squared_metric_ ? Square(f_dist_ratio) : f_dist_ratio);

    vec_putative_matches.reserve(vec_nn_ratio_idx.size());
    for (size_t k=0; k < vec_nn_ratio_idx.size(); ++k)
    {
      const size_t index = vec_nn_ratio_idx[k];
      vec_putative_matches.emplace_back(vec_nIndice[index*NNN__]._j, vec_nIndice[index*NNN__]._i);
    }

    // Remove duplicates
    matching::IndMatch::getDeduplicated(vec_putative_matches);

    // Remove matches that have the same (X,Y) coordinates
    matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
      regions_->GetRegionsPositions(), queryregions_.GetRegionsPositions());
    matchDeduplicator.getDeduplicated(vec_putative_matches);

    return (!vec_putative_matches.empty());
  }
  
  const features::Regions* getDatabaseRegions() const { return regions_; } 
  
};

}  // namespace matching
}  // namespace openMVG
