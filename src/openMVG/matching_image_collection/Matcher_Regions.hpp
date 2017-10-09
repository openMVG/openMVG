// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_MATCHER_REGIONS_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_MATCHER_REGIONS_HPP

#include <functional>
#include <memory>

#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/matching/regions_matcher.hpp"

namespace openMVG { namespace matching { class PairWiseMatchesContainer; } }
namespace openMVG { namespace sfm { struct Regions_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace matching_image_collection {

/// Implementation of an Image Collection Matcher
/// Computes the putative matches between a collection of picture regions.
///
class Matcher_Regions : public Matcher
{
  public:
  /// Functor to create on the fly the required matching::RegionsMatcher
  using CreateMatcherFunctor = std::function
    <std::unique_ptr<matching::RegionsMatcher>(const features::Regions &)>;
  /// Functor to call the Matching function of the created matching::RegionsMatcher
  using CallMatcherFunctor = std::function
    <void(matching::RegionsMatcher *,
          const features::Regions &,
          matching::IndMatches &)>;
  /// Functor to post_process the found matched (can be empty or do features depuplication...)
  using CallPostProcessMatchFunctor = std::function
    <const matching::IndMatches(
      const features::Regions &,
      const features::Regions &,
      const matching::IndMatches &,
      const std::pair<int, int> &,
      const std::pair<int, int> &)>;

  Matcher_Regions
  (
    const CreateMatcherFunctor & create_matcher,
    const CallMatcherFunctor & call_matcher,
    const CallPostProcessMatchFunctor & post_process_match
  ): create_matcher_(create_matcher),
     call_matcher_(call_matcher),
     call_post_process_match_functor(post_process_match)
  {};

  /// Find corresponding points between some pair of view Ids
  void Match
  (
    const sfm::SfM_Data & sfm_data,
    const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
    const Pair_Set & pairs,
    matching::PairWiseMatchesContainer & map_PutativesMatches, // the pairwise photometric corresponding points
    C_Progress *  progress = nullptr
  ) const override;

  private:
  const CreateMatcherFunctor create_matcher_;
  const CallMatcherFunctor call_matcher_;
  const CallPostProcessMatchFunctor call_post_process_match_functor;

};

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_MATCHER_REGIONS_HPP
