// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_MATCHER_REGIONS_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_MATCHER_REGIONS_HPP

#include <memory>

#include "openMVG/matching/matcher_type.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

namespace openMVG { namespace matching { class PairWiseMatchesContainer; } }
namespace openMVG { namespace sfm { struct Regions_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace matching_image_collection {

/// Implementation of an Image Collection Matcher
/// Compute putative matches between a collection of pictures
/// Spurious correspondences are discarded by using the
///  a threshold over the distance ratio of the 2 nearest neighbours.
///
class Matcher_Regions : public Matcher
{
  public:
  Matcher_Regions
  (
    float dist_ratio,
    matching::EMatcherType eMatcherType
  );

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
  // Distance ratio used to discard spurious correspondence
  float f_dist_ratio_;
  // Matcher Type
  matching::EMatcherType eMatcherType_;
};

} // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_MATCHER_REGIONS_HPP
