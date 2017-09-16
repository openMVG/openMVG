// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_GEOMETRIC_FILTER_UTILS_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_GEOMETRIC_FILTER_UTILS_HPP

#include "openMVG/matching/indMatch.hpp"
#include <openMVG/features/feature_container.hpp>
#include <openMVG/numeric/eigen_alias_definition.hpp>

namespace openMVG { namespace cameras { struct IntrinsicBase; } }
namespace openMVG { namespace sfm { struct Regions_Provider; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }
namespace openMVG { namespace sfm { struct View; } }

namespace openMVG {

namespace matching_image_collection {

/**
* @brief Get "image perfect" features (un-distorted feature positions)
* @param[in] putativeMatches Selected corresponding features id (match)
* @param[in] cam_I Inth Camera interface
* @param[in] feature_I Inth view features
* @param[in] cam_J Jnth Camera interface
* @param[in] feature_J Jnth view features
* @param[out] x_I Pixel perfect features from the Inth image putativeMatches matches
* @param[out] x_J Pixel perfect features from the Jnth image putativeMatches matches
*/
void MatchesPointsToMat
(
  const matching::IndMatches & putativeMatches,
  const cameras::IntrinsicBase * cam_I,
  const features::PointFeatures & feature_I,
  const cameras::IntrinsicBase * cam_J,
  const features::PointFeatures & feature_J,
  Mat & x_I,
  Mat & x_J
);

/**
* @brief Get un-distorted feature positions for the pair pairIndex from the Regions_Provider interface
* @param[in] pairIndex Pair from which you need to extract the corresponding points
* @param[in] putativeMatches Matches of the 'pairIndex' pair
* @param[in] sfm_data SfM_Data scene container
* @param[in] regions_provider Interface that provides the features positions
* @param[out] x_I Pixel perfect features from the Inth image putativeMatches matches
* @param[out] x_J Pixel perfect features from the Jnth image putativeMatches matches
*/
void MatchesPairToMat
(
  const Pair pairIndex,
  const matching::IndMatches & putativeMatches,
  const sfm::SfM_Data * sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  Mat & x_I,
  Mat & x_J
);

/**
* @brief Get un-distorted feature positions for the pair pairIndex from the Features_Provider interface
* @param[in] pairIndex Pair from which you need to extract the corresponding points
* @param[in] putativeMatches Matches of the 'pairIndex' pair
* @param[in] sfm_data SfM_Data scene container
* @param[in] features_provider Interface that provides the features positions
* @param[out] x_I Pixel perfect features from the Inth image putativeMatches matches
* @param[out] x_J Pixel perfect features from the Jnth image putativeMatches matches
*/
void MatchesPairToMat
(
  const Pair pairIndex,
  const matching::IndMatches & putativeMatches,
  const sfm::SfM_Data * sfm_data,
  const std::shared_ptr<sfm::Features_Provider> & features_provider,
  Mat & x_I,
  Mat & x_J
);

} //namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_GEOMETRIC_FILTER_UTILS_HPP
