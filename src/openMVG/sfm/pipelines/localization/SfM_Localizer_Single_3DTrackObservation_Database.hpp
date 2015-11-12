
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"
#include "openMVG/matching/regions_matcher.hpp"

namespace openMVG {
namespace sfm {

// Implementation of a naive method:
// - init the database of descriptor from the structure and the observations.
// - create a large array with all the used descriptors and init a Matcher with it
// - to localize an input image compare it's regions to the database and robust estimate
//   the pose from found 2d-3D correspondences

class SfM_Localization_Single_3DTrackObservation_Database : public SfM_Localizer
{
public:

  SfM_Localization_Single_3DTrackObservation_Database();

  /**
  * @brief Build the retrieval database (3D points descriptors)
  *
  * @param[in] sfm_data the SfM scene that have to be described
  * @param[in] region_provider regions provider
  * @return True if the database has been correctly setup
  */
  bool Init
  (
    const SfM_Data & sfm_data,
    const Regions_Provider & regions_provider
  );

  /**
  * @brief Try to localize an image in the database
  *
  * @param[in] image_size the w,h image size
  * @param[in] optional_intrinsics camera intrinsic if known (else nullptr)
  * @param[in] query_regions the image regions (type must be the same as the database)
  * @param[out] pose found pose
  * @param[out] resection_data matching data (2D-3D and inliers; optional)
  * @return True if a putative pose has been estimated
  */
  bool Localize
  (
    const Pair & image_size,
    const cameras::IntrinsicBase * optional_intrinsics,
    const features::Regions & query_regions,
    geometry::Pose3 & pose,
    Image_Localizer_Match_Data * resection_data_ptr = NULL // optional
  ) const;

private:
  // Reference to the scene
  const SfM_Data * sfm_data_;
  /// Association of a regions to a landmark observation
  std::unique_ptr<features::Regions> landmark_observations_descriptors_;
  /// Association of a track observation to a track Id (used for retrieval)
  std::vector<IndexT> index_to_landmark_id_;
  /// A matching interface to find matches between 2D descriptor matches
  ///  and 3D points observation descriptors
  std::unique_ptr<matching::Matcher_Regions_Database> matching_interface_;
};

} // namespace sfm
} // namespace openMVG
