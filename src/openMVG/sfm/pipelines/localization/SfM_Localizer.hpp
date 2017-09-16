// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_PIPELINES_LOCALIZATION_SFM_LOCALIZER_HPP
#define OPENMVG_SFM_PIPELINES_LOCALIZATION_SFM_LOCALIZER_HPP

#include <limits>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/types.hpp"

namespace openMVG { namespace cameras { struct IntrinsicBase; } }
namespace openMVG { namespace features { class Regions; } }
namespace openMVG { namespace geometry { class Pose3; } }
namespace openMVG { namespace sfm { struct Regions_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

struct Image_Localizer_Match_Data
{
  Mat34 projection_matrix;
  Mat pt3D;
  Mat pt2D;
  std::vector<uint32_t> vec_inliers;
  // Upper bound pixel(s) tolerance for residual errors
  double error_max = std::numeric_limits<double>::infinity();
  uint32_t max_iteration = 4096;
};

class SfM_Localizer
{
public:
  virtual ~SfM_Localizer() = default;

  /**
  * @brief Build the retrieval database (3D points descriptors)
  *
  * @param[in] sfm_data the SfM scene that have to be described
  * @param[in] region_provider regions provider
  * @return True if the database has been correctly setup
  */
  virtual bool Init
  (
    const sfm::SfM_Data & sfm_data,
    const sfm::Regions_Provider & region_provider
  ) = 0;

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
  virtual bool Localize
  (
    const Pair & image_size,
    const cameras::IntrinsicBase * optional_intrinsics,
    const features::Regions & query_regions,
    geometry::Pose3 & pose,
    Image_Localizer_Match_Data * resection_data = nullptr // optional
  ) const = 0;


  /**
  * @brief Try to localize an image from known 2D-3D matches
  *
  * @param[in] image_size the w,h image size
  * @param[in] optional_intrinsics camera intrinsic if known (else nullptr)
  * @param[in,out] resection_data matching data (with filled 2D-3D correspondences)
  * @param[out] pose found pose
  * @return True if a putative pose has been estimated
  */
  static bool Localize
  (
    const Pair & image_size,
    const cameras::IntrinsicBase * optional_intrinsics,
    Image_Localizer_Match_Data & resection_data,
    geometry::Pose3 & pose
  );

  /**
  * @brief Refine a pose according 2D-3D matching & camera model data
  *
  * @param[in,out] intrinsics Camera model
  * @param[in,out] pose Camera pose
  * @param[in] matching_data Corresponding 2D-3D data
  * @param[in] b_refine_pose tell if pose must be refined
  * @param[in] b_refine_intrinsic tell if intrinsics must be refined
  * @return True if the refinement decreased the RMSE pixel residual error
  */
  static bool RefinePose
  (
    cameras::IntrinsicBase * intrinsics,
    geometry::Pose3 & pose,
    Image_Localizer_Match_Data & matching_data,
    bool b_refine_pose,
    bool b_refine_intrinsic
  );
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_PIPELINES_LOCALIZATION_SFM_LOCALIZER_HPP
