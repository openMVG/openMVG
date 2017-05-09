// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_VIEW_PRIORS_HPP
#define OPENMVG_SFM_SFM_VIEW_PRIORS_HPP

#include <string>
#include <vector>

#include "openMVG/geometry/pose3.hpp"

namespace openMVG {
namespace sfm {

struct View;

using namespace openMVG::geometry;

/**
* @brief Define a View that contains optional Pose priors (pose and/or rotation)
*/
struct ViewPriors : public View
{
  ViewPriors
  (
    const std::string & sImgPath = "",
    IndexT view_id = UndefinedIndexT,
    IndexT intrinsic_id = UndefinedIndexT,
    IndexT pose_id = UndefinedIndexT,
    IndexT width = UndefinedIndexT,
    IndexT height = UndefinedIndexT
  )
  : View
    (
      sImgPath,
      view_id,
      intrinsic_id,
      pose_id,
      width,
      height
    ),
    b_use_pose_center_(false),
    b_use_pose_rotation_(false)
  {
  }

  ~ViewPriors() override = default;

  void SetPoseCenterPrior
  (
    const Vec3 & center,
    const Vec3 & weight
  )
  {
    b_use_pose_center_  = true;
    center_weight_      = weight;
    pose_center_        = center;
  }

  void SetPoseRotationPrior
  (
    const Mat3 & rotation,
    const double weight
  )
  {
    rotation_weight_ = weight;
    pose_rotation_   = rotation;
  }

  /**
  * Serialization out
  * @param ar Archive
  */
  template <class Archive>
  void save( Archive & ar ) const;

  /**
  * @brief Serialization in
  * @param ar Archive
  */
  template <class Archive>
  void load( Archive & ar );

  // Pose center prior
  bool b_use_pose_center_ = false; // Tell if the pose prior must be used
  Vec3 center_weight_ = Vec3::Constant(1.0);
  Vec3 pose_center_ = Vec3::Zero();

  // Pose rotation prior
  bool b_use_pose_rotation_ = false; // Tell if the rotation prior must be used
  double rotation_weight_ = 1.0;
  Mat3 pose_rotation_ = Mat3::Identity();
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_VIEW_PRIORS_HPP
