// Copyright (c) 2016 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_VIEW_PRIORS_HPP
#define OPENMVG_SFM_SFM_VIEW_PRIORS_HPP

#include "openMVG/geometry/pose3.hpp"

#include <cereal/types/polymorphic.hpp>

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
  void save( Archive & ar ) const
  {
    View::save(ar);

    // Pose center prior
    if (b_use_pose_center_)
    {
      ar( cereal::make_nvp( "use_pose_center_prior", b_use_pose_center_ ) );
      const std::vector<double> vec_weights = { center_weight_( 0 ), center_weight_( 1 ), center_weight_( 2 ) };
      ar( cereal::make_nvp( "center_weight", vec_weights ) );
      const std::vector<double> vec = { pose_center_( 0 ), pose_center_( 1 ), pose_center_( 2 ) };
      ar( cereal::make_nvp( "center", vec ) );
    }

    // Pose rotation prior
    /*
    if (b_use_pose_rotation_)
    {
      ar( cereal::make_nvp( "use_pose_rotation_prior", b_use_pose_rotation_ ) );
      ar( cereal::make_nvp( "rotation_weight", rotation_weight_ ) );
      const std::vector<std::vector<double>> mat =
      {
        { pose_rotation_( 0, 0 ), pose_rotation_( 0, 1 ), pose_rotation_( 0, 2 ) },
        { pose_rotation_( 1, 0 ), pose_rotation_( 1, 1 ), pose_rotation_( 1, 2 ) },
        { pose_rotation_( 2, 0 ), pose_rotation_( 2, 1 ), pose_rotation_( 2, 2 ) }
      };
      ar( cereal::make_nvp( "rotation", mat ) );
    }
    */
  }

  /**
  * @brief Serialization in
  * @param ar Archive
  */
  template <class Archive>
  void load( Archive & ar )
  {
    View::load(ar);

    // Pose center prior
    try
    {
      ar( cereal::make_nvp( "use_pose_center_prior", b_use_pose_center_ ) );
      std::vector<double> vec( 3 );
      ar( cereal::make_nvp( "center_weight", vec ) );
      center_weight_ = Eigen::Map<const Vec3>( &vec[0] );
      ar( cereal::make_nvp( "center", vec ) );
      pose_center_ = Eigen::Map<const Vec3>( &vec[0] );
    }
    catch( cereal::Exception e )
    {
      // if it fails just use a default settings
      b_use_pose_center_ = false;
    }

    // Pose rotation prior
    /*
    try
    {
      ar( cereal::make_nvp( "use_pose_rotation_prior", b_use_pose_rotation_ ) );
      ar( cereal::make_nvp( "rotation_weight", rotation_weight_ ) );
      std::vector<std::vector<double>> mat( 3, std::vector<double>( 3 ) );
      ar( cereal::make_nvp( "rotation", mat ) );
      // copy back to the rotation
      pose_rotation_.row( 0 ) = Eigen::Map<const Vec3>( &( mat[0][0] ) );
      pose_rotation_.row( 1 ) = Eigen::Map<const Vec3>( &( mat[1][0] ) );
      pose_rotation_.row( 2 ) = Eigen::Map<const Vec3>( &( mat[2][0] ) );
    }
    catch( const cereal::Exception & e )
    {
      // if it fails just use a default settings
      b_use_pose_rotation_ = false;
    }
    */
  }

  // Pose center prior
  bool b_use_pose_center_ = false; // Tell if the pose prior must be used
  Vec3 center_weight_ = Vec3(1.0,1.0,1.0);
  Vec3 pose_center_ = Vec3::Zero();

  // Pose rotation prior
  bool b_use_pose_rotation_ = false; // Tell if the rotation prior must be used
  double rotation_weight_ = 1.0;
  Mat3 pose_rotation_ = Mat3::Identity();
};

} // namespace sfm
} // namespace openMVG

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::sfm::ViewPriors, "view_priors" );
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::sfm::View, openMVG::sfm::ViewPriors);

#endif // OPENMVG_SFM_SFM_VIEW_PRIORS_HPP
