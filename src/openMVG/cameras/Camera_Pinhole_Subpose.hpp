
// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "ceres/rotation.h"

#include <vector>

namespace openMVG
{
namespace cameras
{

/**
* @brief Define an ideal Pinhole camera intrinsics with a subpose  (used for camera embedded in a rigid rig configuration)
* @note This is an ideal Pinhole camera because it doesn't handle skew and distortion
* @note The camera does only handle one focal length (ie: \f$ f_x = f_y = f \f$ )
* @note Subpose is hidden and not used directly in the class
* @note Use intrinsic->subpose() * globalPose when you send a pose to this class.
*/
class Pinhole_Intrinsic_Subpose : public Pinhole_Intrinsic
{
  typedef Pinhole_Intrinsic_Subpose class_type;

  protected:

    geometry::Pose3 subpose_;

  public:

    /**
    * @brief Constructor
    * @param w Width of the image plane
    * @param h Height of the image plane
    * @param focal_length_pix Focal length (in pixel) of the camera
    * @param ppx Principal point on x-axis
    * @param ppy Principal point on y-axis
    * @param subpose The camera subpose
    */
    Pinhole_Intrinsic_Subpose(
      unsigned int w = 0,
      unsigned int h = 0,
      double focal_length_pix = 0.0,
      double ppx = 0.0,
      double ppy = 0.0,
      geometry::Pose3 subpose = geometry::Pose3())
      : Pinhole_Intrinsic( w, h, focal_length_pix, ppx, ppy),
        subpose_(subpose)
    {
    }

    /**
    * @brief Destructor
    */
    virtual ~Pinhole_Intrinsic_Subpose() override = default;

    /**
    * @brief Get type of the intrinsic
    * @retval PINHOLE_CAMERA_SUBPOSE
    */
    EINTRINSIC getType() const override
    {
      return PINHOLE_CAMERA_SUBPOSE;
    }

    /**
    * @brief Return if the camera use a subpose
    */
    bool use_subpose() const override
    {
      return true;
    }

    /**
    * @brief Return the used subpose
    */
    geometry::Pose3 get_subpose() const override
    {
      return subpose_;
    }

    /**
    * @brief Return the projection matrix (interior & exterior) as a simplified projective projection
    * @param pose Extrinsic matrix
    * @return Concatenation of intrinsic matrix and extrinsic matrix
    * @note use the argument pose (use intrinsic->subpose() * globalPose in order to use full pose)
    */
    Mat34 get_projective_equivalent( const geometry::Pose3 & pose ) const override
    {
      Mat34 P;
      P_From_KRt( K(), pose.rotation(), pose.translation(), &P );
      return P;
    }

    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    std::vector<double> getParams() const override
    {
      std::vector<double> params = Pinhole_Intrinsic::getParams();
      const Mat3 R = subpose_.rotation();
      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      params.push_back(angleAxis[0]);
      params.push_back(angleAxis[1]);
      params.push_back(angleAxis[2]);
      const Vec3 translation = subpose_.translation();
      params.push_back(translation[0]);
      params.push_back(translation[1]);
      params.push_back(translation[2]);
      return params;
    }


    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    bool updateFromParams(const std::vector<double> & params) override
    {
      if (params.size() == 9) {
        const double angleAxis[3] = { params[3], params[4], params[5] }; // rx, ry, rz
        Mat3 R;
        ceres::AngleAxisToRotationMatrix(&angleAxis[0], R.data());
        Vec3 translation;
        translation << params[6], params[7], params[8]; // tx, ty, tz
        *this = Pinhole_Intrinsic_Subpose(
          w_, h_,
          params[0], params[1], params[2], // focal, ppx, ppy
          geometry::Pose3(R, -R.transpose() * translation));
        return true;
      }
      else
      {
        return false;
      }
    }

    /**
    * @brief Return the list of parameter indexes that must be held constant
    * @param parametrization The given parametrization
    */
    std::vector<int> subsetParameterization
    (
      const Intrinsic_Parameter_Type & parametrization) const override
    {
      std::vector<int> constant_index;
      const int param = static_cast<int>(parametrization);
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH)
           || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(0);
      }
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT)
          || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        constant_index.push_back(1);
        constant_index.push_back(2);
      }
      if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_SUBPOSE)
          || param & (int)Intrinsic_Parameter_Type::NONE )
      {
        // Rotation
        constant_index.push_back(3);
        constant_index.push_back(4);
        constant_index.push_back(5);
        // Translation
        constant_index.push_back(6);
        constant_index.push_back(7);
        constant_index.push_back(8);
      }
      return constant_index;
    }

    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    void save( Archive & ar ) const
    {
      Pinhole_Intrinsic::save(ar);
      ar(cereal::make_nvp("subpose", subpose_));
    }


    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar )
    {
      Pinhole_Intrinsic::load(ar);
      ar(cereal::make_nvp("subpose", subpose_));
    }

    /**
    * @brief Clone the object
    * @return A clone (copy of the stored object)
    */
    IntrinsicBase * clone( void ) const override
    {
      return new class_type( *this );
    }
};

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Subpose, "pinhole_subpose" );
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::Pinhole_Intrinsic, openMVG::cameras::Pinhole_Intrinsic_Subpose)

