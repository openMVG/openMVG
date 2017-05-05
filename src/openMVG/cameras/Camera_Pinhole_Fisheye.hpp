// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Romain Janvier <romain.janvier~AT~univ-orleans.fr> for the given adaptation

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This is an adaptation of the Fisheye distortion model implemented in OpenCV
// https://github.com/Itseez/opencv/blob/master/modules/calib3d/src/fisheye.cpp

#ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_FISHEYE_HPP
#define OPENMVG_CAMERAS_CAMERA_PINHOLE_FISHEYE_HPP

#include "openMVG/cameras/Camera_Pinhole.hpp"

namespace openMVG
{
namespace cameras
{

/**
* @brief Implement a simple Fish-eye camera model
*/
class Pinhole_Intrinsic_Fisheye : public Pinhole_Intrinsic
{
  using class_type = Pinhole_Intrinsic_Fisheye;

  protected:

    /// center of distortion is applied by the Intrinsics class
    std::vector<double> params_; // K1, K2, K3, K4


  public:

    /**
    * @brief Constructor
    * @param w Width of image plane
    * @param h Height of image plane
    * @param focal Focal distance in pixel
    * @param ppx Principal point on X-axis
    * @param ppy Principal point on Y-axis
    * @param k1 Distortion coefficient
    * @param k2 Distortion coefficient
    * @param k3 Distortion coefficient
    * @param k4 Distortion coefficient
    */
    Pinhole_Intrinsic_Fisheye(
      int w = 0, int h = 0,
      double focal = 0.0, double ppx = 0, double ppy = 0,
      double k1 = 0.0, double k2 = 0.0, double k3 = 0.0, double k4 = 0.0 );

    ~Pinhole_Intrinsic_Fisheye() override = default;

    /**
    * @brief Tell from which type the embed camera is
    * @retval PINHOLE_CAMERA_FISHEYE
    */
    EINTRINSIC getType() const override;

    /**
    * @brief Does the camera model handle a distortion field?
    * @retval true
    */
    bool have_disto() const override;

    /**
    * @brief Add the distortion field to a point (that is in normalized camera frame)
    * @param p Point before distortion computation (in normalized camera frame)
    * @return point with distortion
    */
    Vec2 add_disto( const Vec2 & p ) const override;

    /**
    * @brief Remove the distortion to a camera point (that is in normalized camera frame)
    * @param p Point with distortion
    * @return Point without distortion
    */
    Vec2 remove_disto( const Vec2 & p ) const override;

    /**
    * @brief Data wrapper for non linear optimization (get data)
    * @return vector of parameter of this intrinsic
    */
    std::vector<double> getParams() const override;

    /**
    * @brief Data wrapper for non linear optimization (update from data)
    * @param params List of params used to update this intrinsic
    * @retval true if update is correct
    * @retval false if there was an error during update
    */
    bool updateFromParams( const std::vector<double> & params ) override;

    /**
    * @brief Return the list of parameter indexes that must be held constant
    * @param parametrization The given parametrization
    */
    std::vector<int> subsetParameterization
    (
      const Intrinsic_Parameter_Type & parametrization) const override;

    /**
    * @brief Return the un-distorted pixel (with removed distortion)
    * @param p Input distorted pixel
    * @return Point without distortion
    */
    Vec2 get_ud_pixel( const Vec2& p ) const override;

    /**
    * @brief Return the distorted pixel (with added distortion)
    * @param p Input pixel
    * @return Distorted pixel
    */
    Vec2 get_d_pixel( const Vec2& p ) const override;

    /**
    * @brief Serialization out
    * @param ar Archive
    */
    template <class Archive>
    void save( Archive & ar ) const;

    /**
    * @brief  Serialization in
    * @param ar Archive
    */
    template <class Archive>
    void load( Archive & ar );

    /**
    * @brief Clone the object
    * @return A clone (copy of the stored object)
    */
    IntrinsicBase * clone( void ) const override;
};


} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_PINHOLE_FISHEYE_HPP
