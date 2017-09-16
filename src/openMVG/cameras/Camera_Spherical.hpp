// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_SPHERICAL_HPP
#define OPENMVG_CAMERAS_CAMERA_SPHERICAL_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Intrinsics.hpp"

namespace openMVG
{
namespace cameras
{

/**
 * @brief Implement a Spherical camera model
 */
class Intrinsic_Spherical : public IntrinsicBase {

using  class_type = Intrinsic_Spherical;

  public:

  /**
  * @brief Constructor
  * @param w Width of the image plane
  * @param h Height of the image plane
  */
  Intrinsic_Spherical
  (
    unsigned int w = 0,
    unsigned int h = 0
  )
  : IntrinsicBase(w, h)
  {
  }

  ~Intrinsic_Spherical() override = default;

  /**
  * @brief Tell from which type the embed camera is
  * @retval CAMERA_SPHERICAL
  */
  virtual EINTRINSIC getType() const override
  {
    return CAMERA_SPHERICAL;
  }

  /**
  * @brief Data wrapper for non linear optimization (get data)
  * @return an empty vector of parameter since a spherical camera does not have any intrinsic parameter
  */
  virtual std::vector<double> getParams() const override
  {
    return {};
  }

  /**
  * @brief Data wrapper for non linear optimization (update from data)
  * @param params List of params used to update this intrinsic
  * @retval true if update is correct
  * @retval false if there was an error during update
  */
  virtual bool updateFromParams(const std::vector<double> &params) override
  {
    return true;
  }

  /**
  * @brief Return the list of parameter indexes that must be held constant
  * @param parametrization The given parametrization
  */
  virtual std::vector<int> subsetParameterization
  (
    const Intrinsic_Parameter_Type & parametrization
  ) const override
  {
    return {};
  }

  /**
  * @brief Transform a point from the camera plane to the image plane
  * @param p Camera plane point
  * @return Point on image plane
  */
  virtual Vec2 cam2ima(const Vec2 &p) const override
  {
    const size_t size = std::max(w(), h());
    return {
      p.x() * size + w() / 2.0,
      p.y() * size + h() / 2.0 };
  }

  /**
  * @brief Transform a point from the image plane to the camera plane
  * @param p Image plane point
  * @return camera plane point
  */
  virtual Vec2 ima2cam(const Vec2 &p) const override
  {
    const size_t size = std::max(w(), h());
    return {
      (p.x() - w() / 2.0) / size,
      (p.y() - h() / 2.0) / size };
  }

  /**
  * @brief Get bearing vector of a point given an image coordinate
  * @return bearing vector
  */
  virtual Vec3 operator () ( const Vec2& p ) const override
  {
    const Vec2 uv = ima2cam(p);

    const double
      lon = uv.x() * 2 * M_PI,
      lat = uv.y() * 2 * M_PI;

    return {
      cos(lat) * sin(lon),
      -sin(lat),
      cos(lat) * cos(lon)};
  }

  /**
  * @brief Compute projection of a 3D point into the image plane
  * (Apply pose, disto (if any) and Intrinsics)
  * @param pose Pose used to compute projection
  * @param pt3D 3D-point to project on image plane
  * @return Projected (2D) point on image plane
  */
  Vec2 project(
    const geometry::Pose3 & pose,
    const Vec3 & pt3D ) const override
  {
    const Vec3 X = pose( pt3D ); // apply pose
    const double lon = atan2(X.x(), X.z()); // Horizontal normalization of the  X-Z component
    const double lat = atan2(-X.y(), sqrt(X.x()*X.x() + X.z()*X.z())); // Tilt angle
    // denormalization (angle to pixel value)
    return cam2ima({lon / (2 * M_PI), lat / (2 * M_PI)});
  }

  /**
  * @brief Does the camera model handle a distortion field?
  * @retval false
  */
  virtual bool have_disto() const override { return false; }

  /**
  * @brief Add the distortion field to a point (that is in normalized camera frame)
  * @param p Point before distortion computation (in normalized camera frame)
  * @return the initial point p (spherical camera does not have distortion field)
  */
  virtual Vec2 add_disto(const Vec2 &p) const override { return p; }

  /**
  * @brief Remove the distortion to a camera point (that is in normalized camera frame)
  * @param p Point with distortion
  * @return the initial point p (spherical camera does not have distortion field)
  */
  virtual Vec2 remove_disto(const Vec2 &p) const override { return p; }

  /**
  * @brief Return the un-distorted pixel (with removed distortion)
  * @param p Input distorted pixel
  * @return Point without distortion
  */
  virtual Vec2 get_ud_pixel(const Vec2 &p) const override { return p; }

  /**
  * @brief Return the distorted pixel (with added distortion)
  * @param p Input pixel
  * @return Distorted pixel
  */
  virtual Vec2 get_d_pixel(const Vec2 &p) const override { return p; }

  /**
  * @brief Normalize a given unit pixel error to the camera plane
  * @param value Error in image plane
  * @return error of passing from the image plane to the camera plane
  */
  virtual double imagePlane_toCameraPlaneError(double value) const override { return value; }

  /**
  * @brief Return the projection matrix (interior & exterior) as a simplified projective projection
  * @param pose Extrinsic matrix
  * @return Concatenation of intrinsic matrix and extrinsic matrix
  */
  virtual Mat34 get_projective_equivalent(const geometry::Pose3 &pose) const override
  {
    return HStack(pose.rotation(), pose.translation());
  }

  /**
  * @brief Serialization out
  * @param ar Archive
  */
  template <class Archive>
  inline void save( Archive & ar ) const;

  /**
  * @brief  Serialization in
  * @param ar Archive
  */
  template <class Archive>
  inline void load( Archive & ar );

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

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_SPHERICAL_HPP
