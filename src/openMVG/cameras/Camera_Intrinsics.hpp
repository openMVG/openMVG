// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_INTRINSICS_HPP
#define OPENMVG_CAMERAS_CAMERA_INTRINSICS_HPP

#include <vector>

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/stl/hash.hpp"

namespace openMVG
{
namespace cameras
{

/**
* @brief Struct used to force "clonability"
*/
template< typename T>
struct Clonable
{
  virtual T * clone() const = 0;
};

/**
* @brief Base class used to store common intrinsics parameters
*/
struct IntrinsicBase : public Clonable<IntrinsicBase>
{
  /// Width of image
  unsigned int w_;
  /// Height of image
  unsigned int h_;

  /**
  * @brief Constructor
  * @param w Width of the image
  * @param h Height of the image
  */
  IntrinsicBase( unsigned int w = 0, unsigned int h = 0 )
    : w_( w ),
      h_( h )
  {

  }

  /**
  * @brief Destructor
  */
  virtual ~IntrinsicBase() = default;

  /**
  * @brief Get width of the image
  * @return width of the image
  */
  unsigned int w() const
  {
    return w_;
  }

  /**
  * @brief Get height of the image
  * @return height of the image
  */
  unsigned int h() const
  {
    return h_;
  }

  /**
  * @brief Compute projection of a 3D point into the image plane
  * (Apply pose, disto (if any) and Intrinsics)
  * @param pose Pose used to compute projection
  * @param pt3D 3D-point to project on image plane
  * @return Projected (2D) point on image plane
  */
  virtual Vec2 project(
    const geometry::Pose3 & pose,
    const Vec3 & pt3D ) const
  {
    const Vec3 X = pose( pt3D ); // apply pose
    if ( this->have_disto() ) // apply disto & intrinsics
    {
      return this->cam2ima( this->add_disto( X.hnormalized() ) );
    }
    else // apply intrinsics
    {
      return this->cam2ima( X.hnormalized() );
    }
  }

  /**
  * @brief Compute the residual between the 3D projected point and an image observation
  * @param pose Pose used to project point on camera plane
  * @param X 3d point to project on camera plane
  * @param x image observation
  * @brief Relative 2d distance between projected and observed points
  */
  Vec2 residual(
    const geometry::Pose3 & pose,
    const Vec3 & X,
    const Vec2 & x ) const
  {
    const Vec2 proj = this->project( pose, X );
    return x - proj;
  }

  // --
  // Virtual members
  // --

  /**
  * @brief Tell from which type the embed camera is
  * @return Corresponding intrinsic
  */
  virtual EINTRINSIC getType() const = 0;

  /**
  * @brief Data wrapper for non linear optimization (get data)
  * @return vector of parameter of this intrinsic
  */
  virtual std::vector<double> getParams() const = 0;

  /**
  * @brief Data wrapper for non linear optimization (update from data)
  * @param params List of params used to update this intrinsic
  * @retval true if update is correct
  * @retval false if there was an error during update
  */
  virtual bool updateFromParams( const std::vector<double> & params ) = 0;

  /**
  * @brief Return the list of parameter indexes that must be held constant
  * @param parametrization The given parametrization
  */
  virtual std::vector<int> subsetParameterization(
    const Intrinsic_Parameter_Type & parametrization) const = 0;

  /**
  * @brief Get bearing vector of a point given an image coordinate
  * @return bearing vector
  */
  virtual Vec3 operator () ( const Vec2& p ) const = 0;

  /**
  * @brief Transform a point from the camera plane to the image plane
  * @param p Camera plane point
  * @return Point on image plane
  */
  virtual Vec2 cam2ima( const Vec2& p ) const = 0;

  /**
  * @brief Transform a point from the image plane to the camera plane
  * @param p Image plane point
  * @return camera plane point
  */
  virtual Vec2 ima2cam( const Vec2& p ) const = 0;

  /**
  * @brief Does the camera model handle a distortion field?
  * @retval true if intrinsic holds distortion
  * @retval false if intrinsic does not hold distortion
  */
  virtual bool have_disto() const
  {
    return false;
  }

  /**
  * @brief Add the distortion field to a point (that is in normalized camera frame)
  * @param p Point before distortion computation (in normalized camera frame)
  * @return point with distortion
  */
  virtual Vec2 add_disto( const Vec2& p ) const = 0;

  /**
  * @brief Remove the distortion to a camera point (that is in normalized camera frame)
  * @param p Point with distortion
  * @return Point without distortion
  */
  virtual Vec2 remove_disto( const Vec2& p ) const  = 0;

  /**
  * @brief Return the un-distorted pixel (with removed distortion)
  * @param p Input distorted pixel
  * @return Point without distortion
  */
  virtual Vec2 get_ud_pixel( const Vec2& p ) const = 0;

  /**
  * @brief Return the distorted pixel (with added distortion)
  * @param p Input pixel
  * @return Distorted pixel
  */
  virtual Vec2 get_d_pixel( const Vec2& p ) const = 0;

  /**
  * @brief Normalize a given unit pixel error to the camera plane
  * @param value Error in image plane
  * @return error of passing from the image plane to the camera plane
  */
  virtual double imagePlane_toCameraPlaneError( double value ) const = 0;

  /**
  * @brief Return the projection matrix (interior & exterior) as a simplified projective projection
  * @param pose Extrinsic matrix
  * @return Concatenation of intrinsic matrix and extrinsic matrix
  */
  virtual Mat34 get_projective_equivalent( const geometry::Pose3 & pose ) const = 0;
  
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
  * @brief Generate a unique Hash from the camera parameters (used for grouping)
  * @return Hash value
  */
  virtual std::size_t hashValue() const
  {
    size_t seed = 0;
    stl::hash_combine( seed, static_cast<int>( this->getType() ) );
    stl::hash_combine( seed, w_ );
    stl::hash_combine( seed, h_ );
    const std::vector<double> params = this->getParams();
    for ( const auto & param : params )
      stl::hash_combine( seed , param );
    return seed;
  }
};


/**
* @brief Compute angle between two bearing rays
* Bearing rays are computed from position on image plane in each cameras
*
* @param pose1 Pose of the first camera
* @param intrinsic1 Intrinsic of the first camera
* @param pose2 Pose of the second camera
* @param intrinsic2 Intrinsic of the second camera
* @param x1 Image coordinate of a point in first camera
* @param x2 Image coordinate of a point in the second camera
*
* @return Angle (in degree) between the two rays
*/
inline double AngleBetweenRay(
  const geometry::Pose3 & pose1,
  const IntrinsicBase * intrinsic1,
  const geometry::Pose3 & pose2,
  const IntrinsicBase * intrinsic2,
  const Vec2 & x1, const Vec2 & x2 )
{
  // x = (u, v, 1.0)  // image coordinates
  // X = R.t() * K.inv() * x + C // Camera world point
  // getting the ray:
  // ray = X - C = R.t() * K.inv() * x
  const Vec3 ray1 = ( pose1.rotation().transpose() * intrinsic1->operator()( x1 ) ).normalized();
  const Vec3 ray2 = ( pose2.rotation().transpose() * intrinsic2->operator()( x2 ) ).normalized();
  const double mag = ray1.norm() * ray2.norm();
  const double dotAngle = ray1.dot( ray2 );
  return R2D( acos( clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 ) ) );
}

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_INTRINSICS_HPP
