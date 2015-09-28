
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_INTRINSICS_H
#define OPENMVG_CAMERA_INTRINSICS_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/stl/hash.hpp"
#include <vector>
#include <cereal/cereal.hpp> // Serialization

namespace openMVG {
namespace cameras {

/// Basis class for all intrinsic parameters of a camera
/// Store the image size & define all basis optical modelization of a camera
struct IntrinsicBase
{
  unsigned int _w, _h;

  IntrinsicBase(unsigned int w = 0, unsigned int h = 0):_w(w), _h(h) {}
  virtual ~IntrinsicBase() {}

  unsigned int w() const {return _w;}
  unsigned int h() const {return _h;}

  /// Projection of a 3D point into the camera plane (Apply pose, disto (if any) and Intrinsics)
  Vec2 project(
    const geometry::Pose3 & pose,
    const Vec3 & pt3D) const
  {
    const Vec3 X = pose(pt3D); // apply pose
    if (this->have_disto()) // apply disto & intrinsics
      return this->cam2ima( this->add_disto(X.head<2>()/X(2)) );
    else // apply intrinsics
      return this->cam2ima( X.head<2>()/X(2) );
  }

  /// Compute the residual between the 3D projected point X and an image observation x
  Vec2 residual(
    const geometry::Pose3 & pose,
    const Vec3 & X,
    const Vec2 & x) const
  {
    const Vec2 proj = this->project(pose, X);
    return x - proj;
  }

  // --
  // Virtual members
  // --

  /// Tell from which type the embed camera is
  virtual EINTRINSIC getType() const = 0;

  /// Data wrapper for non linear optimization (get data)
  virtual std::vector<double> getParams() const = 0;

  /// Data wrapper for non linear optimization (update from data)
  virtual bool updateFromParams(const std::vector<double> & params) = 0;

  /// Get bearing vector of p point (image coord)
  virtual Vec3 operator () (const Vec2& p) const = 0;

  /// Transform a point from the camera plane to the image plane
  virtual Vec2 cam2ima(const Vec2& p) const = 0;

  /// Transform a point from the image plane to the camera plane
  virtual Vec2 ima2cam(const Vec2& p) const = 0;

  /// Does the camera model handle a distortion field?
  virtual bool have_disto() const {return false;}

  /// Add the distortion field to a point (that is in normalized camera frame)
  virtual Vec2 add_disto(const Vec2& p) const = 0;

  /// Remove the distortion to a camera point (that is in normalized camera frame)
  virtual Vec2 remove_disto(const Vec2& p) const  = 0;

  /// Return the un-distorted pixel (with removed distortion)
  virtual Vec2 get_ud_pixel(const Vec2& p) const = 0;

  /// Return the distorted pixel (with added distortion)
  virtual Vec2 get_d_pixel(const Vec2& p) const = 0;

  /// Normalize a given unit pixel error to the camera plane
  virtual double imagePlane_toCameraPlaneError(double value) const = 0;

  /// Return the intrinsic (interior & exterior) as a simplified projective projection
  virtual Mat34 get_projective_equivalent(const geometry::Pose3 & pose) const = 0;

  /// Serialization out
  template <class Archive>
  void save( Archive & ar) const
  {
    ar(cereal::make_nvp("width", _w));
    ar(cereal::make_nvp("height", _h));
  }

  /// Serialization in
  template <class Archive>
  void load( Archive & ar)
  {
    ar(cereal::make_nvp("width", _w));
    ar(cereal::make_nvp("height", _h));
  }

  /// Generate an unique Hash from the camera parameters (used for grouping)
  virtual std::size_t hashValue() const
  {
    size_t seed = 0;
    stl::hash_combine(seed, static_cast<int>(this->getType()));
    stl::hash_combine(seed, _w);
    stl::hash_combine(seed, _h);
    const std::vector<double> params = this->getParams();
    for (size_t i=0; i < params.size(); ++i)
      stl::hash_combine(seed, params[i]);
    return seed;
  }
};

/// Return the angle (degree) between two bearing vector rays
static double AngleBetweenRay(
  const geometry::Pose3 & pose1,
  const IntrinsicBase * intrinsic1,
  const geometry::Pose3 & pose2,
  const IntrinsicBase * intrinsic2,
  const Vec2 & x1, const Vec2 & x2)
{
  // x = (u, v, 1.0)  // image coordinates
  // X = R.t() * K.inv() * x + C // Camera world point
  // getting the ray:
  // ray = X - C = R.t() * K.inv() * x
  const Vec3 ray1 = (pose1.rotation().transpose() * intrinsic1->operator()(x1)).normalized();
  const Vec3 ray2 = (pose2.rotation().transpose() * intrinsic2->operator()(x2)).normalized();
  const double mag = ray1.norm() * ray2.norm();
  const double dotAngle = ray1.dot(ray2);
  return R2D(acos(clamp(dotAngle/mag, -1.0 + 1.e-8, 1.0 - 1.e-8)));
}

} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERA_INTRINSICS_H

