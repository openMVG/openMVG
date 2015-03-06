
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_INTRINSICS_H
#define OPENMVG_CAMERA_INTRINSICS_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/stl/hash.hpp"
#include <cereal/cereal.hpp>

#include <vector>

namespace openMVG{

enum EINTRINSIC
{
  PINHOLE_UNDEFINED = 0,      // camera with no parameters (only the W & H is known)
  PINHOLE_CAMERA    = 1,      // No distortion
  PINHOLE_CAMERA_RADIAL1 = 2, // radial distortion K1
  PINHOLE_CAMERA_RADIAL3 = 3, // radial distortion K1,K2,K3
  PINHOLE_CAMERA_BROWN   = 4  // radial + tangential
};

/// Basis class for all intrinsic parameters of a camera
/// Store the image size & define all basis optical modelization of a camera
struct IntrinsicBase
{
  unsigned int _w, _h;

  IntrinsicBase(unsigned int w = 0, unsigned int h = 0):_w(w), _h(h) {}

  const unsigned int w() const {return _w;}
  const unsigned int h() const {return _h;}

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

  /// Apply the distortion field to a point (that is in normalized camera frame)
  virtual Vec2 apply(const Vec2& p) const = 0;

  /// Remove the distortion to a camera point (that is in normalized camera frame)
  virtual Vec2 remove(const Vec2& p) const  = 0;

  /// Generate a unique Hash from the camera parameters (used for grouping)
  virtual std::size_t hashValue() const = 0;

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
};

/// Define a classic Pinhole camera (store a K 3x3 matrix) 
///  with intrinsic parameters defining the K calibration matrix
class Pinhole_Intrinsic : public IntrinsicBase
{
  protected:
    // Focal & principal point are embed into the calibration matrix K
    Mat3 _K, _Kinv;

  public:
  Pinhole_Intrinsic(
    unsigned int w = 0, unsigned int h = 0,
    double focal_length_pix = 0.0,
    double ppx = 0.0, double ppy = 0.0)
    :IntrinsicBase(w,h)
  {
    _K << focal_length_pix, 0., ppx, 0., focal_length_pix, ppy, 0., 0., 1.;
    _Kinv = _K.inverse();
  }

  virtual EINTRINSIC getType() const { return PINHOLE_CAMERA; }

  const Mat3& K() const { return _K; }
  const Mat3& Kinv() const { return _Kinv; }
  /// Return the value of the focal in pixels
  const double focal() const {return _K(0,0);}
  const Vec2 principal_point() const {return Vec2(_K(0,2), _K(1,2));}

  // Get bearing vector of p point (image coord)
  Vec3 operator () (const Vec2& p) const
  {
    Vec3 p3(p(0),p(1),1.0);
    return (_Kinv * p3).normalized();
  }

  // Transform a point from the camera plane to the image plane
  Vec2 cam2ima(const Vec2& p) const
  {
    // (focal * p) + principal point
    return _K(0,0) * p + Vec2(_K(0,2), _K(1,2));
  }

  // Transform a point from the image plane to the camera plane
  Vec2 ima2cam(const Vec2& p) const
  {
    // (p - principal_point) / focal
    return ( p -  Vec2(_K(0,2), _K(1,2)) ) / _K(0,0);
  }

  virtual bool have_disto() const {  return false; }

  virtual Vec2 apply(const Vec2& p) const  { return p; }

  virtual Vec2 remove(const Vec2& p) const  { return p; }

  virtual std::size_t hashValue() const
  {
    size_t seed;
    std::hash_combine(seed, static_cast<int>(this->getType()));
    std::hash_combine(seed, _w);
    std::hash_combine(seed, _h);
    const std::vector<double> params = this->getParams();
    for (size_t i=0; i < params.size(); ++i)
      std::hash_combine(seed, params[i]);
    return seed;
  }

  // Data wrapper for non linear optimization (get data)
  virtual std::vector<double> getParams() const
  {
    const std::vector<double> params = {_K(0,0), _K(0,2), _K(1,2)};
    return params;
  }

  // Data wrapper for non linear optimization (update from data)
  virtual bool updateFromParams(const std::vector<double> & params)
  {
    if (params.size() == 3) {
      *this = Pinhole_Intrinsic(_w, _h, params[0], params[1], params[2]);
      return true;
    }
    else  {
      return false;
    }
  }

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
    IntrinsicBase::save(ar);
    ar(cereal::make_nvp("focal_length", _K(0,0) ));
    const std::vector<double> pp = {_K(0,2), _K(1,2)};
    ar(cereal::make_nvp("principal_point", pp));
  }

  // Serialization
  template <class Archive>
  void load( Archive & ar)
  {
    IntrinsicBase::load(ar);
    double focal_length;
    ar(cereal::make_nvp("focal_length", focal_length ));
    std::vector<double> pp(2);
    ar(cereal::make_nvp("principal_point", pp));
    *this = Pinhole_Intrinsic(_w, _h, focal_length, pp[0], pp[1]);
  }
};

} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::Pinhole_Intrinsic, "pinhole");

#endif // #ifndef OPENMVG_CAMERA_INTRINSICS_H

