
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_PINHOLE_HPP
#define OPENMVG_CAMERA_PINHOLE_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/geometry/pose3.hpp"

#include <vector>

namespace openMVG {
namespace cameras {

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

  virtual ~Pinhole_Intrinsic() {}

  virtual EINTRINSIC getType() const { return PINHOLE_CAMERA; }

  const Mat3& K() const { return _K; }
  const Mat3& Kinv() const { return _Kinv; }
  /// Return the value of the focal in pixels
  inline double focal() const {return _K(0,0);}
  inline Vec2 principal_point() const {return Vec2(_K(0,2), _K(1,2));}

  // Get bearing vector of p point (image coord)
  Vec3 operator () (const Vec2& p) const
  {
    Vec3 p3(p(0),p(1),1.0);
    return (_Kinv * p3).normalized();
  }

  // Transform a point from the camera plane to the image plane
  Vec2 cam2ima(const Vec2& p) const
  {
    return focal() * p + principal_point();
  }

  // Transform a point from the image plane to the camera plane
  Vec2 ima2cam(const Vec2& p) const
  {
    return ( p -  principal_point() ) / focal();
  }

  virtual bool have_disto() const {  return false; }

  virtual Vec2 add_disto(const Vec2& p) const  { return p; }

  virtual Vec2 remove_disto(const Vec2& p) const  { return p; }

  virtual double imagePlane_toCameraPlaneError(double value) const
  {
    return value / focal();
  }

  virtual Mat34 get_projective_equivalent(const geometry::Pose3 & pose) const
  {
    Mat34 P;
    P_From_KRt(K(), pose.rotation(), pose.translation(), &P);
    return P;
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

  /// Return the un-distorted pixel (with removed distortion)
  virtual Vec2 get_ud_pixel(const Vec2& p) const {return p;}

  /// Return the distorted pixel (with added distortion)
  virtual Vec2 get_d_pixel(const Vec2& p) const {return p;}

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

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic, "pinhole");

#endif // #ifndef OPENMVG_CAMERA_PINHOLE_HPP

