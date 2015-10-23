// This is an adaptation of the Fisheye distortion model implemented in OpenCV
// https://github.com/Itseez/opencv/blob/master/modules/calib3d/src/fisheye.cpp

// Copyright (c) 2015 Romain Janvier <romain.janvier~AT~univ-orleans.fr> for the given adaptation

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_PINHOLE_FISHEYE_HPP
#define OPENMVG_CAMERA_PINHOLE_FISHEYE_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Common.hpp"

#include <vector>

namespace openMVG {
namespace cameras {

/// Implement a simple Fish-eye camera model
class Pinhole_Intrinsic_Fisheye : public Pinhole_Intrinsic
{
  protected:
  // center of distortion is applied by the Intrinsics class
  std::vector<double> _params; // K1, K2, K3, K4


  public:

  Pinhole_Intrinsic_Fisheye(
    int w = 0, int h = 0,
    double focal = 0.0, double ppx = 0, double ppy = 0,
    double k1 = 0.0, double k2 = 0.0, double k3 = 0.0, double k4 = 0.0)
        :Pinhole_Intrinsic(w, h, focal, ppx, ppy)
  {
    _params = {k1, k2, k3, k4};
  }

  EINTRINSIC getType() const { return PINHOLE_CAMERA_FISHEYE; }

  virtual bool have_disto() const { return true;}

  virtual Vec2 add_disto(const Vec2 & p) const{
    const double eps = 1e-8;
    const double k1 = _params[0], k2 = _params[1], k3 = _params[2], k4 = _params[3];
    const double r = std::sqrt(p(0)*p(0) + p(1)*p(1));
    const double theta = std::atan(r);
    const double
      theta2 = theta*theta,
      theta3 = theta2*theta,
      theta4 = theta2*theta2,
      theta5 = theta4*theta,
      theta6 = theta3*theta3,
      theta7 = theta6*theta,
      theta8 = theta4*theta4,
      theta9 = theta8*theta;
    const double theta_dist = theta + k1*theta3 + k2*theta5 + k3*theta7 + k4*theta9;
    const double inv_r = r > eps ? 1.0/r : 1.0;
    const double cdist = r > eps ? theta_dist * inv_r : 1.0;
    return  p*cdist;
  }

  virtual Vec2 remove_disto(const Vec2 & p) const{
    const double eps = 1e-8;
    double scale = 1.0;
    const double theta_dist = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    if (theta_dist > eps)
    {
      double theta = theta_dist;
      for (int j = 0; j < 10; ++j)
      {
        const double
          theta2 = theta*theta,
          theta4 = theta2*theta2,
          theta6 = theta4*theta2,
          theta8 = theta6*theta2;
        theta = theta_dist /
          (1 + _params[0] * theta2
             + _params[1] * theta4
             + _params[2] * theta6
             + _params[3] * theta8);
      }
      scale = std::tan(theta) / theta_dist;
    }
    return p * scale;
  }

  // Data wrapper for non linear optimization (get data)
  virtual std::vector<double> getParams() const
  {
    std::vector<double> params = Pinhole_Intrinsic::getParams();
    params.push_back(_params[0]);
    params.push_back(_params[1]);
    params.push_back(_params[2]);
    params.push_back(_params[3]);
    return params;
  }

  // Data wrapper for non linear optimization (update from data)
  virtual bool updateFromParams(const std::vector<double> & params)
  {
    if (params.size() == 7) {
      *this = Pinhole_Intrinsic_Fisheye(
      _w, _h,
      params[0], params[1], params[2], // focal, ppx, ppy
      params[3], params[4], params[5], params[6]); // k1, k2, k3, k4
      return true;
    }
    else  {
      return false;
    }
  }

  /// Return the un-distorted pixel (with removed distortion)
  virtual Vec2 get_ud_pixel(const Vec2& p) const
  {
    return cam2ima( remove_disto(ima2cam(p)) );
  }

  /// Return the distorted pixel (with added distortion)
  virtual Vec2 get_d_pixel(const Vec2& p) const
  {
    return cam2ima( add_disto(ima2cam(p)) );
  }

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
    Pinhole_Intrinsic::save(ar);
    ar(cereal::make_nvp("fisheye", _params));
  }

  // Serialization
  template <class Archive>
  void load( Archive & ar)
  {
    Pinhole_Intrinsic::load(ar);
    ar(cereal::make_nvp("fisheye", _params));
  }
};


} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_Fisheye, "fisheye");

#endif // #ifndef OPENMVG_CAMERA_PINHOLE_FISHEYE_HPP
