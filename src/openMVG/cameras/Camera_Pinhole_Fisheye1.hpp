// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Common.hpp"

#include <vector>

namespace openMVG {
namespace cameras {

/**
 * Implement a simple Fish-eye camera model with only one parameter
 * 
 * Fredreric Devernay and Olivier Faugeras. 2001. Straight lines have to be 
 * straight: automatic calibration and removal of distortion from scenes of 
 * structured environments. Mach. Vision Appl. 13, 1 (August 2001), 14-24. 
 * DOI: 10.1007/PL00013269 https://hal.inria.fr/inria-00267247/document
 */
class Pinhole_Intrinsic_Fisheye1 : public Pinhole_Intrinsic
{
protected:
  // center of distortion is applied by the Intrinsics class
  std::vector<double> _distortionParams; // K1

public:

  Pinhole_Intrinsic_Fisheye1(
    int w = 0, int h = 0,
    double focal = 0.0, double ppx = 0, double ppy = 0,
    double k1 = 0.0)
        :Pinhole_Intrinsic(w, h, focal, ppx, ppy)
  {
    _distortionParams = {k1};
  }

  Pinhole_Intrinsic_Fisheye1* clone() const { return new Pinhole_Intrinsic_Fisheye1(*this); }
  void assign(const IntrinsicBase& other) { *this = dynamic_cast<const Pinhole_Intrinsic_Fisheye1&>(other); }

  EINTRINSIC getType() const { return PINHOLE_CAMERA_FISHEYE1; }

  virtual bool have_disto() const { return true;}

  virtual Vec2 add_disto(const Vec2 & p) const
  {
    const double k1 = _distortionParams[0];
    const double r = std::sqrt(p(0)*p(0) + p(1)*p(1));
    const double coef = (std::atan(2.0 * r * std::tan(0.5 * k1)) / k1) / r;
    return  p * coef;
  }

  virtual Vec2 remove_disto(const Vec2 & p) const
  {
    const double k1 = _distortionParams[0];
    const double r = std::sqrt(p(0)*p(0) + p(1)*p(1));
    const double coef = 0.5 * std::tan(r * k1) / (std::tan(0.5 * k1) * r);
    return  p * coef;
  }

  // Data wrapper for non linear optimization (get data)
  virtual std::vector<double> getParams() const
  {
    std::vector<double> params = Pinhole_Intrinsic::getParams();
    params.push_back(_distortionParams[0]);
    return params;
  }

  virtual std::vector<double> getDistortionParams() const
  {
    return _distortionParams;
  }

  // Data wrapper for non linear optimization (update from data)
  virtual bool updateFromParams(const std::vector<double> & params)
  {
    if (params.size() == 7)
    {
      this->setK(params[0], params[1], params[2]);
      _distortionParams = {params[3]};
      return true;
    }
    return false;
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
    ar(cereal::make_nvp("fisheye1", _distortionParams));
  }

  // Serialization
  template <class Archive>
  void load( Archive & ar)
  {
    Pinhole_Intrinsic::load(ar);
    ar(cereal::make_nvp("fisheye1", _distortionParams));
  }
};

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_Fisheye1, "fisheye1");
