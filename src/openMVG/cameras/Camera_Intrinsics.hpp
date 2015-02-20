
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERA_INTRINSICS_H
#define OPENMVG_CAMERA_INTRINSICS_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"
#include <vector>

namespace openMVG{

enum EINTRINSIC
{
  PINHOLE_CAMERA, // No distortion
  PINHOLE_CAMERA_RADIAL1, // radial distortion K1
  PINHOLE_CAMERA_RADIAL3, // radial distortion K1,K2,K3
  PINHOLE_CAMERA_BROWN // radial + tangential
};

class Intrinsic
{
  protected:
    int _w, _h;
    Mat3 _K, _Kinv;

  public:
  Intrinsic(int w = 0, int h = 0, double focal = 0.0, double ppx = 0, double ppy = 0):
  _w(w), _h(h)
  {
    _K << focal, 0, ppx, 0, focal, ppy, 0, 0, 1;
    _Kinv = _K.inverse();
  }

  virtual EINTRINSIC getType() const { return PINHOLE_CAMERA; }

  const Mat3& K() const { return _K; }
  const Mat3& Kinv() const { return _Kinv; }
  const int w() const {return _w;}
  const int h() const {return _h;}

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

  // Data wrapper for non linear optimization (get data)
  virtual std::vector<double> getParams() const
  {
    std::vector<double> params(3);
    params[0] = _K(0,0);
    params[1] = _K(0,2);
    params[2] = _K(1,2);
    return params;
  }

  // Data wrapper for non linear optimization (update from data)
  virtual bool updateFromParams(const std::vector<double> & params)
  {
    if (params.size() == 3) {
      *this = Intrinsic(_w, _h, params[0], params[1], params[2]);
      return true;
    }
    else  {
      return false;
    }
  }
};

} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERA_INTRINSICS_H

