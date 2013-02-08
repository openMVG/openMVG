
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SIMPLE_CAMERA_H
#define OPENMVG_SFM_SIMPLE_CAMERA_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"

namespace openMVG{

// Modelize a pinhole camera P = K[R|t], t = -RC
struct SimpleCamera
{
  SimpleCamera(
    const Mat3 & K = Mat3::Identity(),
    const Mat3 & R = Mat3::Identity(),
    const Vec3 & t = Vec3::Zero())
    : _K(K), _R(R), _t(t)
  {
    _C = -R.transpose() * t;
    P_From_KRt(_K, _R, _t, &_P);
  }

  /// Projection matrix P = K[R|t]
  Mat34 _P;

  /// Intrinsic parameter (Focal, principal point)
  Mat3 _K;

  /// Extrinsic Rotation
  Mat3 _R;

  /// Extrinsic translation
  Vec3 _t;

  /// Camera center
  Vec3 _C;
};

} // namespace openMVG

#endif // OPENMVG_SFM_SIMPLE_CAMERA_H

