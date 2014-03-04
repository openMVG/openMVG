
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_BUNDLE_ADJUSTMENT_HELPER_TRANSLATION_ONLY_H
#define OPENMVG_SFM_BUNDLE_ADJUSTMENT_HELPER_TRANSLATION_ONLY_H

#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"
#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace openMVG{
namespace bundle_adjustment{

struct PinholeReprojectionError_t {
  PinholeReprojectionError_t(const double* const pos_2dpoint, double focal, const double * r)
    : m_focal(focal)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];

    m_rotation[0] = r[0];
    m_rotation[1] = r[1];
    m_rotation[2] = r[2];
  }

  template <typename T>
  bool operator()(const T* const cam_t,
                  const T* const pos_3dpoint,
                  T* out_residuals) const {

    T rot[3] = {T(m_rotation[0]), T(m_rotation[1]), T(m_rotation[2])};
    T focal (m_focal);
    pinhole_reprojectionError::computeResidual(
      &rot[0], // => cam_R
      &cam_t[0], // => cam_t
      &focal,
      pos_3dpoint,
      m_pos_2dpoint,
      out_residuals );

    return true;
  }

  // Static data
  double m_pos_2dpoint[2];  // observation
  double m_focal;
  double m_rotation[3];
};

struct PinholeReprojectionError_Rt {
  PinholeReprojectionError_Rt(const double* const pos_2dpoint, double focal)
    : m_focal(focal)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  template <typename T>
  bool operator()(
    const T* const cam_Rt, // [R;t]
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    T focal (m_focal);
    pinhole_reprojectionError::computeResidual(
      &cam_Rt[0], // => cam_R
      &cam_Rt[3], // => cam_t
      &focal,
      pos_3dpoint,
      m_pos_2dpoint,
      out_residuals );

    return true;
  }

  // Static data
  double m_pos_2dpoint[2]; // observation
  double m_focal;
};

} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_SFM_BUNDLE_ADJUSTMENT_HELPER_TRANSLATION_ONLY_H
