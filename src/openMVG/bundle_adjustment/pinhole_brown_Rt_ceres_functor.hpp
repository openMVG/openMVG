// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_BROWN_RT_CERES_FUNCTOR_HPP
#define OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_BROWN_RT_CERES_FUNCTOR_HPP

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "openMVG/cameras/BrownPinholeCamera.hpp"

// Definition of a Functor for minimization of the reprojection error
//  according the Brown's distortion model.
// A decentered distortion model.
// (||x^j_i - disto(f,ppx,ppx,k1,k2,k3,[R|t]_{ij}(X_j))||)
// A 3D point X_j is projected on a image plane i and compare to the observation
//  x^j_i in the original image.

namespace openMVG{
namespace bundle_adjustment{

// Enum to order the intrinsics parameters into a single parameter block
enum {
  OFFSET_FOCAL_LENGTH,
  OFFSET_PRINCIPAL_POINT_X,
  OFFSET_PRINCIPAL_POINT_Y,
  OFFSET_K1,
  OFFSET_K2,
  OFFSET_K3//,
//  OFFSET_P1,
//  OFFSET_P2,
};

/// Templated functor for pinhole camera model to be used with Ceres.
/// The camera is parameterized using two blocks :
///   - first a block of intrinsic values [focal;ppx;ppy;k1;k2;k3]
///   - second for camera position and orientation [R;t]
///     - 3 for rotation(angle axis), 3 for translation.
///   Principal point is assumed being applied on observed points.
///
struct Pinhole_brown_Rt_ReprojectionError {
  Pinhole_brown_Rt_ReprojectionError(double observed_x, double observed_y)
      : m_observed_x(observed_x), m_observed_y(observed_y) {}

  /// Compute the residual error after reprojection
  /// residual = observed - BrownDistortion([R|t] X)
  template <typename T>
  bool operator()(
    const T* const intrinsics,
    const T* const Rt,
    const T* const X,
    T* residuals) const
  {
    // Rt[0,1,2] : angle-axis camera rotation.
    T x[3];
    ceres::AngleAxisRotatePoint(Rt, X, x);

    // Rt[3,4,5] : the camera translation.
    x[0] += Rt[3];
    x[1] += Rt[4];
    x[2] += Rt[5];

    // Point from homogeneous to euclidean.
    T xe = x[0] / x[2];
    T ye = x[1] / x[2];

    // Unpack the intrinsics.
    const T& focal_length      = intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1                = intrinsics[OFFSET_K1];
    const T& k2                = intrinsics[OFFSET_K2];
    const T& k3                = intrinsics[OFFSET_K3];
//    const T& p1                = intrinsics[OFFSET_P1];
//    const T& p2                = intrinsics[OFFSET_P2];

    T predicted_x, predicted_y;

    // Apply distortion to the normalized points.
    ApplyRadialDistortionIntrinsics(
      focal_length,
      focal_length,
      principal_point_x,
      principal_point_y,
      k1, k2, k3,
      xe, ye,
      &predicted_x,
      &predicted_y);

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - T(m_observed_x);
    residuals[1] = predicted_y - T(m_observed_y);

    return true;
  }

  double m_observed_x, m_observed_y;
};

} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_BROWN_RT_CERES_FUNCTOR_HPP
