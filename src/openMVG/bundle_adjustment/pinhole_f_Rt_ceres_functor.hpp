// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_F_RT_CERES_FUNCTOR_HPP
#define OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_F_RT_CERES_FUNCTOR_HPP

#include "ceres/ceres.h"
#include "ceres/rotation.h"

// Definition of a Functor for minimization of the reprojection error
// (||x^j_i - P_{ij}(X_j)||)
// A 3D point X_j is projected on a image plane i and compare to the observation
//  x^j_i.

namespace openMVG{
namespace bundle_adjustment{

/// Templated functor for pinhole camera model to be used with Ceres.
/// The camera is parameterized using two blocks :
///   - first a block of one value for intrinsic [focal]
///   - second for camera position and orientation [R;t]
///     - 3 for rotation(angle axis), 3 for translation.
///   Principal point is assumed being applied on observed points.
///
struct Pinhole_f_Rt_ReprojectionError {
  Pinhole_f_Rt_ReprojectionError(double observed_x, double observed_y)
      : m_observed_x(observed_x), m_observed_y(observed_y) {}

  /// Compute the residual error after reprojection
  /// residual = observed - euclidean( f * [R|t] X)
  template <typename T>
  bool operator()(
    const T* const f,
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

    // Compute final projected point position.
    T predicted_x = (*f) * xe;
    T predicted_y = (*f) * ye;

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - T(m_observed_x);
    residuals[1] = predicted_y - T(m_observed_y);

    return true;
  }

  double m_observed_x, m_observed_y;
};

} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_F_RT_CERES_FUNCTOR_HPP
