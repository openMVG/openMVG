// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_RTF_CERES_FUNCTOR_HPP
#define OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_RTF_CERES_FUNCTOR_HPP

#include "ceres/ceres.h"
#include "ceres/rotation.h"

// Definition of a Functor for minimization of the reprojection error
// (||x^j_i - P_{ij}(X_j)||)
// A 3D point X_j is projected on a image plane i and compare to the observation
//  x^j_i.

namespace openMVG{
namespace bundle_adjustment{

/// Templated functor to refine a pinhole camera model and 3D points
///   to be used with Ceres.
/// The camera is parameterized using one block of 7 parameters [R;t;f]:
///   - 3 for rotation(angle axis), 3 for translation, 1 for the focal length.
///   Principal point is assumed being applied on observed points.
///
struct Pinhole_Rtf_ReprojectionError {
  Pinhole_Rtf_ReprojectionError(double observed_x, double observed_y)
      : m_observed_x(observed_x), m_observed_y(observed_y) {}

  /// Compute the residual error after reprojection
  /// residual = observed - euclidean( f * [R|t] X)
  template <typename T>
  bool operator()(
    const T* const Rtf, // [R;t;f]
    const T* const X,
    T* residuals) const
  {
    // Rtf[0,1,2] : angle-axis camera rotation
    T x[3];
    ceres::AngleAxisRotatePoint(Rtf, X, x);

    // Rtf[3,4,5] : the camera translation
    x[0] += Rtf[3];
    x[1] += Rtf[4];
    x[2] += Rtf[5];

    // Point from homogeneous to euclidean
    T xe = x[0] / x[2];
    T ye = x[1] / x[2];

    // Rtf[6] : the focal length
    const T& focal = Rtf[6];
    T predicted_x = focal * xe;
    T predicted_y = focal * ye;

    // The error is the difference between the predicted and observed position
    residuals[0] = predicted_x - T(m_observed_x);
    residuals[1] = predicted_y - T(m_observed_y);

    return true;
  }

  double m_observed_x, m_observed_y;
};

} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_RTF_CERES_FUNCTOR_HPP
