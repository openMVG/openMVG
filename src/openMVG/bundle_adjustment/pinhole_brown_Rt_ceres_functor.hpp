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

namespace pinhole_brown_reprojectionError {

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

/**
 * @brief Ceres functor to refine a pinhole camera using a brown radial
 *  distortion model, 3D points and focal.
 * 
 *  - first the intrinsic data block [focal;ppx;ppy;k1;k2;k3]
 *  - second the camera extrinsic block(camera orientation and position) [R;t]
 *    - 3 for rotation (angle axis), 3 for translation.
 *  - third the 3D point data block
 * 
 * @see ApplyRadialDistortionIntrinsics
 */

struct ErrorFunc_Refine_Camera_3DPoints
{
  ErrorFunc_Refine_Camera_3DPoints(const double* const pos_2dpoint)
  {
    m_pos_2dpoint[0] = pos_2dpoint[0];
    m_pos_2dpoint[1] = pos_2dpoint[1];
  }

  /// Compute the residual error after reprojection
  /// residual = observed - BrownDistortion([R|t] X)
  /**
   * @param[in] cam_intrinsics: Camera intrinsics
   * @param[in] cam_Rt: Camera extrinsics using one block of 6 parameters [R;t]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_Rt,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    // Rt[0,1,2] : angle-axis camera rotation.
    T x[3];
    ceres::AngleAxisRotatePoint(cam_Rt, pos_3dpoint, x);

    // Rt[3,4,5] : the camera translation.
    x[0] += cam_Rt[3];
    x[1] += cam_Rt[4];
    x[2] += cam_Rt[5];

    // Point from homogeneous to euclidean.
    T xe = x[0] / x[2];
    T ye = x[1] / x[2];

    // Unpack the intrinsics.
    const T& focal_length      = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1                = cam_intrinsics[OFFSET_K1];
    const T& k2                = cam_intrinsics[OFFSET_K2];
    const T& k3                = cam_intrinsics[OFFSET_K3];
//  const T& p1                = cam_intrinsics[OFFSET_P1];
//  const T& p2                = cam_intrinsics[OFFSET_P2];

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
    out_residuals[0] = predicted_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = predicted_y - T(m_pos_2dpoint[1]);

    return true;
  }

  double m_pos_2dpoint[2]; // The 2D observation
};

} // namespace pinhole_brown_reprojectionError
} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_BUNDLE_ADJUSTMENT_PINHOLE_BROWN_RT_CERES_FUNCTOR_HPP
