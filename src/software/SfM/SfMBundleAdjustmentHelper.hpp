
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_BUNDLE_ADJUSTMENT_HELPER_H
#define OPENMVG_SFM_BUNDLE_ADJUSTMENT_HELPER_H

#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace openMVG{
namespace bundleAdjustment{

/// Container for a Bundle Adjustment dataset
///  Allows to refine NCamParam per camera and the structure (3Dpts)
///  Used with 7 parameter by default (Rotation(angle,axis), t, focal)
template<unsigned char NCamParam = 7>
class BAProblem {
 public:

  size_t num_observations()    const { return num_observations_; }
  const double* observations() const { return &observations_[0]; }
  double* mutable_cameras()          { return &parameters_[0];   }
  double* mutable_points()           { return &parameters_[0] + NCamParam * num_cameras_; }

  double* mutable_camera_for_observation(size_t i) {
    return mutable_cameras() + camera_index_[i] * NCamParam;
  }
  double* mutable_point_for_observation(size_t i) {
    return mutable_points() + point_index_[i] * 3;
  }
  size_t num_cameras_;      // # of cameras
  size_t num_points_;       // # of 3D points
  size_t num_observations_; // # of observations
  size_t num_parameters_;   // # of parameters ( NCamParam * #Cam + 3 * #Points)

  std::vector<size_t> point_index_;  // re-projection linked to the Inth 2d point
  std::vector<size_t> camera_index_; // re-projection linked to the Inth camera
  std::vector<double> observations_; // 3D points
  // Camera parametrization
  std::vector<double> parameters_;
};

/// Templated functor for pinhole camera model to be used with Ceres.
/// The camera is parameterized using 7 parameters:
///   - 3 for rotation(angle axis), 3 for translation, 1 for the focal
///    length. Principal point is assumed be located at the image center
///
struct PinholeReprojectionError {
  PinholeReprojectionError(double observed_x, double observed_y)
      : observed_x(observed_x), observed_y(observed_y) {}

  /// Compute the residual error after reprojection
  /// residual = observed - euclidean( focal * [R|t] point)
  template <typename T>
  bool operator()(const T* const camera,
                  const T* const point,
                  T* residuals) const {
    // camera[0,1,2] are the angle-axis rotation.
    T p[3];
    ceres::AngleAxisRotatePoint(camera, point, p);

    // camera[3,4,5] are the translation.
    p[0] += camera[3];
    p[1] += camera[4];
    p[2] += camera[5];

    // Point from homogeneous to euclidean
    T xp = p[0] / p[2];
    T yp = p[1] / p[2];

    // Compute final projected point position.
    const T& focal = camera[6];
    T predicted_x = focal * xp;
    T predicted_y = focal * yp;

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - T(observed_x);
    residuals[1] = predicted_y - T(observed_y);

    return true;
  }

  double observed_x;
  double observed_y;
};

} // namespace bundleAdjustment
} // namespace openMVG

#endif // OPENMVG_SFM_BUNDLE_ADJUSTMENT_HELPER_H
