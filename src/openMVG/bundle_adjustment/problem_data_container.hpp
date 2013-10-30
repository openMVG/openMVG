// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_BUNDLE_ADJUSTMENT_PROBLEM_DATA_CONTAINER_HPP
#define OPENMVG_BUNDLE_ADJUSTMENT_PROBLEM_DATA_CONTAINER_HPP

#include <vector>

namespace openMVG{
namespace bundle_adjustment{

/// Container for a Bundle Adjustment dataset
///  Allows to refine NCamParam per camera and the structure (3Dpts)
///  Used with 7 parameter by default (Rotation(angle,axis), t, focal)
template<unsigned char NCamParam = 7>
class BA_Problem_data_container {
 public:
  /// Return the number of observed 3D points
  size_t num_observations()    const { return num_observations_; }
  /// Return a pointer to observed points [X_0, ... ,X_n]
  const double* observations() const { return &observations_[0]; }
  /// Return pointer to camera data
  double* mutable_cameras() {
    return &parameters_[0];}
  /// Return point to points data
  double* mutable_points()  {
    return &parameters_[0] + NCamParam * num_cameras_;}

  /// Return a pointer to the camera that observe the Inth observation
  double* mutable_camera_for_observation(size_t i) {
    return mutable_cameras() + camera_index_[i] * NCamParam;
  }
  /// Return a pointer to the point that observe the Inth observation
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

  std::vector<double> parameters_;  // Camera parametrization
};

} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_BUNDLE_ADJUSTMENT_PROBLEM_DATA_CONTAINER_HPP
