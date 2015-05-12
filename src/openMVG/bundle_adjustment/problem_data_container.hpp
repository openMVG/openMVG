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
///  Used with 7 parameters by default:
///   (Rotation(angle,axis), t, focal)
template<unsigned char NCamParam = 7>
class BA_Problem_data {
 public:

  // Number of camera parameters
  static const unsigned char NCAMPARAM = NCamParam;

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

/// Container for a Bundle Adjustment dataset
/// Allows to refine cameras, shared intrinsics and the structure
/// Camera can be parametrized by the number of desired values
///   External parameter => 6: [Rotation(angle,axis), t]
///   Intrinsic => 1: [Focal]
template<
  unsigned char NExternalParam = 6,
  unsigned char NIntrinsicParam = 1>
class BA_Problem_data_camMotionAndIntrinsic {
 public:

  // Number of camera parameters
  static const unsigned char NEXTERNALPARAM = NExternalParam;
  static const unsigned char NINTRINSICPARAM = NIntrinsicParam;

  /// Return the number of observed 3D points
  size_t num_observations()    const { return num_observations_; }
  /// Return a pointer to observed points [X_0, ... ,X_n]
  const double* observations() const { return &observations_[0]; }

  /// Return the number of extrinsic groups
  size_t num_extrinsics()  const { return num_cameras_ ;}
  /// Return the number of intrinsic groups
  size_t num_intrinsics()  const { return num_intrinsics_;}

  /// Return a pointer to external camera data
  double* mutable_cameras_extrinsic() {
    return &parameters_[0];}
  /// Return a pointer to intrinsic camera data
  double* mutable_cameras_intrinsic() {
    return &parameters_[0] + NExternalParam * num_cameras_;}
  /// Return a pointer to 3D points data
  double* mutable_points()  {
    return &parameters_[0]
      + NExternalParam * num_cameras_
      + NIntrinsicParam * num_intrinsics_;}

  /// Return a pointer to the camera extrinsic that observe the Inth observation
  double* mutable_camera_extrinsic_for_observation(size_t i) {
    return mutable_cameras_extrinsic() + camera_index_extrinsic[i] * NExternalParam;
  }
  /// Return a pointer to the camera intrinsic that observe the Inth observation
  double* mutable_camera_intrinsic_for_observation(size_t i) {
    return mutable_cameras_intrinsic() + camera_index_intrinsic[i] * NIntrinsicParam;
  }

  /// Return a pointer to the Inth intrinsic parameters
  double* mutable_cameras_intrinsic(size_t i) {
    return mutable_cameras_intrinsic() + i * NIntrinsicParam;
  }
  /// Return a pointer to the Inth extrinsic parameters
  double* mutable_cameras_extrinsic(size_t i) {
    return mutable_cameras_extrinsic() + i * NExternalParam;
  }


  /// Return a pointer to the point that observe the Inth observation
  double* mutable_point_for_observation(size_t i) {
    return mutable_points() + point_index_[i] * 3;
  }

  size_t num_cameras_;      // # of cameras
  size_t num_intrinsics_;    // # of intrinsic groups
  size_t num_points_;       // # of 3D points
  size_t num_observations_; // # of observations
  // # of parameters: NIntrinsicParam * num_cameras_ + NIntrinsicParam * num_intrinsics_ + 3 * num_points_
  size_t num_parameters_;

  std::vector<size_t> point_index_;  // re-projection linked to the Inth 2d point
  std::vector<size_t> camera_index_extrinsic; // re-projection linked to the Inth camera extrinsic
  std::vector<size_t> camera_index_intrinsic; // re-projection linked to the Inth camera intrinsic
  std::vector<double> observations_; // 3D points

  // Camera parametrization ([R|t]_0,...,[R|t]_n,[f]_0,...,[f]_n)
  std::vector<double> parameters_;

};

} // namespace bundle_adjustment
} // namespace openMVG

#endif // OPENMVG_BUNDLE_ADJUSTMENT_PROBLEM_DATA_CONTAINER_HPP
