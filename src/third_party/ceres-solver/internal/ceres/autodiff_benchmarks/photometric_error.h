// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2020 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: nikolaus@nikolaus-demmel.de (Nikolaus Demmel)
//
//
#ifndef CERES_INTERNAL_AUTODIFF_BENCHMARK_PHOTOMETRIC_ERROR_H_
#define CERES_INTERNAL_AUTODIFF_BENCHMARK_PHOTOMETRIC_ERROR_H_

#include <Eigen/Dense>

#include "ceres/cubic_interpolation.h"

namespace ceres {

// Photometric residual that computes the intensity difference for a patch
// between host and target frame. The point is parameterized with inverse
// distance relative to the host frame. The relative pose between host and
// target frame is computed from their respective absolute poses.
//
// The residual is similar to the one defined by Engel et al. [1]. Differences
// include:
//
// 1. Use of a camera model based on spherical projection, namely the enhanced
// unified camera model [2][3]. This is intended to bring some variability to
// the benchmark compared to the SnavelyReprojection that uses a
// polynomial-based distortion model.
//
// 2. To match the camera model, inverse distance parameterization is used for
// points instead of inverse depth [4].
//
// 3. For simplicity, camera intrinsics are assumed constant, and thus host
// frame points are passed as (unprojected) bearing vectors, which avoids the
// need for an 'unproject' function.
//
// 4. Some details of the residual in [1] are omitted for simplicity: The
// brightness transform parameters [a,b], the constant pre-weight w, and the
// per-pixel robust norm.
//
// [1] J. Engel, V. Koltun and D. Cremers, "Direct Sparse Odometry," in IEEE
// Transactions on Pattern Analysis and Machine Intelligence, vol. 40, no. 3,
// pp. 611-625, 1 March 2018.
//
// [2] B. Khomutenko, G. Garcia and P. Martinet, "An Enhanced Unified Camera
// Model," in IEEE Robotics and Automation Letters, vol. 1, no. 1, pp. 137-144,
// Jan. 2016.
//
// [3] V. Usenko, N. Demmel and D. Cremers, "The Double Sphere Camera Model,"
// 2018 International Conference on 3D Vision (3DV), Verona, 2018, pp. 552-560.
//
// [4] H. Matsuki, L. von Stumberg, V. Usenko, J. Stückler and D. Cremers,
// "Omnidirectional DSO: Direct Sparse Odometry With Fisheye Cameras," in IEEE
// Robotics and Automation Letters, vol. 3, no. 4, pp. 3693-3700, Oct. 2018.
template <int PATCH_SIZE_ = 8>
struct PhotometricError {
  static constexpr int PATCH_SIZE = PATCH_SIZE_;
  static constexpr int POSE_SIZE = 7;
  static constexpr int POINT_SIZE = 1;

  using Grid = Grid2D<uint8_t, 1>;
  using Interpolator = BiCubicInterpolator<Grid>;
  using Intrinsics = Eigen::Array<double, 6, 1>;

  template <typename T>
  using Patch = Eigen::Array<T, PATCH_SIZE, 1>;

  template <typename T>
  using PatchVectors = Eigen::Matrix<T, 3, PATCH_SIZE>;

  PhotometricError(const Patch<double>& intensities_host,
                   const PatchVectors<double>& bearings_host,
                   const Interpolator& image_target,
                   const Intrinsics& intrinsics)
      : intensities_host_(intensities_host),
        bearings_host_(bearings_host),
        image_target_(image_target),
        intrinsics_(intrinsics) {}

  template <typename T>
  inline bool Project(Eigen::Matrix<T, 2, 1>& proj,
                      const Eigen::Matrix<T, 3, 1>& p) const {
    const double& fx = intrinsics_[0];
    const double& fy = intrinsics_[1];
    const double& cx = intrinsics_[2];
    const double& cy = intrinsics_[3];
    const double& alpha = intrinsics_[4];
    const double& beta = intrinsics_[5];

    const T rho2 = beta * (p.x() * p.x() + p.y() * p.y()) + p.z() * p.z();
    const T rho = sqrt(rho2);

    // Check if 3D point is in domain of projection function.
    // See (8) and (17) in [3].
    constexpr double NUMERIC_EPSILON = 1e-10;
    const double w =
        alpha > 0.5 ? (1.0 - alpha) / alpha : alpha / (1.0 - alpha);
    if (p.z() <= -w * rho + NUMERIC_EPSILON) {
      return false;
    }

    const T norm = alpha * rho + (1.0 - alpha) * p.z();
    const T norm_inv = 1.0 / norm;

    const T mx = p.x() * norm_inv;
    const T my = p.y() * norm_inv;

    proj[0] = fx * mx + cx;
    proj[1] = fy * my + cy;

    return true;
  }

  template <typename T>
  inline bool operator()(const T* const pose_host_ptr,
                         const T* const pose_target_ptr,
                         const T* const idist_ptr,
                         T* residuals_ptr) const {
    Eigen::Map<const Eigen::Quaternion<T>> q_w_h(pose_host_ptr);
    Eigen::Map<const Eigen::Matrix<T, 3, 1>> t_w_h(pose_host_ptr + 4);
    Eigen::Map<const Eigen::Quaternion<T>> q_w_t(pose_target_ptr);
    Eigen::Map<const Eigen::Matrix<T, 3, 1>> t_w_t(pose_target_ptr + 4);
    const T& idist = *idist_ptr;
    Eigen::Map<Patch<T>> residuals(residuals_ptr);

    // Compute relative pose from host to target frame.
    const Eigen::Quaternion<T> q_t_h = q_w_t.conjugate() * q_w_h;
    const Eigen::Matrix<T, 3, 3> R_t_h = q_t_h.toRotationMatrix();
    const Eigen::Matrix<T, 3, 1> t_t_h = q_w_t.conjugate() * (t_w_h - t_w_t);

    // Transform points from host to target frame. 3D point in target frame is
    // scaled by idist for numerical stability when idist is close to 0
    // (projection is invariant to scaling).
    PatchVectors<T> p_target_scaled =
        (R_t_h * bearings_host_).colwise() + idist * t_t_h;

    // Project points and interpolate image.
    Patch<T> intensities_target;
    for (int i = 0; i < p_target_scaled.cols(); ++i) {
      Eigen::Matrix<T, 2, 1> uv;
      if (!Project(uv, Eigen::Matrix<T, 3, 1>(p_target_scaled.col(i)))) {
        // If any point of the patch is outside the domain of the projection
        // function, the residual cannot be evaluated. For the benchmark we want
        // to avoid this case and thus throw an exception to indicate
        // immediately if it does actually happen after possible future changes.
        throw std::runtime_error("Benchmark data leads to invalid projection.");
      }

      // Mind the order of u and v: Evaluate takes (row, column), but u is
      // left-to-right and v top-to-bottom image axis.
      image_target_.Evaluate(uv[1], uv[0], &intensities_target[i]);
    }

    // Residual is intensity difference between host and target frame.
    residuals = intensities_target - intensities_host_;

    return true;
  }

 private:
  const Patch<double>& intensities_host_;
  const PatchVectors<double>& bearings_host_;
  const Interpolator& image_target_;
  const Intrinsics& intrinsics_;
};
}  // namespace ceres
#endif  // CERES_INTERNAL_AUTODIFF_BENCHMARK_PHOTOMETRIC_ERROR_H_
