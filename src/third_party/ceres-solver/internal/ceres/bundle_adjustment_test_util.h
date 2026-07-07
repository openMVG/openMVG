// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
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
// Author: keir@google.com (Keir Mierle)
//
// End-to-end bundle adjustment test utilities for Ceres. This base is used in
// the generated bundle adjustment test binaries. The reason to split the
// bundle tests into separate binaries is so the tests can get parallelized.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "ceres/autodiff_cost_function.h"
#include "ceres/internal/export.h"
#include "ceres/ordered_groups.h"
#include "ceres/problem.h"
#include "ceres/rotation.h"
#include "ceres/solver.h"
#include "ceres/stringprintf.h"
#include "ceres/test_util.h"
#include "ceres/types.h"
#include "glog/logging.h"

namespace ceres {
namespace internal {

const bool kAutomaticOrdering = true;
const bool kUserOrdering = false;

// This class implements the SystemTestProblem interface and provides
// access to a bundle adjustment problem. It is based on
// examples/bundle_adjustment_example.cc. Currently a small 16 camera
// problem is hard coded in the constructor.
class BundleAdjustmentProblem {
 public:
  BundleAdjustmentProblem(const std::string input_file) {
    ReadData(input_file);
    BuildProblem();
  }
  BundleAdjustmentProblem() {
    const std::string input_file =
        TestFileAbsolutePath("problem-16-22106-pre.txt");
    ReadData(input_file);
    BuildProblem();
  }

  ~BundleAdjustmentProblem() {
    delete[] point_index_;
    delete[] camera_index_;
    delete[] observations_;
    delete[] parameters_;
  }

  Problem* mutable_problem() { return &problem_; }
  Solver::Options* mutable_solver_options() { return &options_; }

  // clang-format off
  int num_cameras()                const { return num_cameras_; }
  int num_points()                 const { return num_points_; }
  int num_observations()           const { return num_observations_; }
  const int* point_index()         const { return point_index_; }
  const int* camera_index()        const { return camera_index_; }
  const double* observations()     const { return observations_; }
  double* mutable_cameras()              { return parameters_; }
  double* mutable_points()               { return parameters_ + 9 * num_cameras_; }
  const Solver::Options& options() const { return options_; }
  // clang-format on

  static double kResidualTolerance;

 private:
  void ReadData(const std::string& filename) {
    FILE* fptr = fopen(filename.c_str(), "r");

    if (!fptr) {
      LOG(FATAL) << "File Error: unable to open file " << filename;
    }

    // This will die horribly on invalid files. Them's the breaks.
    FscanfOrDie(fptr, "%d", &num_cameras_);
    FscanfOrDie(fptr, "%d", &num_points_);
    FscanfOrDie(fptr, "%d", &num_observations_);

    VLOG(1) << "Header: " << num_cameras_ << " " << num_points_ << " "
            << num_observations_;

    point_index_ = new int[num_observations_];
    camera_index_ = new int[num_observations_];
    observations_ = new double[2 * num_observations_];

    num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
    parameters_ = new double[num_parameters_];

    for (int i = 0; i < num_observations_; ++i) {
      FscanfOrDie(fptr, "%d", camera_index_ + i);
      FscanfOrDie(fptr, "%d", point_index_ + i);
      for (int j = 0; j < 2; ++j) {
        FscanfOrDie(fptr, "%lf", observations_ + 2 * i + j);
      }
    }

    for (int i = 0; i < num_parameters_; ++i) {
      FscanfOrDie(fptr, "%lf", parameters_ + i);
    }

    fclose(fptr);
  }

  void BuildProblem() {
    double* points = mutable_points();
    double* cameras = mutable_cameras();

    for (int i = 0; i < num_observations(); ++i) {
      // Each Residual block takes a point and a camera as input and
      // outputs a 2 dimensional residual.
      CostFunction* cost_function =
          new AutoDiffCostFunction<BundlerResidual, 2, 9, 3>(
              new BundlerResidual(observations_[2 * i + 0],
                                  observations_[2 * i + 1]));

      // Each observation corresponds to a pair of a camera and a point
      // which are identified by camera_index()[i] and
      // point_index()[i] respectively.
      double* camera = cameras + 9 * camera_index_[i];
      double* point = points + 3 * point_index()[i];
      problem_.AddResidualBlock(cost_function, nullptr, camera, point);
    }

    options_.linear_solver_ordering =
        std::make_shared<ParameterBlockOrdering>();

    // The points come before the cameras.
    for (int i = 0; i < num_points_; ++i) {
      options_.linear_solver_ordering->AddElementToGroup(points + 3 * i, 0);
    }

    for (int i = 0; i < num_cameras_; ++i) {
      options_.linear_solver_ordering->AddElementToGroup(cameras + 9 * i, 1);
    }

    options_.linear_solver_type = DENSE_SCHUR;
    options_.max_num_iterations = 25;
    options_.function_tolerance = 1e-10;
    options_.gradient_tolerance = 1e-10;
    options_.parameter_tolerance = 1e-10;
  }

  template <typename T>
  void FscanfOrDie(FILE* fptr, const char* format, T* value) {
    int num_scanned = fscanf(fptr, format, value);
    if (num_scanned != 1) {
      LOG(FATAL) << "Invalid UW data file.";
    }
  }

  // Templated pinhole camera model.  The camera is parameterized
  // using 9 parameters. 3 for rotation, 3 for translation, 1 for
  // focal length and 2 for radial distortion. The principal point is
  // not modeled (i.e. it is assumed to be located at the image
  // center).
  struct BundlerResidual {
    // (u, v): the position of the observation with respect to the image
    // center point.
    BundlerResidual(double u, double v) : u(u), v(v) {}

    template <typename T>
    bool operator()(const T* const camera,
                    const T* const point,
                    T* residuals) const {
      T p[3];
      AngleAxisRotatePoint(camera, point, p);

      // Add the translation vector
      p[0] += camera[3];
      p[1] += camera[4];
      p[2] += camera[5];

      const T& focal = camera[6];
      const T& l1 = camera[7];
      const T& l2 = camera[8];

      // Compute the center of distortion.  The sign change comes from
      // the camera model that Noah Snavely's Bundler assumes, whereby
      // the camera coordinate system has a negative z axis.
      T xp = -focal * p[0] / p[2];
      T yp = -focal * p[1] / p[2];

      // Apply second and fourth order radial distortion.
      T r2 = xp * xp + yp * yp;
      T distortion = T(1.0) + r2 * (l1 + l2 * r2);

      residuals[0] = distortion * xp - u;
      residuals[1] = distortion * yp - v;

      return true;
    }

    double u;
    double v;
  };

  Problem problem_;
  Solver::Options options_;

  int num_cameras_;
  int num_points_;
  int num_observations_;
  int num_parameters_;

  int* point_index_;
  int* camera_index_;
  double* observations_;
  // The parameter vector is laid out as follows
  // [camera_1, ..., camera_n, point_1, ..., point_m]
  double* parameters_;
};

double BundleAdjustmentProblem::kResidualTolerance = 1e-4;
using BundleAdjustmentTest = SystemTest<BundleAdjustmentProblem>;

}  // namespace internal
}  // namespace ceres
