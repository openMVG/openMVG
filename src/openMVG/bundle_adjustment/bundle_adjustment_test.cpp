
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


// An example of a minimal, self-contained bundle adjuster using Ceres
// It refines Focal, Rotation and Translation of the cameras.
// => A synthetic scene is used:
//    a random noise between [-.5,.5] is added on observed data points

#include "testing/testing.h"

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/multiview/projection.hpp"

// Bundle Adjustment includes
#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"
#include "openMVG/bundle_adjustment/problem_data_container.hpp"

using namespace openMVG;
using namespace openMVG::bundle_adjustment;

#include <cmath>
#include <cstdio>
#include <iostream>

TEST(BUNDLE_ADJUSTMENT, EffectiveMinimization_RTf) {

  int nviews = 3;
  int npoints = 6;
  NViewDataSet d = NRealisticCamerasRing(nviews, npoints);

  // Setup a BA problem
  BA_Problem_data<9> ba_problem;

  // Configure the size of the problem
  ba_problem.num_cameras_ = nviews;
  ba_problem.num_points_ = npoints;
  ba_problem.num_observations_ = nviews * npoints;

  ba_problem.point_index_.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_.reserve(ba_problem.num_observations_);
  ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

  ba_problem.num_parameters_ =
    9 * ba_problem.num_cameras_ // #[Rotation|translation|K] = [3x1]|[3x1]|[3x1]
    + 3 * ba_problem.num_points_; // #[X] = [3x1]
  ba_problem.parameters_.reserve(ba_problem.num_parameters_);

  double ppx = 500, ppy = 500;
  // Fill it with data (tracks and points coords)
  for (int i = 0; i < npoints; ++i) {
    // Collect the image of point i in each frame.
    for (int j = 0; j < nviews; ++j) {
      ba_problem.camera_index_.push_back(j);
      ba_problem.point_index_.push_back(i);
      const Vec2 & pt = d._x[j].col(i);
      // => random noise between [-.5,.5] is added
      ba_problem.observations_.push_back( pt(0) - ppx + rand()/RAND_MAX - .5);
      ba_problem.observations_.push_back( pt(1) - ppy + rand()/RAND_MAX - .5);
    }
  }

  // Add camera parameters (R, t, K)
  for (int j = 0; j < nviews; ++j) {
    // Rotation matrix to angle axis
    std::vector<double> angleAxis(3);
    ceres::RotationMatrixToAngleAxis((const double*)d._R[j].data(), &angleAxis[0]);
    // translation
    Vec3 t = d._t[j];
    double focal = d._K[j](0,0);
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
    ba_problem.parameters_.push_back(focal);
    ba_problem.parameters_.push_back(ppx);
    ba_problem.parameters_.push_back(ppy);
  }

  // Add 3D points coordinates parameters
  for (int i = 0; i < npoints; ++i) {
    Vec3 pt3D = d._X.col(i);
    ba_problem.parameters_.push_back(pt3D[0]);
    ba_problem.parameters_.push_back(pt3D[1]);
    ba_problem.parameters_.push_back(pt3D[2]);
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  for (int i = 0; i < ba_problem.num_observations(); ++i) {
    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<pinhole_reprojectionError::ErrorFunc_Refine_Camera_3DPoints, 2, 9, 3>(
            new pinhole_reprojectionError::ErrorFunc_Refine_Camera_3DPoints(
                & ba_problem.observations()[2 * i]));

    problem.AddResidualBlock(cost_function,
                             NULL, // squared loss
                             ba_problem.mutable_camera_for_observation(i),
                             ba_problem.mutable_point_for_observation(i));
  }

  // Make Ceres automatically detect the bundle structure. Note that the
  // standard solver, SPARSE_NORMAL_CHOLESKY, also works fine but it is slower
  // for standard bundle adjustment problems.
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::SPARSE_SCHUR;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
  else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse backend for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << "\n";

  double dResidual_before = std::sqrt( summary.initial_cost / (ba_problem.num_observations_*2.));
  double dResidual_after = std::sqrt( summary.final_cost / (ba_problem.num_observations_*2.));

  std::cout << std::endl
    << " Initial RMSE : " << dResidual_before << "\n"
    << " Final RMSE : " << dResidual_after << "\n"
    << std::endl;

  CHECK(summary.IsSolutionUsable());

  EXPECT_TRUE( dResidual_before > dResidual_after);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
