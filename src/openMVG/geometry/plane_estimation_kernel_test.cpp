// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/plane_estimation_kernel.hpp"
#include "openMVG/numeric/numeric.h"
#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"
#include <random>

using namespace openMVG;
using namespace openMVG::geometry;
using namespace openMVG::geometry::plane;

TEST(Kernel, SolverAndError)
{
  halfPlane::Half_planes half_plane;
  Mat3X points(3, 3);
  points << 0, 5, 0,
            0, 0, 5,
            0, 0, 0;
  PlaneSolver::Solve(points, &half_plane);

  EXPECT_FALSE(half_plane.empty());
  const Vec3 gt_normal(Vec3::UnitZ());
  EXPECT_NEAR(0.0, gt_normal.cross(half_plane[0].normal()).lpNorm<Eigen::Infinity>(), 1e-8);

  // Distance error check:
  for (int i : {0, 1, 2})
  {
    EXPECT_NEAR(0.0, AbsDistanceError::Error(half_plane[0], points.col(i)), 1e-8);
  }
  // Distance error check:
  EXPECT_NEAR(10.0, AbsDistanceError::Error(half_plane[0], Vec3(0, 0,  10)), 1e-8);
  EXPECT_NEAR(10.0, AbsDistanceError::Error(half_plane[0], Vec3(0, 0, -10)), 1e-8);
}

//
// Testing robust estimation of a plane thanks to LMeds
//

#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"

/// Generate a grid of points on the x, y plane (z = 0)
Mat3X GeneratePointGrid
(
  const int subdivision = 4
)
{
  Mat3X points(3, subdivision * subdivision);
  int cpt = 0;
  for (int x = 0; x < subdivision; ++x)
    for (int y = 0; y < subdivision; ++y)
      points.col(cpt++) << x, y, 0;
  return points;
}

TEST(Kernel, RotatedGrid_RobustEstimation_Distance_Metric)
{
  // Generate point on a grid
  const int subdivision = 4;
  Mat3X points = GeneratePointGrid(subdivision);

  // Add 5 outliers
  {
    std::mt19937 gen(std::mt19937::default_seed);
    std::uniform_int_distribution<> dis(4, 10);
    for (int i : {0,1,2,3,4})
      points.col(i).z() += dis(gen);
  }

  // Apply a known rotation
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_real_distribution<> dis(0, M_PI);
  const Mat3 rotation = RotationAroundX(dis(gen)) * RotationAroundY(dis(gen));
  points = rotation * points;

  // Robust estimation of plane according an point to plane distance parametrization
  {
    HaflPlaneKernel kernel(points);
    HaflPlaneKernel::Model half_plane;
    double dThreshold = std::numeric_limits<double>::infinity();
    const double dBestMedian = robust::LeastMedianOfSquares(kernel, &half_plane, &dThreshold);
    // Check that LMedS found a viable solution
    EXPECT_NEAR(0.0, dThreshold, 1e-8);
    const Vec3 gt_normal(rotation * Vec3::UnitZ());
    EXPECT_NEAR(0.0, gt_normal.cross(half_plane.normal()).lpNorm<Eigen::Infinity>(), 1e-8);

    // Check that distance computation are correct
    for (int i = 0; i < points.cols(); ++i)
    {
      if (i < 5)
      {
        EXPECT_TRUE(
          AbsDistanceError::Error(half_plane, points.col(i)) +
          std::numeric_limits<float>::epsilon() >= 4.0);
      }
      else
      {
        EXPECT_NEAR(0.0, AbsDistanceError::Error(half_plane, points.col(i)), 1e-8);
      }
    }
  }
}

TEST(Kernel, RotatedGrid_RobustEstimation_Angular_Metric)
{
  // Generate point on a grid
  const int subdivision = 4;
  Mat3X points = GeneratePointGrid(subdivision);

  // Add 5 outliers
  {
    std::mt19937 gen(std::mt19937::default_seed);
    std::uniform_int_distribution<> dis(4, 10);
    for (int i : {0,1,2,3,4})
      points.col(i).z() += dis(gen);
  }

  // Apply a known rotation
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_real_distribution<> dis(0, M_PI);
  const Mat3 rotation = RotationAroundX(dis(gen)) * RotationAroundY(dis(gen));
  points = rotation * points;

  // Robust estimation of plane according an angular parametrization
  {
    HaflPlaneKernelAngular kernel(points);
    HaflPlaneKernelAngular::Model half_plane;
    double dThreshold = std::numeric_limits<double>::infinity();
    const double dBestMedian = robust::LeastMedianOfSquares(kernel, &half_plane, &dThreshold);
    // Check that LMedS found a viable solution
    EXPECT_NEAR(0.0, dThreshold, 1e-8);
    const Vec3 gt_normal(rotation * Vec3::UnitZ());
    EXPECT_NEAR(0.0, gt_normal.cross(half_plane.first.normal()).lpNorm<Eigen::Infinity>(), 1e-8);

    // Check that distance computation are correct
    for (int i = 0; i < points.cols(); ++i)
    {
      if (i < 5)
      {
        EXPECT_TRUE(AbsAngularError::Error(half_plane, points.col(i)) >= M_PI / 4.0);
      }
      else
      {
        EXPECT_NEAR(0.0, AbsAngularError::Error(half_plane, points.col(i)), 1e-8);
      }
    }
  }
}

//
// Testing robust estimation of a plane thanks to ACRANSAC
//

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"

using namespace openMVG::robust;

/// ACRansac Kernel adaptor for plane estimation
template <typename SolverArg,
  typename ErrorArg,
  typename ModelArg >
class ACRANSACOneViewKernel
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;

  ACRANSACOneViewKernel(const Mat3X &X): X_(X), logalpha0_(log10(1.0/2.0))
  {
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<Model> *models
  )
  const
  {
    const Mat sampled_xs = ExtractColumns(X_, samples);
    Solver::Solve(sampled_xs, models);
  }

  double Error
  (
    uint32_t sample,
    const Model &model
  )
  const
  {
    return ErrorArg::Error(model, X_.col(sample));
  }

  void Errors
  (
    const Model &model,
    std::vector<double> & vec_errors
  )
  const
  {
    for (size_t sample = 0; sample < X_.cols(); ++sample)
      vec_errors[sample] = Square(ErrorArg::Error(model, X_.col(sample)));
  }

  size_t NumSamples() const {return X_.cols();}

  void Unnormalize(Model * model) const { // Model is left unchanged
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1./4.;}

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const {return sqrt(val);}

private:
  Mat X_;
  double logalpha0_;
};

TEST(Kernel, RotatedGrid_RobustEstimation_Angular_Metric_ACRANSAC)
{
  // Generate point on a grid
  const int subdivision = 4;
  Mat3X points = GeneratePointGrid(subdivision);

  // Add 5 outliers
  {
    std::mt19937 gen(std::mt19937::default_seed);
    std::uniform_int_distribution<> dis(4, 10);
    for (int i : {0,1,2,3,4})
      points.col(i).z() += dis(gen);
  }

  // Apply a known rotation
  std::mt19937 gen(std::mt19937::default_seed);
  std::uniform_real_distribution<> dis(0, M_PI);
  const Mat3 rotation = RotationAroundX(dis(gen)) * RotationAroundY(dis(gen));
  points = rotation * points;

  // Robust estimation of plane according an angular parametrization
  {
    ACRANSACOneViewKernel<
      PlaneSolver,
      AbsAngularError,
      std::pair<halfPlane::Half_plane, Vec3>> kernel(points);

    std::vector<uint32_t> vec_inliers;
    std::pair<halfPlane::Half_plane, Vec3> half_plane;
    ACRANSAC(kernel, vec_inliers, 300, &half_plane);

    for (int index : {0,1,2,3,4})
      EXPECT_EQ(0,
                std::count(vec_inliers.cbegin(), vec_inliers.cend(), index));
    // Check that ACRANSAC found a viable solution
    const Vec3 gt_normal(rotation * Vec3::UnitZ());
    EXPECT_NEAR(0.0, gt_normal.cross(half_plane.first.normal()).lpNorm<Eigen::Infinity>(), 1e-8);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
