
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
// Author: mierle@gmail.com (Keir Mierle)

#include "ceres/tiny_solver_autodiff_function.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "ceres/tiny_solver.h"
#include "ceres/tiny_solver_test_util.h"
#include "gtest/gtest.h"

namespace ceres {

struct AutoDiffTestFunctor {
  template <typename T>
  bool operator()(const T* const parameters, T* residuals) const {
    // Shift the parameters so the solution is not at the origin, to prevent
    // accidentally showing "PASS".
    const T& a = parameters[0] - T(1.0);
    const T& b = parameters[1] - T(2.0);
    const T& c = parameters[2] - T(3.0);
    residuals[0] = 2. * a + 0. * b + 1. * c;
    residuals[1] = 0. * a + 4. * b + 6. * c;
    return true;
  }
};

// Leave a factor of 10 slop since these tests tend to mysteriously break on
// other compilers or architectures if the tolerance is too tight.
static double const kTolerance = std::numeric_limits<double>::epsilon() * 10;

TEST(TinySolverAutoDiffFunction, SimpleFunction) {
  using AutoDiffTestFunction =
      TinySolverAutoDiffFunction<AutoDiffTestFunctor, 2, 3>;
  AutoDiffTestFunctor autodiff_test_functor;
  AutoDiffTestFunction f(autodiff_test_functor);

  Eigen::Vector3d x(2.0, 1.0, 4.0);
  Eigen::Vector2d residuals;

  // Check the case with cost-only evaluation.
  residuals.setConstant(555);  // Arbitrary.
  EXPECT_TRUE(f(&x(0), &residuals(0), nullptr));
  EXPECT_NEAR(3.0, residuals(0), kTolerance);
  EXPECT_NEAR(2.0, residuals(1), kTolerance);

  // Check the case with cost and Jacobian evaluation.
  Eigen::Matrix<double, 2, 3> jacobian;
  residuals.setConstant(555);  // Arbitrary.
  jacobian.setConstant(555);
  EXPECT_TRUE(f(&x(0), &residuals(0), &jacobian(0, 0)));

  // Verify cost.
  EXPECT_NEAR(3.0, residuals(0), kTolerance);
  EXPECT_NEAR(2.0, residuals(1), kTolerance);

  // Verify Jacobian Row 1.
  EXPECT_NEAR(2.0, jacobian(0, 0), kTolerance);
  EXPECT_NEAR(0.0, jacobian(0, 1), kTolerance);
  EXPECT_NEAR(1.0, jacobian(0, 2), kTolerance);

  // Verify Jacobian row 2.
  EXPECT_NEAR(0.0, jacobian(1, 0), kTolerance);
  EXPECT_NEAR(4.0, jacobian(1, 1), kTolerance);
  EXPECT_NEAR(6.0, jacobian(1, 2), kTolerance);
}

class DynamicResidualsFunctor {
 public:
  using Scalar = double;
  enum {
    NUM_RESIDUALS = Eigen::Dynamic,
    NUM_PARAMETERS = 3,
  };

  int NumResiduals() const { return 2; }

  template <typename T>
  bool operator()(const T* parameters, T* residuals) const {
    // Jacobian is not evaluated by cost function, but by autodiff.
    T* jacobian = nullptr;
    return EvaluateResidualsAndJacobians(parameters, residuals, jacobian);
  }
};

template <typename Function, typename Vector>
void TestHelper(const Function& f, const Vector& x0) {
  Vector x = x0;
  Eigen::Vector2d residuals;
  f(x.data(), residuals.data(), nullptr);
  EXPECT_GT(residuals.squaredNorm() / 2.0, 1e-10);

  TinySolver<Function> solver;
  solver.Solve(f, &x);
  EXPECT_NEAR(0.0, solver.summary.final_cost, 1e-10);
}

// A test case for when the number of residuals is
// dynamically sized and we use autodiff
TEST(TinySolverAutoDiffFunction, ResidualsDynamicAutoDiff) {
  Eigen::Vector3d x0(0.76026643, -30.01799744, 0.55192142);

  DynamicResidualsFunctor f;
  using AutoDiffCostFunctor = ceres::
      TinySolverAutoDiffFunction<DynamicResidualsFunctor, Eigen::Dynamic, 3>;
  AutoDiffCostFunctor f_autodiff(f);

  Eigen::Vector2d residuals;
  f_autodiff(x0.data(), residuals.data(), nullptr);
  EXPECT_GT(residuals.squaredNorm() / 2.0, 1e-10);

  TinySolver<AutoDiffCostFunctor> solver;
  solver.Solve(f_autodiff, &x0);
  EXPECT_NEAR(0.0, solver.summary.final_cost, 1e-10);
}

}  // namespace ceres
