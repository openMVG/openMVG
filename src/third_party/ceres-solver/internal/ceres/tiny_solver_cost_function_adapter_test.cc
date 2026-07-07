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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/tiny_solver_cost_function_adapter.h"

#include <algorithm>
#include <cmath>
#include <memory>

#include "ceres/cost_function.h"
#include "ceres/sized_cost_function.h"
#include "gtest/gtest.h"

namespace ceres {

class CostFunction2x3 : public SizedCostFunction<2, 3> {
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    double x = parameters[0][0];
    double y = parameters[0][1];
    double z = parameters[0][2];

    residuals[0] = x + 2 * y + 4 * z;
    residuals[1] = y * z;

    if (jacobians && jacobians[0]) {
      jacobians[0][0] = 1;
      jacobians[0][1] = 2;
      jacobians[0][2] = 4;

      jacobians[0][3 + 0] = 0;
      jacobians[0][3 + 1] = z;
      jacobians[0][3 + 2] = y;
    }

    return true;
  }
};

template <int kNumResiduals, int kNumParameters>
void TestHelper() {
  std::unique_ptr<CostFunction> cost_function(new CostFunction2x3);
  using CostFunctionAdapter =
      TinySolverCostFunctionAdapter<kNumResiduals, kNumParameters>;
  CostFunctionAdapter cfa(*cost_function);
  EXPECT_EQ(CostFunctionAdapter::NUM_RESIDUALS, kNumResiduals);
  EXPECT_EQ(CostFunctionAdapter::NUM_PARAMETERS, kNumParameters);

  EXPECT_EQ(cfa.NumResiduals(), 2);
  EXPECT_EQ(cfa.NumParameters(), 3);

  Eigen::Matrix<double, 2, 1> actual_residuals, expected_residuals;
  Eigen::Matrix<double, 2, 3, Eigen::ColMajor> actual_jacobian;
  Eigen::Matrix<double, 2, 3, Eigen::RowMajor> expected_jacobian;

  double xyz[3] = {1.0, -1.0, 2.0};
  double* parameters[1] = {xyz};

  // Check that residual only evaluation works.
  cost_function->Evaluate(parameters, expected_residuals.data(), nullptr);
  cfa(xyz, actual_residuals.data(), nullptr);
  EXPECT_NEAR(
      (expected_residuals - actual_residuals).norm() / actual_residuals.norm(),
      0.0,
      std::numeric_limits<double>::epsilon())
      << "\nExpected residuals: " << expected_residuals.transpose()
      << "\nActual residuals: " << actual_residuals.transpose();

  // Check that residual and jacobian evaluation works.
  double* jacobians[1] = {expected_jacobian.data()};
  cost_function->Evaluate(parameters, expected_residuals.data(), jacobians);
  cfa(xyz, actual_residuals.data(), actual_jacobian.data());

  EXPECT_NEAR(
      (expected_residuals - actual_residuals).norm() / actual_residuals.norm(),
      0.0,
      std::numeric_limits<double>::epsilon())
      << "\nExpected residuals: " << expected_residuals.transpose()
      << "\nActual residuals: " << actual_residuals.transpose();

  EXPECT_NEAR(
      (expected_jacobian - actual_jacobian).norm() / actual_jacobian.norm(),
      0.0,
      std::numeric_limits<double>::epsilon())
      << "\nExpected jacobian: " << expected_jacobian.transpose()
      << "\nActual jacobian: " << actual_jacobian.transpose();
}

TEST(TinySolverCostFunctionAdapter, StaticResidualsStaticParameterBlock) {
  TestHelper<2, 3>();
}

TEST(TinySolverCostFunctionAdapter, DynamicResidualsStaticParameterBlock) {
  TestHelper<Eigen::Dynamic, 3>();
}

TEST(TinySolverCostFunctionAdapter, StaticResidualsDynamicParameterBlock) {
  TestHelper<2, Eigen::Dynamic>();
}

TEST(TinySolverCostFunctionAdapter, DynamicResidualsDynamicParameterBlock) {
  TestHelper<Eigen::Dynamic, Eigen::Dynamic>();
}

}  // namespace ceres
