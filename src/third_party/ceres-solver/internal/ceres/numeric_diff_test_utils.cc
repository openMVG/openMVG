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
//         tbennun@gmail.com (Tal Ben-Nun)

#include "ceres/numeric_diff_test_utils.h"

#include <algorithm>
#include <cmath>

#include "ceres/cost_function.h"
#include "ceres/test_util.h"
#include "ceres/types.h"
#include "gtest/gtest.h"

namespace ceres::internal {

bool EasyFunctor::operator()(const double* x1,
                             const double* x2,
                             double* residuals) const {
  residuals[0] = residuals[1] = residuals[2] = 0;
  for (int i = 0; i < 5; ++i) {
    residuals[0] += x1[i] * x2[i];
    residuals[2] += x2[i] * x2[i];
  }
  residuals[1] = residuals[0] * residuals[0];
  return true;
}

void EasyFunctor::ExpectCostFunctionEvaluationIsNearlyCorrect(
    const CostFunction& cost_function, NumericDiffMethodType method) const {
  // The x1[0] is made deliberately small to test the performance near zero.
  // clang-format off
  double x1[] = { 1e-64, 2.0, 3.0, 4.0, 5.0 };
  double x2[] = { 9.0, 9.0, 5.0, 5.0, 1.0 };
  double *parameters[] = { &x1[0], &x2[0] };
  // clang-format on

  double dydx1[15];  // 3 x 5, row major.
  double dydx2[15];  // 3 x 5, row major.
  double* jacobians[2] = {&dydx1[0], &dydx2[0]};

  double residuals[3] = {-1e-100, -2e-100, -3e-100};

  ASSERT_TRUE(
      cost_function.Evaluate(&parameters[0], &residuals[0], &jacobians[0]));

  double expected_residuals[3];
  EasyFunctor functor;
  functor(x1, x2, expected_residuals);
  EXPECT_EQ(expected_residuals[0], residuals[0]);
  EXPECT_EQ(expected_residuals[1], residuals[1]);
  EXPECT_EQ(expected_residuals[2], residuals[2]);

  double tolerance = 0.0;
  switch (method) {
    default:
    case CENTRAL:
      tolerance = 3e-9;
      break;

    case FORWARD:
      tolerance = 2e-5;
      break;

    case RIDDERS:
      tolerance = 1e-13;
      break;
  }

  for (int i = 0; i < 5; ++i) {
    // clang-format off
    ExpectClose(x2[i],                    dydx1[5 * 0 + i], tolerance);  // y1
    ExpectClose(x1[i],                    dydx2[5 * 0 + i], tolerance);
    ExpectClose(2 * x2[i] * residuals[0], dydx1[5 * 1 + i], tolerance);  // y2
    ExpectClose(2 * x1[i] * residuals[0], dydx2[5 * 1 + i], tolerance);
    ExpectClose(0.0,                      dydx1[5 * 2 + i], tolerance);  // y3
    ExpectClose(2 * x2[i],                dydx2[5 * 2 + i], tolerance);
    // clang-format on
  }
}

bool TranscendentalFunctor::operator()(const double* x1,
                                       const double* x2,
                                       double* residuals) const {
  double x1x2 = 0;
  for (int i = 0; i < 5; ++i) {
    x1x2 += x1[i] * x2[i];
  }
  residuals[0] = sin(x1x2);
  residuals[1] = exp(-x1x2 / 10);
  return true;
}

void TranscendentalFunctor::ExpectCostFunctionEvaluationIsNearlyCorrect(
    const CostFunction& cost_function, NumericDiffMethodType method) const {
  struct TestParameterBlocks {
    double x1[5];
    double x2[5];
  };

  // clang-format off
  std::vector<TestParameterBlocks> kTests =  {
    { { 1.0, 2.0, 3.0, 4.0, 5.0 },  // No zeros.
      { 9.0, 9.0, 5.0, 5.0, 1.0 },
    },
    { { 0.0, 2.0, 3.0, 0.0, 5.0 },  // Some zeros x1.
      { 9.0, 9.0, 5.0, 5.0, 1.0 },
    },
    { { 1.0, 2.0, 3.0, 1.0, 5.0 },  // Some zeros x2.
      { 0.0, 9.0, 0.0, 5.0, 0.0 },
    },
    { { 0.0, 0.0, 0.0, 0.0, 0.0 },  // All zeros x1.
      { 9.0, 9.0, 5.0, 5.0, 1.0 },
    },
    { { 1.0, 2.0, 3.0, 4.0, 5.0 },  // All zeros x2.
      { 0.0, 0.0, 0.0, 0.0, 0.0 },
    },
    { { 0.0, 0.0, 0.0, 0.0, 0.0 },  // All zeros.
      { 0.0, 0.0, 0.0, 0.0, 0.0 },
    },
  };
  // clang-format on

  for (auto& test : kTests) {
    double* x1 = &(test.x1[0]);
    double* x2 = &(test.x2[0]);
    double* parameters[] = {x1, x2};

    double dydx1[10];
    double dydx2[10];
    double* jacobians[2] = {&dydx1[0], &dydx2[0]};

    double residuals[2];

    ASSERT_TRUE(
        cost_function.Evaluate(&parameters[0], &residuals[0], &jacobians[0]));
    double x1x2 = 0;
    for (int i = 0; i < 5; ++i) {
      x1x2 += x1[i] * x2[i];
    }

    double tolerance = 0.0;
    switch (method) {
      default:
      case CENTRAL:
        tolerance = 2e-7;
        break;

      case FORWARD:
        tolerance = 2e-5;
        break;

      case RIDDERS:
        tolerance = 3e-12;
        break;
    }

    for (int i = 0; i < 5; ++i) {
      // clang-format off
      ExpectClose( x2[i] * cos(x1x2),              dydx1[5 * 0 + i], tolerance);
      ExpectClose( x1[i] * cos(x1x2),              dydx2[5 * 0 + i], tolerance);
      ExpectClose(-x2[i] * exp(-x1x2 / 10.) / 10., dydx1[5 * 1 + i], tolerance);
      ExpectClose(-x1[i] * exp(-x1x2 / 10.) / 10., dydx2[5 * 1 + i], tolerance);
      // clang-format on
    }
  }
}

bool ExponentialFunctor::operator()(const double* x1, double* residuals) const {
  residuals[0] = exp(x1[0]);
  return true;
}

void ExponentialFunctor::ExpectCostFunctionEvaluationIsNearlyCorrect(
    const CostFunction& cost_function) const {
  // Evaluating the functor at specific points for testing.
  std::vector<double> kTests = {1.0, 2.0, 3.0, 4.0, 5.0};

  // Minimal tolerance w.r.t. the cost function and the tests.
  const double kTolerance = 2e-14;

  for (double& test : kTests) {
    double* parameters[] = {&test};
    double dydx;
    double* jacobians[1] = {&dydx};
    double residual;

    ASSERT_TRUE(
        cost_function.Evaluate(&parameters[0], &residual, &jacobians[0]));

    double expected_result = exp(test);

    // Expect residual to be close to exp(x).
    ExpectClose(residual, expected_result, kTolerance);

    // Check evaluated differences. dydx should also be close to exp(x).
    ExpectClose(dydx, expected_result, kTolerance);
  }
}

bool RandomizedFunctor::operator()(const double* x1, double* residuals) const {
  double random_value = uniform_distribution_(*prng_);
  residuals[0] = x1[0] * x1[0] + random_value;
  return true;
}

void RandomizedFunctor::ExpectCostFunctionEvaluationIsNearlyCorrect(
    const CostFunction& cost_function) const {
  std::vector<double> kTests = {0.0, 1.0, 3.0, 4.0, 50.0};

  const double kTolerance = 2e-4;

  for (double& test : kTests) {
    double* parameters[] = {&test};
    double dydx;
    double* jacobians[1] = {&dydx};
    double residual;

    ASSERT_TRUE(
        cost_function.Evaluate(&parameters[0], &residual, &jacobians[0]));

    // Expect residual to be close to x^2 w.r.t. noise factor.
    ExpectClose(residual, test * test, noise_factor_);

    // Check evaluated differences. (dy/dx = ~2x)
    ExpectClose(dydx, 2 * test, kTolerance);
  }
}

}  // namespace ceres::internal
