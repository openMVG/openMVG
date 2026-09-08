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

#include "ceres/numeric_diff_first_order_function.h"

#include <memory>

#include "ceres/array_utils.h"
#include "ceres/first_order_function.h"
#include "gtest/gtest.h"

namespace ceres::internal {

class QuadraticCostFunctor {
 public:
  explicit QuadraticCostFunctor(double a) : a_(a) {}
  bool operator()(const double* const x, double* cost) const {
    cost[0] = x[0] * x[1] + x[2] * x[3] - a_;
    return true;
  }

 private:
  double a_;
};

TEST(NumericDiffFirstOrderFunction, BilinearDifferentiationTestStatic) {
  auto function = std::make_unique<
      NumericDiffFirstOrderFunction<QuadraticCostFunctor, CENTRAL, 4>>(
      new QuadraticCostFunctor(1.0));

  double parameters[4] = {1.0, 2.0, 3.0, 4.0};
  double gradient[4];
  double cost;

  function->Evaluate(parameters, &cost, nullptr);
  EXPECT_EQ(cost, 13.0);

  cost = -1.0;
  function->Evaluate(parameters, &cost, gradient);

  EXPECT_EQ(cost, 13.0);

  const double kTolerance = 1e-9;
  EXPECT_NEAR(gradient[0], parameters[1], kTolerance);
  EXPECT_NEAR(gradient[1], parameters[0], kTolerance);
  EXPECT_NEAR(gradient[2], parameters[3], kTolerance);
  EXPECT_NEAR(gradient[3], parameters[2], kTolerance);
}

TEST(NumericDiffFirstOrderFunction, BilinearDifferentiationTestDynamic) {
  auto function = std::make_unique<
      NumericDiffFirstOrderFunction<QuadraticCostFunctor, CENTRAL>>(
      new QuadraticCostFunctor(1.0), 4);

  double parameters[4] = {1.0, 2.0, 3.0, 4.0};
  double gradient[4];
  double cost;

  function->Evaluate(parameters, &cost, nullptr);
  EXPECT_EQ(cost, 13.0);

  cost = -1.0;
  function->Evaluate(parameters, &cost, gradient);

  EXPECT_EQ(cost, 13.0);

  const double kTolerance = 1e-9;
  EXPECT_NEAR(gradient[0], parameters[1], kTolerance);
  EXPECT_NEAR(gradient[1], parameters[0], kTolerance);
  EXPECT_NEAR(gradient[2], parameters[3], kTolerance);
  EXPECT_NEAR(gradient[3], parameters[2], kTolerance);
}

}  // namespace ceres::internal
