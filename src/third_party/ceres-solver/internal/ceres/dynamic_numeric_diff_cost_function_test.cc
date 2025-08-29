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
//         mierle@gmail.com (Keir Mierle)

#include "ceres/dynamic_numeric_diff_cost_function.h"

#include <cstddef>
#include <memory>
#include <vector>

#include "gtest/gtest.h"

namespace ceres::internal {

const double kTolerance = 1e-6;

// Takes 2 parameter blocks:
//     parameters[0] is size 10.
//     parameters[1] is size 5.
// Emits 21 residuals:
//     A: i - parameters[0][i], for i in [0,10)  -- this is 10 residuals
//     B: parameters[0][i] - i, for i in [0,10)  -- this is another 10.
//     C: sum(parameters[0][i]^2 - 8*parameters[0][i]) + sum(parameters[1][i])
class MyCostFunctor {
 public:
  bool operator()(double const* const* parameters, double* residuals) const {
    const double* params0 = parameters[0];
    int r = 0;
    for (int i = 0; i < 10; ++i) {
      residuals[r++] = i - params0[i];
      residuals[r++] = params0[i] - i;
    }

    double c_residual = 0.0;
    for (int i = 0; i < 10; ++i) {
      c_residual += pow(params0[i], 2) - 8.0 * params0[i];
    }

    const double* params1 = parameters[1];
    for (int i = 0; i < 5; ++i) {
      c_residual += params1[i];
    }
    residuals[r++] = c_residual;
    return true;
  }
};

TEST(DynamicNumericdiffCostFunctionTest, TestResiduals) {
  std::vector<double> param_block_0(10, 0.0);
  std::vector<double> param_block_1(5, 0.0);
  DynamicNumericDiffCostFunction<MyCostFunctor> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Test residual computation.
  std::vector<double> residuals(21, -100000);
  std::vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];
  EXPECT_TRUE(
      cost_function.Evaluate(&parameter_blocks[0], residuals.data(), nullptr));
  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(0, residuals.at(20));
}

TEST(DynamicNumericdiffCostFunctionTest, TestJacobian) {
  // Test the residual counting.
  std::vector<double> param_block_0(10, 0.0);
  for (int i = 0; i < 10; ++i) {
    param_block_0[i] = 2 * i;
  }
  std::vector<double> param_block_1(5, 0.0);
  DynamicNumericDiffCostFunction<MyCostFunctor> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Prepare the residuals.
  std::vector<double> residuals(21, -100000);

  // Prepare the parameters.
  std::vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];

  // Prepare the jacobian.
  std::vector<std::vector<double>> jacobian_vect(2);
  jacobian_vect[0].resize(21 * 10, -100000);
  jacobian_vect[1].resize(21 * 5, -100000);
  std::vector<double*> jacobian;
  jacobian.push_back(jacobian_vect[0].data());
  jacobian.push_back(jacobian_vect[1].data());

  // Test jacobian computation.
  EXPECT_TRUE(cost_function.Evaluate(
      parameter_blocks.data(), residuals.data(), jacobian.data()));

  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(+1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(420, residuals.at(20));
  for (int p = 0; p < 10; ++p) {
    // Check "A" Jacobian.
    EXPECT_NEAR(-1.0, jacobian_vect[0][2 * p * 10 + p], kTolerance);
    // Check "B" Jacobian.
    EXPECT_NEAR(+1.0, jacobian_vect[0][(2 * p + 1) * 10 + p], kTolerance);
    jacobian_vect[0][2 * p * 10 + p] = 0.0;
    jacobian_vect[0][(2 * p + 1) * 10 + p] = 0.0;
  }

  // Check "C" Jacobian for first parameter block.
  for (int p = 0; p < 10; ++p) {
    EXPECT_NEAR(4 * p - 8, jacobian_vect[0][20 * 10 + p], kTolerance);
    jacobian_vect[0][20 * 10 + p] = 0.0;
  }
  for (double entry : jacobian_vect[0]) {
    EXPECT_NEAR(0.0, entry, kTolerance);
  }

  // Check "C" Jacobian for second parameter block.
  for (int p = 0; p < 5; ++p) {
    EXPECT_NEAR(1.0, jacobian_vect[1][20 * 5 + p], kTolerance);
    jacobian_vect[1][20 * 5 + p] = 0.0;
  }
  for (double entry : jacobian_vect[1]) {
    EXPECT_NEAR(0.0, entry, kTolerance);
  }
}

TEST(DynamicNumericdiffCostFunctionTest,
     JacobianWithFirstParameterBlockConstant) {  // NOLINT
  // Test the residual counting.
  std::vector<double> param_block_0(10, 0.0);
  for (int i = 0; i < 10; ++i) {
    param_block_0[i] = 2 * i;
  }
  std::vector<double> param_block_1(5, 0.0);
  DynamicNumericDiffCostFunction<MyCostFunctor> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Prepare the residuals.
  std::vector<double> residuals(21, -100000);

  // Prepare the parameters.
  std::vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];

  // Prepare the jacobian.
  std::vector<std::vector<double>> jacobian_vect(2);
  jacobian_vect[0].resize(21 * 10, -100000);
  jacobian_vect[1].resize(21 * 5, -100000);
  std::vector<double*> jacobian;
  jacobian.push_back(nullptr);
  jacobian.push_back(jacobian_vect[1].data());

  // Test jacobian computation.
  EXPECT_TRUE(cost_function.Evaluate(
      parameter_blocks.data(), residuals.data(), jacobian.data()));

  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(+1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(420, residuals.at(20));

  // Check "C" Jacobian for second parameter block.
  for (int p = 0; p < 5; ++p) {
    EXPECT_NEAR(1.0, jacobian_vect[1][20 * 5 + p], kTolerance);
    jacobian_vect[1][20 * 5 + p] = 0.0;
  }
  for (double& i : jacobian_vect[1]) {
    EXPECT_EQ(0.0, i);
  }
}

TEST(DynamicNumericdiffCostFunctionTest,
     JacobianWithSecondParameterBlockConstant) {  // NOLINT
  // Test the residual counting.
  std::vector<double> param_block_0(10, 0.0);
  for (int i = 0; i < 10; ++i) {
    param_block_0[i] = 2 * i;
  }
  std::vector<double> param_block_1(5, 0.0);
  DynamicNumericDiffCostFunction<MyCostFunctor> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Prepare the residuals.
  std::vector<double> residuals(21, -100000);

  // Prepare the parameters.
  std::vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];

  // Prepare the jacobian.
  std::vector<std::vector<double>> jacobian_vect(2);
  jacobian_vect[0].resize(21 * 10, -100000);
  jacobian_vect[1].resize(21 * 5, -100000);
  std::vector<double*> jacobian;
  jacobian.push_back(jacobian_vect[0].data());
  jacobian.push_back(nullptr);

  // Test jacobian computation.
  EXPECT_TRUE(cost_function.Evaluate(
      parameter_blocks.data(), residuals.data(), jacobian.data()));

  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(+1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(420, residuals.at(20));
  for (int p = 0; p < 10; ++p) {
    // Check "A" Jacobian.
    EXPECT_NEAR(-1.0, jacobian_vect[0][2 * p * 10 + p], kTolerance);
    // Check "B" Jacobian.
    EXPECT_NEAR(+1.0, jacobian_vect[0][(2 * p + 1) * 10 + p], kTolerance);
    jacobian_vect[0][2 * p * 10 + p] = 0.0;
    jacobian_vect[0][(2 * p + 1) * 10 + p] = 0.0;
  }

  // Check "C" Jacobian for first parameter block.
  for (int p = 0; p < 10; ++p) {
    EXPECT_NEAR(4 * p - 8, jacobian_vect[0][20 * 10 + p], kTolerance);
    jacobian_vect[0][20 * 10 + p] = 0.0;
  }
  for (double& i : jacobian_vect[0]) {
    EXPECT_EQ(0.0, i);
  }
}

// Takes 3 parameter blocks:
//     parameters[0] (x) is size 1.
//     parameters[1] (y) is size 2.
//     parameters[2] (z) is size 3.
// Emits 7 residuals:
//     A: x[0] (= sum_x)
//     B: y[0] + 2.0 * y[1] (= sum_y)
//     C: z[0] + 3.0 * z[1] + 6.0 * z[2] (= sum_z)
//     D: sum_x * sum_y
//     E: sum_y * sum_z
//     F: sum_x * sum_z
//     G: sum_x * sum_y * sum_z
class MyThreeParameterCostFunctor {
 public:
  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const {
    const T* x = parameters[0];
    const T* y = parameters[1];
    const T* z = parameters[2];

    T sum_x = x[0];
    T sum_y = y[0] + 2.0 * y[1];
    T sum_z = z[0] + 3.0 * z[1] + 6.0 * z[2];

    residuals[0] = sum_x;
    residuals[1] = sum_y;
    residuals[2] = sum_z;
    residuals[3] = sum_x * sum_y;
    residuals[4] = sum_y * sum_z;
    residuals[5] = sum_x * sum_z;
    residuals[6] = sum_x * sum_y * sum_z;
    return true;
  }
};

class ThreeParameterCostFunctorTest : public ::testing::Test {
 protected:
  void SetUp() final {
    // Prepare the parameters.
    x_.resize(1);
    x_[0] = 0.0;

    y_.resize(2);
    y_[0] = 1.0;
    y_[1] = 3.0;

    z_.resize(3);
    z_[0] = 2.0;
    z_[1] = 4.0;
    z_[2] = 6.0;

    parameter_blocks_.resize(3);
    parameter_blocks_[0] = &x_[0];
    parameter_blocks_[1] = &y_[0];
    parameter_blocks_[2] = &z_[0];

    // Prepare the cost function.
    using DynamicMyThreeParameterCostFunction =
        DynamicNumericDiffCostFunction<MyThreeParameterCostFunctor>;
    auto cost_function = std::make_unique<DynamicMyThreeParameterCostFunction>(
        new MyThreeParameterCostFunctor());
    cost_function->AddParameterBlock(1);
    cost_function->AddParameterBlock(2);
    cost_function->AddParameterBlock(3);
    cost_function->SetNumResiduals(7);

    cost_function_ = std::move(cost_function);

    // Setup jacobian data.
    jacobian_vect_.resize(3);
    jacobian_vect_[0].resize(7 * x_.size(), -100000);
    jacobian_vect_[1].resize(7 * y_.size(), -100000);
    jacobian_vect_[2].resize(7 * z_.size(), -100000);

    // Prepare the expected residuals.
    const double sum_x = x_[0];
    const double sum_y = y_[0] + 2.0 * y_[1];
    const double sum_z = z_[0] + 3.0 * z_[1] + 6.0 * z_[2];

    expected_residuals_.resize(7);
    expected_residuals_[0] = sum_x;
    expected_residuals_[1] = sum_y;
    expected_residuals_[2] = sum_z;
    expected_residuals_[3] = sum_x * sum_y;
    expected_residuals_[4] = sum_y * sum_z;
    expected_residuals_[5] = sum_x * sum_z;
    expected_residuals_[6] = sum_x * sum_y * sum_z;

    // Prepare the expected jacobian entries.
    expected_jacobian_x_.resize(7);
    expected_jacobian_x_[0] = 1.0;
    expected_jacobian_x_[1] = 0.0;
    expected_jacobian_x_[2] = 0.0;
    expected_jacobian_x_[3] = sum_y;
    expected_jacobian_x_[4] = 0.0;
    expected_jacobian_x_[5] = sum_z;
    expected_jacobian_x_[6] = sum_y * sum_z;

    expected_jacobian_y_.resize(14);
    expected_jacobian_y_[0] = 0.0;
    expected_jacobian_y_[1] = 0.0;
    expected_jacobian_y_[2] = 1.0;
    expected_jacobian_y_[3] = 2.0;
    expected_jacobian_y_[4] = 0.0;
    expected_jacobian_y_[5] = 0.0;
    expected_jacobian_y_[6] = sum_x;
    expected_jacobian_y_[7] = 2.0 * sum_x;
    expected_jacobian_y_[8] = sum_z;
    expected_jacobian_y_[9] = 2.0 * sum_z;
    expected_jacobian_y_[10] = 0.0;
    expected_jacobian_y_[11] = 0.0;
    expected_jacobian_y_[12] = sum_x * sum_z;
    expected_jacobian_y_[13] = 2.0 * sum_x * sum_z;

    expected_jacobian_z_.resize(21);
    expected_jacobian_z_[0] = 0.0;
    expected_jacobian_z_[1] = 0.0;
    expected_jacobian_z_[2] = 0.0;
    expected_jacobian_z_[3] = 0.0;
    expected_jacobian_z_[4] = 0.0;
    expected_jacobian_z_[5] = 0.0;
    expected_jacobian_z_[6] = 1.0;
    expected_jacobian_z_[7] = 3.0;
    expected_jacobian_z_[8] = 6.0;
    expected_jacobian_z_[9] = 0.0;
    expected_jacobian_z_[10] = 0.0;
    expected_jacobian_z_[11] = 0.0;
    expected_jacobian_z_[12] = sum_y;
    expected_jacobian_z_[13] = 3.0 * sum_y;
    expected_jacobian_z_[14] = 6.0 * sum_y;
    expected_jacobian_z_[15] = sum_x;
    expected_jacobian_z_[16] = 3.0 * sum_x;
    expected_jacobian_z_[17] = 6.0 * sum_x;
    expected_jacobian_z_[18] = sum_x * sum_y;
    expected_jacobian_z_[19] = 3.0 * sum_x * sum_y;
    expected_jacobian_z_[20] = 6.0 * sum_x * sum_y;
  }

 protected:
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> z_;

  std::vector<double*> parameter_blocks_;

  std::unique_ptr<CostFunction> cost_function_;

  std::vector<std::vector<double>> jacobian_vect_;

  std::vector<double> expected_residuals_;

  std::vector<double> expected_jacobian_x_;
  std::vector<double> expected_jacobian_y_;
  std::vector<double> expected_jacobian_z_;
};

TEST_F(ThreeParameterCostFunctorTest, TestThreeParameterResiduals) {
  std::vector<double> residuals(7, -100000);
  EXPECT_TRUE(cost_function_->Evaluate(
      parameter_blocks_.data(), residuals.data(), nullptr));
  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }
}

TEST_F(ThreeParameterCostFunctorTest, TestThreeParameterJacobian) {
  std::vector<double> residuals(7, -100000);

  std::vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(jacobian_vect_[1].data());
  jacobian.push_back(jacobian_vect_[2].data());

  EXPECT_TRUE(cost_function_->Evaluate(
      parameter_blocks_.data(), residuals.data(), jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 7; ++i) {
    EXPECT_NEAR(expected_jacobian_x_[i], jacobian[0][i], kTolerance);
  }

  for (int i = 0; i < 14; ++i) {
    EXPECT_NEAR(expected_jacobian_y_[i], jacobian[1][i], kTolerance);
  }

  for (int i = 0; i < 21; ++i) {
    EXPECT_NEAR(expected_jacobian_z_[i], jacobian[2][i], kTolerance);
  }
}

TEST_F(ThreeParameterCostFunctorTest,
       ThreeParameterJacobianWithFirstAndLastParameterBlockConstant) {
  std::vector<double> residuals(7, -100000);

  std::vector<double*> jacobian;
  jacobian.push_back(nullptr);
  jacobian.push_back(jacobian_vect_[1].data());
  jacobian.push_back(nullptr);

  EXPECT_TRUE(cost_function_->Evaluate(
      parameter_blocks_.data(), residuals.data(), jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 14; ++i) {
    EXPECT_NEAR(expected_jacobian_y_[i], jacobian[1][i], kTolerance);
  }
}

TEST_F(ThreeParameterCostFunctorTest,
       ThreeParameterJacobianWithSecondParameterBlockConstant) {
  std::vector<double> residuals(7, -100000);

  std::vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(nullptr);
  jacobian.push_back(jacobian_vect_[2].data());

  EXPECT_TRUE(cost_function_->Evaluate(
      parameter_blocks_.data(), residuals.data(), jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 7; ++i) {
    EXPECT_NEAR(expected_jacobian_x_[i], jacobian[0][i], kTolerance);
  }

  for (int i = 0; i < 21; ++i) {
    EXPECT_NEAR(expected_jacobian_z_[i], jacobian[2][i], kTolerance);
  }
}

}  // namespace ceres::internal
