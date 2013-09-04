// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2012 Google Inc. All rights reserved.
// http://code.google.com/p/ceres-solver/
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
// Author: thadh@gmail.com (Thad Hughes)
//         mierle@gmail.com (Keir Mierle)
//         sameeragarwal@google.com (Sameer Agarwal)

#include <cstddef>

#include "ceres/dynamic_autodiff_cost_function.h"
#include "ceres/internal/scoped_ptr.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

// Takes 2 parameter blocks:
//     parameters[0] is size 10.
//     parameters[1] is size 5.
// Emits 21 residuals:
//     A: i - parameters[0][i], for i in [0,10)  -- this is 10 residuals
//     B: parameters[0][i] - i, for i in [0,10)  -- this is another 10.
//     C: sum(parameters[0][i]^2 - 8*parameters[0][i]) + sum(parameters[1][i])
class MyCostFunctor {
 public:
  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const {
    const T* params0 = parameters[0];
    int r = 0;
    for (int i = 0; i < 10; ++i) {
      residuals[r++] = T(i) - params0[i];
      residuals[r++] = params0[i] - T(i);
    }

    T c_residual(0.0);
    for (int i = 0; i < 10; ++i) {
      c_residual += pow(params0[i], 2) - T(8) * params0[i];
    }

    const T* params1 = parameters[1];
    for (int i = 0; i < 5; ++i) {
      c_residual += params1[i];
    }
    residuals[r++] = c_residual;
    return true;
  }
};

TEST(DynamicAutodiffCostFunctionTest, TestResiduals) {
  vector<double> param_block_0(10, 0.0);
  vector<double> param_block_1(5, 0.0);
  DynamicAutoDiffCostFunction<MyCostFunctor, 3> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Test residual computation.
  vector<double> residuals(21, -100000);
  vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];
  EXPECT_TRUE(cost_function.Evaluate(&parameter_blocks[0],
                                     residuals.data(),
                                     NULL));
  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(0, residuals.at(20));
}

TEST(DynamicAutodiffCostFunctionTest, TestJacobian) {
  // Test the residual counting.
  vector<double> param_block_0(10, 0.0);
  for (int i = 0; i < 10; ++i) {
    param_block_0[i] = 2 * i;
  }
  vector<double> param_block_1(5, 0.0);
  DynamicAutoDiffCostFunction<MyCostFunctor, 3> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Prepare the residuals.
  vector<double> residuals(21, -100000);

  // Prepare the parameters.
  vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];

  // Prepare the jacobian.
  vector<vector<double> > jacobian_vect(2);
  jacobian_vect[0].resize(21 * 10, -100000);
  jacobian_vect[1].resize(21 * 5, -100000);
  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect[0].data());
  jacobian.push_back(jacobian_vect[1].data());

  // Test jacobian computation.
  EXPECT_TRUE(cost_function.Evaluate(parameter_blocks.data(),
                                     residuals.data(),
                                     jacobian.data()));

  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(+1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(420, residuals.at(20));
  for (int p = 0; p < 10; ++p) {
    // Check "A" Jacobian.
    EXPECT_EQ(-1.0, jacobian_vect[0][2*p * 10 + p]);
    // Check "B" Jacobian.
    EXPECT_EQ(+1.0, jacobian_vect[0][(2*p+1) * 10 + p]);
    jacobian_vect[0][2*p * 10 + p] = 0.0;
    jacobian_vect[0][(2*p+1) * 10 + p] = 0.0;
  }

  // Check "C" Jacobian for first parameter block.
  for (int p = 0; p < 10; ++p) {
    EXPECT_EQ(4 * p - 8, jacobian_vect[0][20 * 10 + p]);
    jacobian_vect[0][20 * 10 + p] = 0.0;
  }
  for (int i = 0; i < jacobian_vect[0].size(); ++i) {
    EXPECT_EQ(0.0, jacobian_vect[0][i]);
  }

  // Check "C" Jacobian for second parameter block.
  for (int p = 0; p < 5; ++p) {
    EXPECT_EQ(1.0, jacobian_vect[1][20 * 5 + p]);
    jacobian_vect[1][20 * 5 + p] = 0.0;
  }
  for (int i = 0; i < jacobian_vect[1].size(); ++i) {
    EXPECT_EQ(0.0, jacobian_vect[1][i]);
  }
}

TEST(DynamicAutodiffCostFunctionTest, JacobianWithFirstParameterBlockConstant) {
  // Test the residual counting.
  vector<double> param_block_0(10, 0.0);
  for (int i = 0; i < 10; ++i) {
    param_block_0[i] = 2 * i;
  }
  vector<double> param_block_1(5, 0.0);
  DynamicAutoDiffCostFunction<MyCostFunctor, 3> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Prepare the residuals.
  vector<double> residuals(21, -100000);

  // Prepare the parameters.
  vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];

  // Prepare the jacobian.
  vector<vector<double> > jacobian_vect(2);
  jacobian_vect[0].resize(21 * 10, -100000);
  jacobian_vect[1].resize(21 * 5, -100000);
  vector<double*> jacobian;
  jacobian.push_back(NULL);
  jacobian.push_back(jacobian_vect[1].data());

  // Test jacobian computation.
  EXPECT_TRUE(cost_function.Evaluate(parameter_blocks.data(),
                                     residuals.data(),
                                     jacobian.data()));

  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(+1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(420, residuals.at(20));

  // Check "C" Jacobian for second parameter block.
  for (int p = 0; p < 5; ++p) {
    EXPECT_EQ(1.0, jacobian_vect[1][20 * 5 + p]);
    jacobian_vect[1][20 * 5 + p] = 0.0;
  }
  for (int i = 0; i < jacobian_vect[1].size(); ++i) {
    EXPECT_EQ(0.0, jacobian_vect[1][i]);
  }
}

TEST(DynamicAutodiffCostFunctionTest, JacobianWithSecondParameterBlockConstant) {
  // Test the residual counting.
  vector<double> param_block_0(10, 0.0);
  for (int i = 0; i < 10; ++i) {
    param_block_0[i] = 2 * i;
  }
  vector<double> param_block_1(5, 0.0);
  DynamicAutoDiffCostFunction<MyCostFunctor, 3> cost_function(
      new MyCostFunctor());
  cost_function.AddParameterBlock(param_block_0.size());
  cost_function.AddParameterBlock(param_block_1.size());
  cost_function.SetNumResiduals(21);

  // Prepare the residuals.
  vector<double> residuals(21, -100000);

  // Prepare the parameters.
  vector<double*> parameter_blocks(2);
  parameter_blocks[0] = &param_block_0[0];
  parameter_blocks[1] = &param_block_1[0];

  // Prepare the jacobian.
  vector<vector<double> > jacobian_vect(2);
  jacobian_vect[0].resize(21 * 10, -100000);
  jacobian_vect[1].resize(21 * 5, -100000);
  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect[0].data());
  jacobian.push_back(NULL);

  // Test jacobian computation.
  EXPECT_TRUE(cost_function.Evaluate(parameter_blocks.data(),
                                     residuals.data(),
                                     jacobian.data()));

  for (int r = 0; r < 10; ++r) {
    EXPECT_EQ(-1.0 * r, residuals.at(r * 2));
    EXPECT_EQ(+1.0 * r, residuals.at(r * 2 + 1));
  }
  EXPECT_EQ(420, residuals.at(20));
  for (int p = 0; p < 10; ++p) {
    // Check "A" Jacobian.
    EXPECT_EQ(-1.0, jacobian_vect[0][2*p * 10 + p]);
    // Check "B" Jacobian.
    EXPECT_EQ(+1.0, jacobian_vect[0][(2*p+1) * 10 + p]);
    jacobian_vect[0][2*p * 10 + p] = 0.0;
    jacobian_vect[0][(2*p+1) * 10 + p] = 0.0;
  }

  // Check "C" Jacobian for first parameter block.
  for (int p = 0; p < 10; ++p) {
    EXPECT_EQ(4 * p - 8, jacobian_vect[0][20 * 10 + p]);
    jacobian_vect[0][20 * 10 + p] = 0.0;
  }
  for (int i = 0; i < jacobian_vect[0].size(); ++i) {
    EXPECT_EQ(0.0, jacobian_vect[0][i]);
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
  virtual void SetUp() {
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
    typedef DynamicAutoDiffCostFunction<MyThreeParameterCostFunctor, 3>
      DynamicMyThreeParameterCostFunction;
    DynamicMyThreeParameterCostFunction * cost_function =
      new DynamicMyThreeParameterCostFunction(
        new MyThreeParameterCostFunctor());
    cost_function->AddParameterBlock(1);
    cost_function->AddParameterBlock(2);
    cost_function->AddParameterBlock(3);
    cost_function->SetNumResiduals(7);

    cost_function_.reset(cost_function);

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
  vector<double> x_;
  vector<double> y_;
  vector<double> z_;

  vector<double*> parameter_blocks_;

  scoped_ptr<CostFunction> cost_function_;

  vector<vector<double> > jacobian_vect_;

  vector<double> expected_residuals_;

  vector<double> expected_jacobian_x_;
  vector<double> expected_jacobian_y_;
  vector<double> expected_jacobian_z_;
};

TEST_F(ThreeParameterCostFunctorTest, TestThreeParameterResiduals) {
  vector<double> residuals(7, -100000);
  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       NULL));
  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }
}

TEST_F(ThreeParameterCostFunctorTest, TestThreeParameterJacobian) {
  vector<double> residuals(7, -100000);

  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(jacobian_vect_[1].data());
  jacobian.push_back(jacobian_vect_[2].data());

  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_jacobian_x_[i], jacobian[0][i]);
  }

  for (int i = 0; i < 14; ++i) {
    EXPECT_EQ(expected_jacobian_y_[i], jacobian[1][i]);
  }

  for (int i = 0; i < 21; ++i) {
    EXPECT_EQ(expected_jacobian_z_[i], jacobian[2][i]);
  }
}

TEST_F(ThreeParameterCostFunctorTest,
       ThreeParameterJacobianWithFirstAndLastParameterBlockConstant) {
  vector<double> residuals(7, -100000);

  vector<double*> jacobian;
  jacobian.push_back(NULL);
  jacobian.push_back(jacobian_vect_[1].data());
  jacobian.push_back(NULL);

  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 14; ++i) {
    EXPECT_EQ(expected_jacobian_y_[i], jacobian[1][i]);
  }
}

TEST_F(ThreeParameterCostFunctorTest,
       ThreeParameterJacobianWithSecondParameterBlockConstant) {
  vector<double> residuals(7, -100000);

  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(NULL);
  jacobian.push_back(jacobian_vect_[2].data());

  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_jacobian_x_[i], jacobian[0][i]);
  }

  for (int i = 0; i < 21; ++i) {
    EXPECT_EQ(expected_jacobian_z_[i], jacobian[2][i]);
  }
}

// Takes 6 parameter blocks all of size 1:
//     x0, y0, y1, z0, z1, z2
// Same 7 residuals as MyThreeParameterCostFunctor.
// Naming convention for tests is (V)ariable and (C)onstant.
class MySixParameterCostFunctor {
 public:
  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const {
    const T* x0 = parameters[0];
    const T* y0 = parameters[1];
    const T* y1 = parameters[2];
    const T* z0 = parameters[3];
    const T* z1 = parameters[4];
    const T* z2 = parameters[5];

    T sum_x = x0[0];
    T sum_y = y0[0] + 2.0 * y1[0];
    T sum_z = z0[0] + 3.0 * z1[0] + 6.0 * z2[0];

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

class SixParameterCostFunctorTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    // Prepare the parameters.
    x0_ = 0.0;
    y0_ = 1.0;
    y1_ = 3.0;
    z0_ = 2.0;
    z1_ = 4.0;
    z2_ = 6.0;

    parameter_blocks_.resize(6);
    parameter_blocks_[0] = &x0_;
    parameter_blocks_[1] = &y0_;
    parameter_blocks_[2] = &y1_;
    parameter_blocks_[3] = &z0_;
    parameter_blocks_[4] = &z1_;
    parameter_blocks_[5] = &z2_;

    // Prepare the cost function.
    typedef DynamicAutoDiffCostFunction<MySixParameterCostFunctor, 3>
      DynamicMySixParameterCostFunction;
    DynamicMySixParameterCostFunction * cost_function =
      new DynamicMySixParameterCostFunction(
        new MySixParameterCostFunctor());
    for (int i = 0; i < 6; ++i) {
      cost_function->AddParameterBlock(1);
    }
    cost_function->SetNumResiduals(7);

    cost_function_.reset(cost_function);

    // Setup jacobian data.
    jacobian_vect_.resize(6);
    for (int i = 0; i < 6; ++i) {
      jacobian_vect_[i].resize(7, -100000);
    }

    // Prepare the expected residuals.
    const double sum_x = x0_;
    const double sum_y = y0_ + 2.0 * y1_;
    const double sum_z = z0_ + 3.0 * z1_ + 6.0 * z2_;

    expected_residuals_.resize(7);
    expected_residuals_[0] = sum_x;
    expected_residuals_[1] = sum_y;
    expected_residuals_[2] = sum_z;
    expected_residuals_[3] = sum_x * sum_y;
    expected_residuals_[4] = sum_y * sum_z;
    expected_residuals_[5] = sum_x * sum_z;
    expected_residuals_[6] = sum_x * sum_y * sum_z;

    // Prepare the expected jacobian entries.
    expected_jacobians_.resize(6);
    expected_jacobians_[0].resize(7);
    expected_jacobians_[0][0] = 1.0;
    expected_jacobians_[0][1] = 0.0;
    expected_jacobians_[0][2] = 0.0;
    expected_jacobians_[0][3] = sum_y;
    expected_jacobians_[0][4] = 0.0;
    expected_jacobians_[0][5] = sum_z;
    expected_jacobians_[0][6] = sum_y * sum_z;

    expected_jacobians_[1].resize(7);
    expected_jacobians_[1][0] = 0.0;
    expected_jacobians_[1][1] = 1.0;
    expected_jacobians_[1][2] = 0.0;
    expected_jacobians_[1][3] = sum_x;
    expected_jacobians_[1][4] = sum_z;
    expected_jacobians_[1][5] = 0.0;
    expected_jacobians_[1][6] = sum_x * sum_z;

    expected_jacobians_[2].resize(7);
    expected_jacobians_[2][0] = 0.0;
    expected_jacobians_[2][1] = 2.0;
    expected_jacobians_[2][2] = 0.0;
    expected_jacobians_[2][3] = 2.0 * sum_x;
    expected_jacobians_[2][4] = 2.0 * sum_z;
    expected_jacobians_[2][5] = 0.0;
    expected_jacobians_[2][6] = 2.0 * sum_x * sum_z;

    expected_jacobians_[3].resize(7);
    expected_jacobians_[3][0] = 0.0;
    expected_jacobians_[3][1] = 0.0;
    expected_jacobians_[3][2] = 1.0;
    expected_jacobians_[3][3] = 0.0;
    expected_jacobians_[3][4] = sum_y;
    expected_jacobians_[3][5] = sum_x;
    expected_jacobians_[3][6] = sum_x * sum_y;

    expected_jacobians_[4].resize(7);
    expected_jacobians_[4][0] = 0.0;
    expected_jacobians_[4][1] = 0.0;
    expected_jacobians_[4][2] = 3.0;
    expected_jacobians_[4][3] = 0.0;
    expected_jacobians_[4][4] = 3.0 * sum_y;
    expected_jacobians_[4][5] = 3.0 * sum_x;
    expected_jacobians_[4][6] = 3.0 * sum_x * sum_y;

    expected_jacobians_[5].resize(7);
    expected_jacobians_[5][0] = 0.0;
    expected_jacobians_[5][1] = 0.0;
    expected_jacobians_[5][2] = 6.0;
    expected_jacobians_[5][3] = 0.0;
    expected_jacobians_[5][4] = 6.0 * sum_y;
    expected_jacobians_[5][5] = 6.0 * sum_x;
    expected_jacobians_[5][6] = 6.0 * sum_x * sum_y;
  }

 protected:
  double x0_;
  double y0_;
  double y1_;
  double z0_;
  double z1_;
  double z2_;

  vector<double*> parameter_blocks_;

  scoped_ptr<CostFunction> cost_function_;

  vector<vector<double> > jacobian_vect_;

  vector<double> expected_residuals_;
  vector<vector<double> > expected_jacobians_;
};

TEST_F(SixParameterCostFunctorTest, TestSixParameterResiduals) {
  vector<double> residuals(7, -100000);
  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       NULL));
  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }
}

TEST_F(SixParameterCostFunctorTest, TestSixParameterJacobian) {
  vector<double> residuals(7, -100000);

  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(jacobian_vect_[1].data());
  jacobian.push_back(jacobian_vect_[2].data());
  jacobian.push_back(jacobian_vect_[3].data());
  jacobian.push_back(jacobian_vect_[4].data());
  jacobian.push_back(jacobian_vect_[5].data());

  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 7; ++j) {
      EXPECT_EQ(expected_jacobians_[i][j], jacobian[i][j]);
    }
  }
}

TEST_F(SixParameterCostFunctorTest, TestSixParameterJacobianVVCVVC) {
  vector<double> residuals(7, -100000);

  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(jacobian_vect_[1].data());
  jacobian.push_back(NULL);
  jacobian.push_back(jacobian_vect_[3].data());
  jacobian.push_back(jacobian_vect_[4].data());
  jacobian.push_back(NULL);

  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 6; ++i) {
    // Skip the constant variables.
    if (i == 2 || i == 5) {
      continue;
    }

    for (int j = 0; j < 7; ++j) {
      EXPECT_EQ(expected_jacobians_[i][j], jacobian[i][j]);
    }
  }
}

TEST_F(SixParameterCostFunctorTest, TestSixParameterJacobianVCCVCV) {
  vector<double> residuals(7, -100000);

  vector<double*> jacobian;
  jacobian.push_back(jacobian_vect_[0].data());
  jacobian.push_back(NULL);
  jacobian.push_back(NULL);
  jacobian.push_back(jacobian_vect_[3].data());
  jacobian.push_back(NULL);
  jacobian.push_back(jacobian_vect_[5].data());

  EXPECT_TRUE(cost_function_->Evaluate(parameter_blocks_.data(),
                                       residuals.data(),
                                       jacobian.data()));

  for (int i = 0; i < 7; ++i) {
    EXPECT_EQ(expected_residuals_[i], residuals[i]);
  }

  for (int i = 0; i < 6; ++i) {
    // Skip the constant variables.
    if (i == 1 || i == 2 || i == 4) {
      continue;
    }

    for (int j = 0; j < 7; ++j) {
      EXPECT_EQ(expected_jacobians_[i][j], jacobian[i][j]);
    }
  }
}

}  // namespace internal
}  // namespace ceres
