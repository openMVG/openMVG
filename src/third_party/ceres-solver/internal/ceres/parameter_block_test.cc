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

#include "ceres/parameter_block.h"

#include "ceres/internal/eigen.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(ParameterBlock, SetManifoldDiesOnSizeMismatch) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset_wrong_size(4, indices);
  EXPECT_DEATH_IF_SUPPORTED(parameter_block.SetManifold(&subset_wrong_size),
                            "ambient");
}

TEST(ParameterBlock, SetManifoldWithSameExistingManifold) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset(3, indices);
  parameter_block.SetManifold(&subset);
  parameter_block.SetManifold(&subset);
}

TEST(ParameterBlock, SetManifoldAllowsResettingToNull) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset(3, indices);
  parameter_block.SetManifold(&subset);
  EXPECT_EQ(parameter_block.manifold(), &subset);
  parameter_block.SetManifold(nullptr);
  EXPECT_EQ(parameter_block.manifold(), nullptr);
  EXPECT_EQ(parameter_block.PlusJacobian(), nullptr);
}

TEST(ParameterBlock, SetManifoldAllowsResettingToDifferentManifold) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset(3, indices);
  parameter_block.SetManifold(&subset);
  EXPECT_EQ(parameter_block.manifold(), &subset);

  SubsetManifold subset_different(3, indices);
  parameter_block.SetManifold(&subset_different);
  EXPECT_EQ(parameter_block.manifold(), &subset_different);
}

TEST(ParameterBlock, SetManifoldAndNormalOperation) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset(3, indices);
  parameter_block.SetManifold(&subset);

  // Ensure the manifold plus jacobian result is correctly computed.
  ConstMatrixRef manifold_jacobian(parameter_block.PlusJacobian(), 3, 2);
  ASSERT_EQ(1.0, manifold_jacobian(0, 0));
  ASSERT_EQ(0.0, manifold_jacobian(0, 1));
  ASSERT_EQ(0.0, manifold_jacobian(1, 0));
  ASSERT_EQ(0.0, manifold_jacobian(1, 1));
  ASSERT_EQ(0.0, manifold_jacobian(2, 0));
  ASSERT_EQ(1.0, manifold_jacobian(2, 1));

  // Check that updating works as expected.
  double x_plus_delta[3];
  double delta[2] = {0.5, 0.3};
  parameter_block.Plus(x, delta, x_plus_delta);
  ASSERT_EQ(1.5, x_plus_delta[0]);
  ASSERT_EQ(2.0, x_plus_delta[1]);
  ASSERT_EQ(3.3, x_plus_delta[2]);
}

struct TestManifold : public Manifold {
 public:
  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const final {
    LOG(FATAL) << "Shouldn't get called.";
    return true;
  }

  bool PlusJacobian(const double* x, double* jacobian) const final {
    jacobian[0] = *x * 2;
    return true;
  }

  bool Minus(const double* y, const double* x, double* y_minus_x) const final {
    LOG(FATAL) << "Shouldn't get called";
    return true;
  }

  bool MinusJacobian(const double* x, double* jacobian) const final {
    jacobian[0] = *x * 2;
    return true;
  }

  int AmbientSize() const final { return 1; }
  int TangentSize() const final { return 1; }
};

TEST(ParameterBlock, SetStateUpdatesPlusJacobian) {
  TestManifold test_manifold;
  double x[1] = {1.0};
  ParameterBlock parameter_block(x, 1, -1, &test_manifold);

  EXPECT_EQ(2.0, *parameter_block.PlusJacobian());

  x[0] = 5.5;
  parameter_block.SetState(x);
  EXPECT_EQ(11.0, *parameter_block.PlusJacobian());
}

TEST(ParameterBlock, PlusWithNoManifold) {
  double x[2] = {1.0, 2.0};
  ParameterBlock parameter_block(x, 2, -1);

  double delta[2] = {0.2, 0.3};
  double x_plus_delta[2];
  parameter_block.Plus(x, delta, x_plus_delta);
  EXPECT_EQ(1.2, x_plus_delta[0]);
  EXPECT_EQ(2.3, x_plus_delta[1]);
}

// Stops computing the plus_jacobian after the first time.
class BadManifold : public Manifold {
 public:
  BadManifold() = default;

  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const final {
    *x_plus_delta = *x + *delta;
    return true;
  }

  bool PlusJacobian(const double* x, double* jacobian) const final {
    if (calls_ == 0) {
      jacobian[0] = 0;
    }
    ++calls_;
    return true;
  }

  bool Minus(const double* y, const double* x, double* y_minus_x) const final {
    LOG(FATAL) << "Shouldn't get called";
    return true;
  }

  bool MinusJacobian(const double* x, double* jacobian) const final {
    jacobian[0] = *x * 2;
    return true;
  }

  int AmbientSize() const final { return 1; }
  int TangentSize() const final { return 1; }

 private:
  mutable int calls_{0};
};

TEST(ParameterBlock, DetectBadManifold) {
  double x = 1;
  BadManifold bad_manifold;
  ParameterBlock parameter_block(&x, 1, -1, &bad_manifold);
  double y = 2;
  EXPECT_FALSE(parameter_block.SetState(&y));
}

TEST(ParameterBlock, DefaultBounds) {
  double x[2];
  ParameterBlock parameter_block(x, 2, -1, nullptr);
  EXPECT_EQ(parameter_block.UpperBoundForParameter(0),
            std::numeric_limits<double>::max());
  EXPECT_EQ(parameter_block.UpperBoundForParameter(1),
            std::numeric_limits<double>::max());
  EXPECT_EQ(parameter_block.LowerBoundForParameter(0),
            -std::numeric_limits<double>::max());
  EXPECT_EQ(parameter_block.LowerBoundForParameter(1),
            -std::numeric_limits<double>::max());
}

TEST(ParameterBlock, SetBounds) {
  double x[2];
  ParameterBlock parameter_block(x, 2, -1, nullptr);
  parameter_block.SetLowerBound(0, 1);
  parameter_block.SetUpperBound(1, 1);

  EXPECT_EQ(parameter_block.LowerBoundForParameter(0), 1.0);
  EXPECT_EQ(parameter_block.LowerBoundForParameter(1),
            -std::numeric_limits<double>::max());

  EXPECT_EQ(parameter_block.UpperBoundForParameter(0),
            std::numeric_limits<double>::max());
  EXPECT_EQ(parameter_block.UpperBoundForParameter(1), 1.0);
}

TEST(ParameterBlock, PlusWithBoundsConstraints) {
  double x[] = {1.0, 0.0};
  double delta[] = {2.0, -10.0};
  ParameterBlock parameter_block(x, 2, -1, nullptr);
  parameter_block.SetUpperBound(0, 2.0);
  parameter_block.SetLowerBound(1, -1.0);
  double x_plus_delta[2];
  parameter_block.Plus(x, delta, x_plus_delta);
  EXPECT_EQ(x_plus_delta[0], 2.0);
  EXPECT_EQ(x_plus_delta[1], -1.0);
}

TEST(ParameterBlock, ResetManifoldToNull) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset(3, indices);
  parameter_block.SetManifold(&subset);
  EXPECT_EQ(parameter_block.manifold(), &subset);
  parameter_block.SetManifold(nullptr);
  EXPECT_EQ(parameter_block.manifold(), nullptr);
}

TEST(ParameterBlock, ResetManifoldToNotNull) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  std::vector<int> indices;
  indices.push_back(1);
  SubsetManifold subset(3, indices);
  parameter_block.SetManifold(&subset);
  EXPECT_EQ(parameter_block.manifold(), &subset);

  SubsetManifold subset_different(3, indices);
  parameter_block.SetManifold(&subset_different);
  EXPECT_EQ(parameter_block.manifold(), &subset_different);
}

TEST(ParameterBlock, SetNullManifold) {
  double x[3] = {1.0, 2.0, 3.0};
  ParameterBlock parameter_block(x, 3, -1);
  EXPECT_EQ(parameter_block.manifold(), nullptr);

  parameter_block.SetManifold(nullptr);
  EXPECT_EQ(parameter_block.manifold(), nullptr);
}

}  // namespace internal
}  // namespace ceres
