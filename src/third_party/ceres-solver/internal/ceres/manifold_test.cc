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

#include "ceres/manifold.h"

#include <cmath>
#include <limits>
#include <memory>
#include <utility>

#include "Eigen/Geometry"
#include "ceres/constants.h"
#include "ceres/dynamic_numeric_diff_cost_function.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/port.h"
#include "ceres/line_manifold.h"
#include "ceres/manifold_test_utils.h"
#include "ceres/numeric_diff_options.h"
#include "ceres/product_manifold.h"
#include "ceres/rotation.h"
#include "ceres/sphere_manifold.h"
#include "ceres/types.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

constexpr int kNumTrials = 1000;
constexpr double kTolerance = 1e-9;

TEST(EuclideanManifold, StaticNormalFunctionTest) {
  EuclideanManifold<3> manifold;
  EXPECT_EQ(manifold.AmbientSize(), 3);
  EXPECT_EQ(manifold.TangentSize(), 3);

  Vector zero_tangent = Vector::Zero(manifold.TangentSize());
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(manifold.AmbientSize());
    const Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());
    Vector x_plus_delta = Vector::Zero(manifold.AmbientSize());

    manifold.Plus(x.data(), delta.data(), x_plus_delta.data());
    EXPECT_NEAR((x_plus_delta - x - delta).norm() / (x + delta).norm(),
                0.0,
                kTolerance);

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(EuclideanManifold, DynamicNormalFunctionTest) {
  EuclideanManifold<DYNAMIC> manifold(3);
  EXPECT_EQ(manifold.AmbientSize(), 3);
  EXPECT_EQ(manifold.TangentSize(), 3);

  Vector zero_tangent = Vector::Zero(manifold.TangentSize());
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(manifold.AmbientSize());
    const Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());
    Vector x_plus_delta = Vector::Zero(manifold.AmbientSize());

    manifold.Plus(x.data(), delta.data(), x_plus_delta.data());
    EXPECT_NEAR((x_plus_delta - x - delta).norm() / (x + delta).norm(),
                0.0,
                kTolerance);

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(SubsetManifold, EmptyConstantParameters) {
  SubsetManifold manifold(3, {});
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(3);
    const Vector y = Vector::Random(3);
    Vector delta = Vector::Random(3);
    Vector x_plus_delta = Vector::Zero(3);

    manifold.Plus(x.data(), delta.data(), x_plus_delta.data());
    EXPECT_NEAR((x_plus_delta - x - delta).norm() / (x + delta).norm(),
                0.0,
                kTolerance);
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(SubsetManifold, NegativeParameterIndexDeathTest) {
  EXPECT_DEATH_IF_SUPPORTED(SubsetManifold manifold(2, {-1}),
                            "greater than equal to zero");
}

TEST(SubsetManifold, GreaterThanSizeParameterIndexDeathTest) {
  EXPECT_DEATH_IF_SUPPORTED(SubsetManifold manifold(2, {2}),
                            "less than the size");
}

TEST(SubsetManifold, DuplicateParametersDeathTest) {
  EXPECT_DEATH_IF_SUPPORTED(SubsetManifold manifold(2, {1, 1}), "duplicates");
}

TEST(SubsetManifold, NormalFunctionTest) {
  const int kAmbientSize = 4;
  const int kTangentSize = 3;

  for (int i = 0; i < kAmbientSize; ++i) {
    SubsetManifold manifold_with_ith_parameter_constant(kAmbientSize, {i});
    for (int trial = 0; trial < kNumTrials; ++trial) {
      const Vector x = Vector::Random(kAmbientSize);
      Vector y = Vector::Random(kAmbientSize);
      // x and y must have the same i^th coordinate to be on the manifold.
      y[i] = x[i];
      Vector delta = Vector::Random(kTangentSize);
      Vector x_plus_delta = Vector::Zero(kAmbientSize);

      x_plus_delta.setZero();
      manifold_with_ith_parameter_constant.Plus(
          x.data(), delta.data(), x_plus_delta.data());
      int k = 0;
      for (int j = 0; j < kAmbientSize; ++j) {
        if (j == i) {
          EXPECT_EQ(x_plus_delta[j], x[j]);
        } else {
          EXPECT_EQ(x_plus_delta[j], x[j] + delta[k++]);
        }
      }

      EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(
          manifold_with_ith_parameter_constant, x, delta, y, kTolerance);
    }
  }
}

TEST(ProductManifold, Size2) {
  SubsetManifold manifold1(5, {2});
  SubsetManifold manifold2(3, {0, 1});
  ProductManifold<SubsetManifold, SubsetManifold> manifold(manifold1,
                                                           manifold2);

  EXPECT_EQ(manifold.AmbientSize(),
            manifold1.AmbientSize() + manifold2.AmbientSize());
  EXPECT_EQ(manifold.TangentSize(),
            manifold1.TangentSize() + manifold2.TangentSize());
}

TEST(ProductManifold, Size3) {
  SubsetManifold manifold1(5, {2});
  SubsetManifold manifold2(3, {0, 1});
  SubsetManifold manifold3(4, {1});

  ProductManifold<SubsetManifold, SubsetManifold, SubsetManifold> manifold(
      manifold1, manifold2, manifold3);

  EXPECT_EQ(manifold.AmbientSize(),
            manifold1.AmbientSize() + manifold2.AmbientSize() +
                manifold3.AmbientSize());
  EXPECT_EQ(manifold.TangentSize(),
            manifold1.TangentSize() + manifold2.TangentSize() +
                manifold3.TangentSize());
}

TEST(ProductManifold, Size4) {
  SubsetManifold manifold1(5, {2});
  SubsetManifold manifold2(3, {0, 1});
  SubsetManifold manifold3(4, {1});
  SubsetManifold manifold4(2, {0});

  ProductManifold<SubsetManifold,
                  SubsetManifold,
                  SubsetManifold,
                  SubsetManifold>
      manifold(manifold1, manifold2, manifold3, manifold4);

  EXPECT_EQ(manifold.AmbientSize(),
            manifold1.AmbientSize() + manifold2.AmbientSize() +
                manifold3.AmbientSize() + manifold4.AmbientSize());
  EXPECT_EQ(manifold.TangentSize(),
            manifold1.TangentSize() + manifold2.TangentSize() +
                manifold3.TangentSize() + manifold4.TangentSize());
}

TEST(ProductManifold, NormalFunctionTest) {
  SubsetManifold manifold1(5, {2});
  SubsetManifold manifold2(3, {0, 1});
  SubsetManifold manifold3(4, {1});
  SubsetManifold manifold4(2, {0});

  ProductManifold<SubsetManifold,
                  SubsetManifold,
                  SubsetManifold,
                  SubsetManifold>
      manifold(manifold1, manifold2, manifold3, manifold4);

  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());
    Vector x_plus_delta = Vector::Zero(manifold.AmbientSize());
    Vector x_plus_delta_expected = Vector::Zero(manifold.AmbientSize());

    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), x_plus_delta.data()));

    int ambient_cursor = 0;
    int tangent_cursor = 0;

    EXPECT_TRUE(manifold1.Plus(&x[ambient_cursor],
                               &delta[tangent_cursor],
                               &x_plus_delta_expected[ambient_cursor]));
    ambient_cursor += manifold1.AmbientSize();
    tangent_cursor += manifold1.TangentSize();

    EXPECT_TRUE(manifold2.Plus(&x[ambient_cursor],
                               &delta[tangent_cursor],
                               &x_plus_delta_expected[ambient_cursor]));
    ambient_cursor += manifold2.AmbientSize();
    tangent_cursor += manifold2.TangentSize();

    EXPECT_TRUE(manifold3.Plus(&x[ambient_cursor],
                               &delta[tangent_cursor],
                               &x_plus_delta_expected[ambient_cursor]));
    ambient_cursor += manifold3.AmbientSize();
    tangent_cursor += manifold3.TangentSize();

    EXPECT_TRUE(manifold4.Plus(&x[ambient_cursor],
                               &delta[tangent_cursor],
                               &x_plus_delta_expected[ambient_cursor]));
    ambient_cursor += manifold4.AmbientSize();
    tangent_cursor += manifold4.TangentSize();

    for (int i = 0; i < x.size(); ++i) {
      EXPECT_EQ(x_plus_delta[i], x_plus_delta_expected[i]);
    }

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(
        manifold, x, delta, x_plus_delta, kTolerance);
  }
}

TEST(ProductManifold, ZeroTangentSizeAndEuclidean) {
  SubsetManifold subset_manifold(1, {0});
  EuclideanManifold<2> euclidean_manifold;
  ProductManifold<SubsetManifold, EuclideanManifold<2>> manifold(
      subset_manifold, euclidean_manifold);
  EXPECT_EQ(manifold.AmbientSize(), 3);
  EXPECT_EQ(manifold.TangentSize(), 2);

  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(3);
    Vector y = Vector::Random(3);
    y[0] = x[0];
    Vector delta = Vector::Random(2);
    Vector x_plus_delta = Vector::Zero(3);

    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), x_plus_delta.data()));

    EXPECT_EQ(x_plus_delta[0], x[0]);
    EXPECT_EQ(x_plus_delta[1], x[1] + delta[0]);
    EXPECT_EQ(x_plus_delta[2], x[2] + delta[1]);

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(ProductManifold, EuclideanAndZeroTangentSize) {
  SubsetManifold subset_manifold(1, {0});
  EuclideanManifold<2> euclidean_manifold;
  ProductManifold<EuclideanManifold<2>, SubsetManifold> manifold(
      euclidean_manifold, subset_manifold);
  EXPECT_EQ(manifold.AmbientSize(), 3);
  EXPECT_EQ(manifold.TangentSize(), 2);

  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(3);
    Vector y = Vector::Random(3);
    y[2] = x[2];
    Vector delta = Vector::Random(2);
    Vector x_plus_delta = Vector::Zero(3);

    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), x_plus_delta.data()));
    EXPECT_EQ(x_plus_delta[0], x[0] + delta[0]);
    EXPECT_EQ(x_plus_delta[1], x[1] + delta[1]);
    EXPECT_EQ(x_plus_delta[2], x[2]);
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

struct CopyableManifold : ceres::Manifold {
  CopyableManifold() = default;
  CopyableManifold(const CopyableManifold&) = default;
  // Do not care about copy-assignment
  CopyableManifold& operator=(const CopyableManifold&) = delete;
  // Not moveable
  CopyableManifold(CopyableManifold&&) = delete;
  CopyableManifold& operator=(CopyableManifold&&) = delete;

  int AmbientSize() const override { return 3; }
  int TangentSize() const override { return 2; }

  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override {
    return true;
  }

  bool PlusJacobian(const double* x, double* jacobian) const override {
    return true;
  }

  bool RightMultiplyByPlusJacobian(const double* x,
                                   const int num_rows,
                                   const double* ambient_matrix,
                                   double* tangent_matrix) const override {
    return true;
  }

  bool Minus(const double* y,
             const double* x,
             double* y_minus_x) const override {
    return true;
  }

  bool MinusJacobian(const double* x, double* jacobian) const override {
    return true;
  }
};

struct MoveableManifold : ceres::Manifold {
  MoveableManifold() = default;
  MoveableManifold(MoveableManifold&&) = default;
  // Do not care about move-assignment
  MoveableManifold& operator=(MoveableManifold&&) = delete;
  // Not copyable
  MoveableManifold(const MoveableManifold&) = delete;
  MoveableManifold& operator=(const MoveableManifold&) = delete;

  int AmbientSize() const override { return 3; }
  int TangentSize() const override { return 2; }

  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override {
    return true;
  }

  bool PlusJacobian(const double* x, double* jacobian) const override {
    return true;
  }

  bool RightMultiplyByPlusJacobian(const double* x,
                                   const int num_rows,
                                   const double* ambient_matrix,
                                   double* tangent_matrix) const override {
    return true;
  }

  bool Minus(const double* y,
             const double* x,
             double* y_minus_x) const override {
    return true;
  }

  bool MinusJacobian(const double* x, double* jacobian) const override {
    return true;
  }
};

TEST(ProductManifold, CopyableOnly) {
  ProductManifold<CopyableManifold, EuclideanManifold<3>> manifold1{
      CopyableManifold{}, EuclideanManifold<3>{}};

  CopyableManifold inner2;
  ProductManifold<CopyableManifold, EuclideanManifold<3>> manifold2{
      inner2, EuclideanManifold<3>{}};

  EXPECT_EQ(manifold1.AmbientSize(), manifold2.AmbientSize());
  EXPECT_EQ(manifold1.TangentSize(), manifold2.TangentSize());
}

TEST(ProductManifold, MoveableOnly) {
  ProductManifold<MoveableManifold, EuclideanManifold<3>> manifold1{
      MoveableManifold{}, EuclideanManifold<3>{}};

  MoveableManifold inner2;
  ProductManifold<MoveableManifold, EuclideanManifold<3>> manifold2{
      std::move(inner2), EuclideanManifold<3>{}};

  EXPECT_EQ(manifold1.AmbientSize(), manifold2.AmbientSize());
  EXPECT_EQ(manifold1.TangentSize(), manifold2.TangentSize());
}

TEST(ProductManifold, CopyableOrMoveable) {
  const CopyableManifold inner12{};
  ProductManifold<MoveableManifold, CopyableManifold> manifold1{
      MoveableManifold{}, inner12};

  MoveableManifold inner21;
  CopyableManifold inner22;
  ProductManifold<MoveableManifold, CopyableManifold> manifold2{
      std::move(inner21), inner22};

  EXPECT_EQ(manifold1.AmbientSize(), manifold2.AmbientSize());
  EXPECT_EQ(manifold1.TangentSize(), manifold2.TangentSize());
}

struct NonDefaultConstructibleManifold : ceres::Manifold {
  NonDefaultConstructibleManifold(int, int) {}
  int AmbientSize() const override { return 4; }
  int TangentSize() const override { return 3; }

  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override {
    return true;
  }

  bool PlusJacobian(const double* x, double* jacobian) const override {
    return true;
  }

  bool RightMultiplyByPlusJacobian(const double* x,
                                   const int num_rows,
                                   const double* ambient_matrix,
                                   double* tangent_matrix) const override {
    return true;
  }

  bool Minus(const double* y,
             const double* x,
             double* y_minus_x) const override {
    return true;
  }

  bool MinusJacobian(const double* x, double* jacobian) const override {
    return true;
  }
};

TEST(ProductManifold, NonDefaultConstructible) {
  ProductManifold<NonDefaultConstructibleManifold, QuaternionManifold>
      manifold1{NonDefaultConstructibleManifold{1, 2}, QuaternionManifold{}};
  ProductManifold<QuaternionManifold, NonDefaultConstructibleManifold>
      manifold2{QuaternionManifold{}, NonDefaultConstructibleManifold{1, 2}};

  EXPECT_EQ(manifold1.AmbientSize(), manifold2.AmbientSize());
  EXPECT_EQ(manifold1.TangentSize(), manifold2.TangentSize());
}

TEST(ProductManifold, DefaultConstructible) {
  ProductManifold<EuclideanManifold<3>, SphereManifold<4>> manifold1;
  ProductManifold<SphereManifold<4>, EuclideanManifold<3>> manifold2;

  EXPECT_EQ(manifold1.AmbientSize(), manifold2.AmbientSize());
  EXPECT_EQ(manifold1.TangentSize(), manifold2.TangentSize());
}

TEST(ProductManifold, Pointers) {
  auto p = std::make_unique<QuaternionManifold>();
  auto q = std::make_shared<EuclideanManifold<3>>();

  ProductManifold<std::unique_ptr<Manifold>,
                  EuclideanManifold<3>,
                  std::shared_ptr<EuclideanManifold<3>>>
      manifold1{
          std::make_unique<QuaternionManifold>(), EuclideanManifold<3>{}, q};
  ProductManifold<QuaternionManifold*,
                  EuclideanManifold<3>,
                  std::shared_ptr<EuclideanManifold<3>>>
      manifold2{p.get(), EuclideanManifold<3>{}, q};

  EXPECT_EQ(manifold1.AmbientSize(), manifold2.AmbientSize());
  EXPECT_EQ(manifold1.TangentSize(), manifold2.TangentSize());
}

TEST(QuaternionManifold, PlusPiBy2) {
  QuaternionManifold manifold;
  Vector x = Vector::Zero(4);
  x[0] = 1.0;

  for (int i = 0; i < 3; ++i) {
    Vector delta = Vector::Zero(3);
    delta[i] = constants::pi / 2;
    Vector x_plus_delta = Vector::Zero(4);
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), x_plus_delta.data()));

    // Expect that the element corresponding to pi/2 is +/- 1. All other
    // elements should be zero.
    for (int j = 0; j < 4; ++j) {
      if (i == (j - 1)) {
        EXPECT_LT(std::abs(x_plus_delta[j]) - 1,
                  std::numeric_limits<double>::epsilon())
            << "\ndelta = " << delta.transpose()
            << "\nx_plus_delta = " << x_plus_delta.transpose()
            << "\n expected the " << j
            << "th element of x_plus_delta to be +/- 1.";
      } else {
        EXPECT_LT(std::abs(x_plus_delta[j]),
                  std::numeric_limits<double>::epsilon())
            << "\ndelta = " << delta.transpose()
            << "\nx_plus_delta = " << x_plus_delta.transpose()
            << "\n expected the " << j << "th element of x_plus_delta to be 0.";
      }
    }
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(
        manifold, x, delta, x_plus_delta, kTolerance);
  }
}

// Compute the expected value of QuaternionManifold::Plus via functions in
// rotation.h and compares it to the one computed by QuaternionManifold::Plus.
MATCHER_P2(QuaternionManifoldPlusIsCorrectAt, x, delta, "") {
  // This multiplication by 2 is needed because AngleAxisToQuaternion uses
  // |delta|/2 as the angle of rotation where as in the implementation of
  // QuaternionManifold for historical reasons we use |delta|.
  const Vector two_delta = delta * 2;
  Vector delta_q(4);
  AngleAxisToQuaternion(two_delta.data(), delta_q.data());

  Vector expected(4);
  QuaternionProduct(delta_q.data(), x.data(), expected.data());
  Vector actual(4);
  EXPECT_TRUE(arg.Plus(x.data(), delta.data(), actual.data()));

  const double n = (actual - expected).norm();
  const double d = expected.norm();
  const double diffnorm = n / d;
  if (diffnorm > kTolerance) {
    *result_listener << "\nx: " << x.transpose()
                     << "\ndelta: " << delta.transpose()
                     << "\nexpected: " << expected.transpose()
                     << "\nactual: " << actual.transpose()
                     << "\ndiff: " << (expected - actual).transpose()
                     << "\ndiffnorm : " << diffnorm;
    return false;
  }
  return true;
}

static Vector RandomQuaternion() {
  Vector x = Vector::Random(4);
  x.normalize();
  return x;
}

TEST(QuaternionManifold, GenericDelta) {
  QuaternionManifold manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    EXPECT_THAT(manifold, QuaternionManifoldPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(QuaternionManifold, SmallDelta) {
  QuaternionManifold manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    delta.normalize();
    delta *= 1e-6;
    EXPECT_THAT(manifold, QuaternionManifoldPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(QuaternionManifold, DeltaJustBelowPi) {
  QuaternionManifold manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    delta.normalize();
    delta *= (constants::pi - 1e-6);
    EXPECT_THAT(manifold, QuaternionManifoldPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

// Compute the expected value of EigenQuaternionManifold::Plus using Eigen and
// compares it to the one computed by QuaternionManifold::Plus.
MATCHER_P2(EigenQuaternionManifoldPlusIsCorrectAt, x, delta, "") {
  // This multiplication by 2 is needed because AngleAxisToQuaternion uses
  // |delta|/2 as the angle of rotation where as in the implementation of
  // Quaternion for historical reasons we use |delta|.
  const Vector two_delta = delta * 2;
  Vector delta_q(4);
  AngleAxisToQuaternion(two_delta.data(), delta_q.data());
  Eigen::Quaterniond delta_eigen_q(
      delta_q[0], delta_q[1], delta_q[2], delta_q[3]);

  Eigen::Map<const Eigen::Quaterniond> x_eigen_q(x.data());

  Eigen::Quaterniond expected = delta_eigen_q * x_eigen_q;
  double actual[4];
  EXPECT_TRUE(arg.Plus(x.data(), delta.data(), actual));
  Eigen::Map<Eigen::Quaterniond> actual_eigen_q(actual);

  const double n = (actual_eigen_q.coeffs() - expected.coeffs()).norm();
  const double d = expected.norm();
  const double diffnorm = n / d;
  if (diffnorm > kTolerance) {
    *result_listener
        << "\nx: " << x.transpose() << "\ndelta: " << delta.transpose()
        << "\nexpected: " << expected.coeffs().transpose()
        << "\nactual: " << actual_eigen_q.coeffs().transpose() << "\ndiff: "
        << (expected.coeffs() - actual_eigen_q.coeffs()).transpose()
        << "\ndiffnorm : " << diffnorm;
    return false;
  }
  return true;
}

TEST(EigenQuaternionManifold, GenericDelta) {
  EigenQuaternionManifold manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    EXPECT_THAT(manifold, EigenQuaternionManifoldPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(EigenQuaternionManifold, SmallDelta) {
  EigenQuaternionManifold manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    delta.normalize();
    delta *= 1e-6;
    EXPECT_THAT(manifold, EigenQuaternionManifoldPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(EigenQuaternionManifold, DeltaJustBelowPi) {
  EigenQuaternionManifold manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    delta.normalize();
    delta *= (constants::pi - 1e-6);
    EXPECT_THAT(manifold, EigenQuaternionManifoldPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

using Eigen::Vector2d;
using Eigen::Vector3d;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Eigen::Vector4d;
using Vector8d = Eigen::Matrix<double, 8, 1>;

TEST(SphereManifold, ZeroTest) {
  Vector4d x{0.0, 0.0, 0.0, 1.0};
  Vector3d delta = Vector3d::Zero();
  Vector4d y = Vector4d::Zero();

  SphereManifold<4> manifold;
  manifold.Plus(x.data(), delta.data(), y.data());
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(SphereManifold, NearZeroTest1) {
  Vector4d x{1e-5, 1e-5, 1e-5, 1.0};
  x.normalize();
  Vector3d delta{0.0, 1.0, 0.0};
  Vector4d y = Vector4d::Zero();

  SphereManifold<4> manifold;
  manifold.Plus(x.data(), delta.data(), y.data());
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(SphereManifold, NearZeroTest2) {
  Vector4d x{0.01, 0.0, 0.0, 0.0};
  Vector3d delta{0.0, 1.0, 0.0};
  Vector4d y = Vector4d::Zero();
  SphereManifold<4> manifold;
  manifold.Plus(x.data(), delta.data(), y.data());
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(SphereManifold, Plus2DTest) {
  Eigen::Vector2d x{0.0, 1.0};
  SphereManifold<2> manifold;

  {
    double delta[1]{constants::pi / 4};
    Eigen::Vector2d y = Eigen::Vector2d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta, y.data()));
    const Eigen::Vector2d gtY(std::sqrt(2.0) / 2.0, std::sqrt(2.0) / 2.0);
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    double delta[1]{constants::pi / 2};
    Eigen::Vector2d y = Eigen::Vector2d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta, y.data()));
    const Eigen::Vector2d gtY = Eigen::Vector2d::UnitX();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    double delta[1]{constants::pi};
    Eigen::Vector2d y = Eigen::Vector2d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta, y.data()));
    const Eigen::Vector2d gtY = -Eigen::Vector2d::UnitY();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    double delta[1]{2.0 * constants::pi};
    Eigen::Vector2d y = Eigen::Vector2d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta, y.data()));
    const Eigen::Vector2d gtY = Eigen::Vector2d::UnitY();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }
}

TEST(SphereManifold, Plus3DTest) {
  Eigen::Vector3d x{0.0, 0.0, 1.0};
  SphereManifold<3> manifold;

  {
    Eigen::Vector2d delta{constants::pi / 2, 0.0};
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = Eigen::Vector3d::UnitX();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta{constants::pi, 0.0};
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = -Eigen::Vector3d::UnitZ();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta{2.0 * constants::pi, 0.0};
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = Eigen::Vector3d::UnitZ();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta{0.0, constants::pi / 2};
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = Eigen::Vector3d::UnitY();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta{0.0, constants::pi};
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = -Eigen::Vector3d::UnitZ();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta{0.0, 2.0 * constants::pi};
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = Eigen::Vector3d::UnitZ();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta =
        Eigen::Vector2d(1, 1).normalized() * constants::pi / 2;
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY(std::sqrt(2.0) / 2.0, std::sqrt(2.0) / 2.0, 0.0);
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta = Eigen::Vector2d(1, 1).normalized() * constants::pi;
    Eigen::Vector3d y = Eigen::Vector3d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    const Eigen::Vector3d gtY = -Eigen::Vector3d::UnitZ();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }
}

TEST(SphereManifold, Minus2DTest) {
  Eigen::Vector2d x{1.0, 0.0};
  SphereManifold<2> manifold;

  {
    double delta[1];
    const Eigen::Vector2d y(std::sqrt(2.0) / 2.0, std::sqrt(2.0) / 2.0);
    const double gtDelta{constants::pi / 4};
    EXPECT_TRUE(manifold.Minus(y.data(), x.data(), delta));
    EXPECT_LT(std::abs(delta[0] - gtDelta), kTolerance);
  }

  {
    double delta[1];
    const Eigen::Vector2d y(-1, 0);
    const double gtDelta{constants::pi};
    EXPECT_TRUE(manifold.Minus(y.data(), x.data(), delta));
    EXPECT_LT(std::abs(delta[0] - gtDelta), kTolerance);
  }
}

TEST(SphereManifold, Minus3DTest) {
  Eigen::Vector3d x{1.0, 0.0, 0.0};
  SphereManifold<3> manifold;

  {
    Eigen::Vector2d delta;
    const Eigen::Vector3d y(std::sqrt(2.0) / 2.0, 0.0, std::sqrt(2.0) / 2.0);
    const Eigen::Vector2d gtDelta(constants::pi / 4, 0.0);
    EXPECT_TRUE(manifold.Minus(y.data(), x.data(), delta.data()));
    EXPECT_LT((delta - gtDelta).norm(), kTolerance);
  }

  {
    Eigen::Vector2d delta;
    const Eigen::Vector3d y(-1, 0, 0);
    const Eigen::Vector2d gtDelta(0.0, constants::pi);
    EXPECT_TRUE(manifold.Minus(y.data(), x.data(), delta.data()));
    EXPECT_LT((delta - gtDelta).norm(), kTolerance);
  }
}

TEST(SphereManifold, DeathTests) {
  EXPECT_DEATH_IF_SUPPORTED(SphereManifold<Eigen::Dynamic> x(1), "size");
}

TEST(SphereManifold, NormalFunctionTest) {
  SphereManifold<4> manifold;
  EXPECT_EQ(manifold.AmbientSize(), 4);
  EXPECT_EQ(manifold.TangentSize(), 3);

  Vector zero_tangent = Vector::Zero(manifold.TangentSize());
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(manifold.AmbientSize());
    Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());

    if (x.norm() == 0.0 || y.norm() == 0.0) {
      continue;
    }

    // X and y need to have the same length.
    y *= x.norm() / y.norm();

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(SphereManifold, NormalFunctionTestDynamic) {
  SphereManifold<ceres::DYNAMIC> manifold(5);
  EXPECT_EQ(manifold.AmbientSize(), 5);
  EXPECT_EQ(manifold.TangentSize(), 4);

  Vector zero_tangent = Vector::Zero(manifold.TangentSize());
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(manifold.AmbientSize());
    Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());

    if (x.norm() == 0.0 || y.norm() == 0.0) {
      continue;
    }

    // X and y need to have the same length.
    y *= x.norm() / y.norm();

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(LineManifold, ZeroTest3D) {
  const Vector6d x = Vector6d::Unit(5);
  const Vector4d delta = Vector4d::Zero();
  Vector6d y = Vector6d::Zero();

  LineManifold<3> manifold;
  EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(LineManifold, ZeroTest4D) {
  const Vector8d x = Vector8d::Unit(7);
  const Vector6d delta = Vector6d::Zero();
  Vector8d y = Vector8d::Zero();

  LineManifold<4> manifold;
  EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(LineManifold, ZeroOriginPointTest3D) {
  const Vector6d x = Vector6d::Unit(5);
  Vector4d delta;
  delta << 0.0, 0.0, 1.0, 2.0;
  Vector6d y = Vector6d::Zero();

  LineManifold<3> manifold;
  EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(LineManifold, ZeroOriginPointTest4D) {
  const Vector8d x = Vector8d::Unit(7);
  Vector6d delta;
  delta << 0.0, 0.0, 0.0, 0.5, 1.0, 1.5;
  Vector8d y = Vector8d::Zero();

  LineManifold<4> manifold;
  EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(LineManifold, ZeroDirTest3D) {
  Vector6d x = Vector6d::Unit(5);
  Vector4d delta;
  delta << 3.0, 2.0, 0.0, 0.0;
  Vector6d y = Vector6d::Zero();

  LineManifold<3> manifold;
  EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(LineManifold, ZeroDirTest4D) {
  Vector8d x = Vector8d::Unit(7);
  Vector6d delta;
  delta << 3.0, 2.0, 1.0, 0.0, 0.0, 0.0;
  Vector8d y = Vector8d::Zero();

  LineManifold<4> manifold;
  EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
  EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
}

TEST(LineManifold, Plus) {
  Vector6d x = Vector6d::Unit(5);
  LineManifold<3> manifold;

  {
    Vector4d delta{0.0, 2.0, constants::pi / 2, 0.0};
    Vector6d y = Vector6d::Random();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    Vector6d gtY;
    gtY << 2.0 * Vector3d::UnitY(), Vector3d::UnitX();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Vector4d delta{3.0, 0.0, 0.0, constants::pi / 2};
    Vector6d y = Vector6d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    Vector6d gtY;
    gtY << 3.0 * Vector3d::UnitX(), Vector3d::UnitY();
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }

  {
    Vector4d delta;
    delta << Vector2d(1.0, 2.0),
        Vector2d(1, 1).normalized() * constants::pi / 2;
    Vector6d y = Vector6d::Zero();
    EXPECT_TRUE(manifold.Plus(x.data(), delta.data(), y.data()));
    Vector6d gtY;
    gtY << Vector3d(1.0, 2.0, 0.0),
        Vector3d(std::sqrt(2.0) / 2.0, std::sqrt(2.0) / 2.0, 0.0);
    EXPECT_LT((y - gtY).norm(), kTolerance);
  }
}

TEST(LineManifold, NormalFunctionTest) {
  LineManifold<3> manifold;
  EXPECT_EQ(manifold.AmbientSize(), 6);
  EXPECT_EQ(manifold.TangentSize(), 4);

  Vector zero_tangent = Vector::Zero(manifold.TangentSize());
  for (int trial = 0; trial < kNumTrials; ++trial) {
    Vector x = Vector::Random(manifold.AmbientSize());
    Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());

    if (x.tail<3>().norm() == 0.0) {
      continue;
    }

    x.tail<3>().normalize();
    manifold.Plus(x.data(), delta.data(), y.data());

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(LineManifold, NormalFunctionTestDynamic) {
  LineManifold<ceres::DYNAMIC> manifold(3);
  EXPECT_EQ(manifold.AmbientSize(), 6);
  EXPECT_EQ(manifold.TangentSize(), 4);

  Vector zero_tangent = Vector::Zero(manifold.TangentSize());
  for (int trial = 0; trial < kNumTrials; ++trial) {
    Vector x = Vector::Random(manifold.AmbientSize());
    Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());

    if (x.tail<3>().norm() == 0.0) {
      continue;
    }

    x.tail<3>().normalize();
    manifold.Plus(x.data(), delta.data(), y.data());

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

}  // namespace ceres::internal
