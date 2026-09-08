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

#include "ceres/autodiff_manifold.h"

#include <cmath>

#include "ceres/constants.h"
#include "ceres/manifold.h"
#include "ceres/manifold_test_utils.h"
#include "ceres/rotation.h"
#include "gtest/gtest.h"

namespace ceres::internal {

namespace {

constexpr int kNumTrials = 1000;
constexpr double kTolerance = 1e-9;

Vector RandomQuaternion() {
  Vector x = Vector::Random(4);
  x.normalize();
  return x;
}

}  // namespace

struct EuclideanFunctor {
  template <typename T>
  bool Plus(const T* x, const T* delta, T* x_plus_delta) const {
    for (int i = 0; i < 3; ++i) {
      x_plus_delta[i] = x[i] + delta[i];
    }
    return true;
  }

  template <typename T>
  bool Minus(const T* y, const T* x, T* y_minus_x) const {
    for (int i = 0; i < 3; ++i) {
      y_minus_x[i] = y[i] - x[i];
    }
    return true;
  }
};

TEST(AutoDiffLManifoldTest, EuclideanManifold) {
  AutoDiffManifold<EuclideanFunctor, 3, 3> manifold;
  EXPECT_EQ(manifold.AmbientSize(), 3);
  EXPECT_EQ(manifold.TangentSize(), 3);

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

struct ScaledFunctor {
  explicit ScaledFunctor(const double s) : s(s) {}

  template <typename T>
  bool Plus(const T* x, const T* delta, T* x_plus_delta) const {
    for (int i = 0; i < 3; ++i) {
      x_plus_delta[i] = x[i] + s * delta[i];
    }
    return true;
  }

  template <typename T>
  bool Minus(const T* y, const T* x, T* y_minus_x) const {
    for (int i = 0; i < 3; ++i) {
      y_minus_x[i] = (y[i] - x[i]) / s;
    }
    return true;
  }

  const double s;
};

TEST(AutoDiffManifoldTest, ScaledManifold) {
  constexpr double kScale = 1.2342;
  AutoDiffManifold<ScaledFunctor, 3, 3> manifold(new ScaledFunctor(kScale));
  EXPECT_EQ(manifold.AmbientSize(), 3);
  EXPECT_EQ(manifold.TangentSize(), 3);

  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = Vector::Random(manifold.AmbientSize());
    const Vector y = Vector::Random(manifold.AmbientSize());
    Vector delta = Vector::Random(manifold.TangentSize());
    Vector x_plus_delta = Vector::Zero(manifold.AmbientSize());

    manifold.Plus(x.data(), delta.data(), x_plus_delta.data());
    EXPECT_NEAR((x_plus_delta - x - delta * kScale).norm() /
                    (x + delta * kScale).norm(),
                0.0,
                kTolerance);

    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

// Templated functor that implements the Plus and Minus operations on the
// Quaternion manifold.
struct QuaternionFunctor {
  template <typename T>
  bool Plus(const T* x, const T* delta, T* x_plus_delta) const {
    const T squared_norm_delta =
        delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

    T q_delta[4];
    if (squared_norm_delta > T(0.0)) {
      T norm_delta = sqrt(squared_norm_delta);
      const T sin_delta_by_delta = sin(norm_delta) / norm_delta;
      q_delta[0] = cos(norm_delta);
      q_delta[1] = sin_delta_by_delta * delta[0];
      q_delta[2] = sin_delta_by_delta * delta[1];
      q_delta[3] = sin_delta_by_delta * delta[2];
    } else {
      // We do not just use q_delta = [1,0,0,0] here because that is a
      // constant and when used for automatic differentiation will
      // lead to a zero derivative. Instead we take a first order
      // approximation and evaluate it at zero.
      q_delta[0] = T(1.0);
      q_delta[1] = delta[0];
      q_delta[2] = delta[1];
      q_delta[3] = delta[2];
    }

    QuaternionProduct(q_delta, x, x_plus_delta);
    return true;
  }

  template <typename T>
  bool Minus(const T* y, const T* x, T* y_minus_x) const {
    T minus_x[4] = {x[0], -x[1], -x[2], -x[3]};
    T ambient_y_minus_x[4];
    QuaternionProduct(y, minus_x, ambient_y_minus_x);
    T u_norm = sqrt(ambient_y_minus_x[1] * ambient_y_minus_x[1] +
                    ambient_y_minus_x[2] * ambient_y_minus_x[2] +
                    ambient_y_minus_x[3] * ambient_y_minus_x[3]);
    if (u_norm > 0.0) {
      T theta = atan2(u_norm, ambient_y_minus_x[0]);
      y_minus_x[0] = theta * ambient_y_minus_x[1] / u_norm;
      y_minus_x[1] = theta * ambient_y_minus_x[2] / u_norm;
      y_minus_x[2] = theta * ambient_y_minus_x[3] / u_norm;
    } else {
      // We do not use [0,0,0] here because even though the value part is
      // a constant, the derivative part is not.
      y_minus_x[0] = ambient_y_minus_x[1];
      y_minus_x[1] = ambient_y_minus_x[2];
      y_minus_x[2] = ambient_y_minus_x[3];
    }
    return true;
  }
};

TEST(AutoDiffManifoldTest, QuaternionPlusPiBy2) {
  AutoDiffManifold<QuaternionFunctor, 4, 3> manifold;

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

// Compute the expected value of Quaternion::Plus via functions in rotation.h
// and compares it to the one computed by Quaternion::Plus.
MATCHER_P2(QuaternionPlusIsCorrectAt, x, delta, "") {
  // This multiplication by 2 is needed because AngleAxisToQuaternion uses
  // |delta|/2 as the angle of rotation where as in the implementation of
  // Quaternion for historical reasons we use |delta|.
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

TEST(AutoDiffManifoldTest, QuaternionGenericDelta) {
  AutoDiffManifold<QuaternionFunctor, 4, 3> manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    EXPECT_THAT(manifold, QuaternionPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(AutoDiffManifoldTest, QuaternionSmallDelta) {
  AutoDiffManifold<QuaternionFunctor, 4, 3> manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    delta.normalize();
    delta *= 1e-6;
    EXPECT_THAT(manifold, QuaternionPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

TEST(AutoDiffManifold, QuaternionDeltaJustBelowPi) {
  AutoDiffManifold<QuaternionFunctor, 4, 3> manifold;
  for (int trial = 0; trial < kNumTrials; ++trial) {
    const Vector x = RandomQuaternion();
    const Vector y = RandomQuaternion();
    Vector delta = Vector::Random(3);
    delta.normalize();
    delta *= (constants::pi - 1e-6);
    EXPECT_THAT(manifold, QuaternionPlusIsCorrectAt(x, delta));
    EXPECT_THAT_MANIFOLD_INVARIANTS_HOLD(manifold, x, delta, y, kTolerance);
  }
}

}  // namespace ceres::internal
