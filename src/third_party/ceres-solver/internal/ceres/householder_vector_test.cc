// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
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
// Author: vitus@google.com (Michael Vitus)

#include "ceres/householder_vector.h"
#include "ceres/internal/eigen.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

void HouseholderTestHelper(const Vector& x) {
  const double kTolerance = 1e-14;

  // Check to ensure that H * x = ||x|| * [0 ... 0 1]'.
  Vector v(x.rows());
  double beta;
  ComputeHouseholderVector(x, &v, &beta);
  Vector result = x - beta * v * (v.transpose() * x);

  Vector expected_result(x.rows());
  expected_result.setZero();
  expected_result(x.rows() - 1) = 1;
  expected_result *= x.norm();

  for (int i = 0; i < x.rows(); ++i) {
    EXPECT_NEAR(expected_result[i], result[i], kTolerance);
  }
}

TEST(HouseholderVector, ZeroPositive) {
  Vector x(3);
  x << 0.0, 0.0, 0.25;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, ZeroNegative) {
  Vector x(3);
  x << 0.0, 0.0, -0.25;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, NearZeroPositive) {
  Vector x(3);
  x << 1e-18, 1e-18, 0.25;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, NearZeroNegative) {
  Vector x(3);
  x << 1e-18, 1e-18, -0.25;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, NonZeroNegative) {
  Vector x(3);
  x << 1.0, 0.0, -3.0;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, NonZeroPositive) {
  Vector x(3);
  x << 1.0, 1.0, 1.0;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, NonZeroPositive_Size4) {
  Vector x(4);
  x << 1.0, 1.0, 0.0, 2.0;

  HouseholderTestHelper(x);
}

TEST(HouseholderVector, LastElementZero) {
  Vector x(4);
  x << 1.0, 1.0, 0.0, 0.0;

  HouseholderTestHelper(x);
}

}  // namespace internal
}  // namespace ceres
