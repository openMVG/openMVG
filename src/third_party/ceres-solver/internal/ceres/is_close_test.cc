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
// Author: dgossow@google.com (David Gossow)
//
// This file contains tests for the IsClose function.

#include "ceres/is_close.h"

#include "gtest/gtest.h"

namespace ceres::internal {

const double kTolerance = 1e-9;

TEST(IsClose, BothParametersPositive) {
  double relative_error = -1;
  double absolute_error = -1;

  // Test cases where both values are positive.
  EXPECT_TRUE(IsClose(9.9, 10.0, 0.011, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(10.0, 9.9, 0.011, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
  relative_error = -1;
  absolute_error = -1;

  EXPECT_FALSE(IsClose(9.9, 10.0, 0.009, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(10.0, 9.9, 0.009, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
}

TEST(IsClose, BothParametersNegative) {
  double relative_error = -1;
  double absolute_error = -1;

  // Test cases where both values are negative.
  EXPECT_TRUE(IsClose(-9.9, -10.0, 0.011, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(-10.0, -9.9, 0.011, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
  relative_error = -1;
  absolute_error = -1;

  EXPECT_FALSE(IsClose(-9.9, -10.0, 0.009, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(-10.0, -9.9, 0.009, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.01, kTolerance);
  EXPECT_NEAR(absolute_error, 0.1, kTolerance);
}

TEST(IsClose, ParametersHaveMixedSigns) {
  double relative_error = -1;
  double absolute_error = -1;

  // Test cases with mixed signs.
  EXPECT_FALSE(IsClose(-0.1, 0.1, 1.99, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 2.0, kTolerance);
  EXPECT_NEAR(absolute_error, 0.2, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(-0.1, 0.1, 2.01, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 2.0, kTolerance);
  EXPECT_NEAR(absolute_error, 0.2, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(0.1, -0.1, 1.99, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 2.0, kTolerance);
  EXPECT_NEAR(absolute_error, 0.2, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(0.1, -0.1, 2.01, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 2.0, kTolerance);
  EXPECT_NEAR(absolute_error, 0.2, kTolerance);
}

TEST(IsClose, OneParameterZero) {
  double relative_error = -1;
  double absolute_error = -1;

  // Test cases where one of the values is zero.
  EXPECT_TRUE(IsClose(0.0, 10.0, 10.1, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(10.0, 0.0, 10.1, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(0.0, -10.0, 10.1, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_TRUE(IsClose(-10.0, 0.0, 10.1, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;

  EXPECT_FALSE(IsClose(0, 10.0, 9.9, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(10.0, 0.0, 9.9, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(0, -10.0, 9.9, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(-10.0, 0.0, 9.9, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 10.0, kTolerance);
  EXPECT_NEAR(absolute_error, 10.0, kTolerance);
}

TEST(IsClose, BothParametersZero) {
  double relative_error = -1;
  double absolute_error = -1;
  EXPECT_TRUE(IsClose(0.0, 0.0, 0.1, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.0, kTolerance);
  EXPECT_NEAR(absolute_error, 0.0, kTolerance);
  relative_error = -1;
  absolute_error = -1;
  EXPECT_FALSE(IsClose(0.0, 0.0, 0.0, &relative_error, &absolute_error));
  EXPECT_NEAR(relative_error, 0.0, kTolerance);
  EXPECT_NEAR(absolute_error, 0.0, kTolerance);
}
}  // namespace ceres::internal
