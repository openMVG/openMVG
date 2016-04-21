// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
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

#include "ceres/array_utils.h"

#include <limits>
#include <cmath>
#include <vector>
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

using std::vector;

TEST(ArrayUtils, IsArrayValid) {
  double x[3];
  x[0] = 0.0;
  x[1] = 1.0;
  x[2] = 2.0;
  EXPECT_TRUE(IsArrayValid(3, x));
  x[1] = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(IsArrayValid(3, x));
  x[1] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(IsArrayValid(3, x));
  x[1] = std::numeric_limits<double>::signaling_NaN();
  EXPECT_FALSE(IsArrayValid(3, x));
  EXPECT_TRUE(IsArrayValid(1, NULL));
  InvalidateArray(3, x);
  EXPECT_FALSE(IsArrayValid(3, x));
}

TEST(ArrayUtils, FindInvalidIndex) {
  double x[3];
  x[0] = 0.0;
  x[1] = 1.0;
  x[2] = 2.0;
  EXPECT_EQ(FindInvalidValue(3, x), 3);
  x[1] = std::numeric_limits<double>::infinity();
  EXPECT_EQ(FindInvalidValue(3, x), 1);
  x[1] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_EQ(FindInvalidValue(3, x), 1);
  x[1] = std::numeric_limits<double>::signaling_NaN();
  EXPECT_EQ(FindInvalidValue(3, x), 1);
  EXPECT_EQ(FindInvalidValue(1, NULL), 1);
  InvalidateArray(3, x);
  EXPECT_EQ(FindInvalidValue(3, x), 0);
}

TEST(MapValuesToContiguousRange, ContiguousEntries) {
  vector<int> array;
  array.push_back(0);
  array.push_back(1);
  vector<int> expected = array;
  MapValuesToContiguousRange(array.size(), &array[0]);
  EXPECT_EQ(array, expected);
  array.clear();

  array.push_back(1);
  array.push_back(0);
  expected = array;
  MapValuesToContiguousRange(array.size(), &array[0]);
  EXPECT_EQ(array, expected);
}

TEST(MapValuesToContiguousRange, NonContiguousEntries) {
  vector<int> array;
  array.push_back(0);
  array.push_back(2);
  vector<int> expected;
  expected.push_back(0);
  expected.push_back(1);
  MapValuesToContiguousRange(array.size(), &array[0]);
  EXPECT_EQ(array, expected);
}

TEST(MapValuesToContiguousRange, NonContiguousRepeatingEntries) {
  vector<int> array;
  array.push_back(3);
  array.push_back(1);
  array.push_back(0);
  array.push_back(0);
  array.push_back(0);
  array.push_back(5);
  vector<int> expected;
  expected.push_back(2);
  expected.push_back(1);
  expected.push_back(0);
  expected.push_back(0);
  expected.push_back(0);
  expected.push_back(3);
  MapValuesToContiguousRange(array.size(), &array[0]);
  EXPECT_EQ(array, expected);
}

}  // namespace internal
}  // namespace ceres
