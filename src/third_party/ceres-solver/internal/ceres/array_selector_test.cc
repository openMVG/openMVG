// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
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
// Author: darius.rueckert@fau.de (Darius Rueckert)
//

#include "ceres/internal/array_selector.h"

#include "gtest/gtest.h"

namespace ceres::internal {

// This test only checks, if the correct array implementations are selected. The
// test for FixedArray is in fixed_array_test.cc. Tests for std::array and
// std::vector are not included in ceres.
TEST(ArraySelector, FixedArray) {
  ArraySelector<int, DYNAMIC, 20> array1(10);
  static_assert(
      std::is_base_of<internal::FixedArray<int, 20>, decltype(array1)>::value);
  EXPECT_EQ(array1.size(), 10);

  ArraySelector<int, DYNAMIC, 10> array2(20);
  static_assert(
      std::is_base_of<internal::FixedArray<int, 10>, decltype(array2)>::value);
  EXPECT_EQ(array2.size(), 20);
}

TEST(ArraySelector, Array) {
  ArraySelector<int, 10, 20> array1(10);
  static_assert(std::is_base_of<std::array<int, 10>, decltype(array1)>::value);
  EXPECT_EQ(array1.size(), 10);

  ArraySelector<int, 20, 20> array2(20);
  static_assert(std::is_base_of<std::array<int, 20>, decltype(array2)>::value);
  EXPECT_EQ(array2.size(), 20);
}

TEST(ArraySelector, Vector) {
  ArraySelector<int, 20, 10> array1(20);
  static_assert(std::is_base_of<std::vector<int>, decltype(array1)>::value);
  EXPECT_EQ(array1.size(), 20);

  ArraySelector<int, 1, 0> array2(1);
  static_assert(std::is_base_of<std::vector<int>, decltype(array2)>::value);
  EXPECT_EQ(array2.size(), 1);
}

}  // namespace ceres::internal
