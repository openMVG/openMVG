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
// Author: wjr@google.com (William Rucklidge)

#include "ceres/parallel_utils.h"

#include "ceres/internal/config.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// Tests that unfolding linear iterations to triangular iterations produces
// indices that are in-range and unique.
TEST(LinearIndexToUpperTriangularIndexTest, UniqueAndValid) {
  for (int n = 0; n < 100; n++) {
    std::set<std::pair<int, int>> seen_pairs;
    int actual_work_items = (n * (n + 1)) / 2;
    for (int k = 0; k < actual_work_items; k++) {
      int i, j;
      LinearIndexToUpperTriangularIndex(k, n, &i, &j);
      EXPECT_GE(i, 0);
      EXPECT_LT(i, n);
      EXPECT_GE(j, i);
      EXPECT_LT(j, n);
      seen_pairs.insert(std::make_pair(i, j));
    }
    EXPECT_EQ(actual_work_items, seen_pairs.size());
  }
}

}  // namespace ceres::internal
