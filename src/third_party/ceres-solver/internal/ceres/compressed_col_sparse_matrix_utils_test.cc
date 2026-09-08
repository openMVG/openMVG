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

#include "ceres/compressed_col_sparse_matrix_utils.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "Eigen/SparseCore"
#include "ceres/internal/export.h"
#include "ceres/triplet_sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(_, BlockPermutationToScalarPermutation) {
  //  Block structure
  //  0  --1-  ---2---  ---3---  4
  // [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  std::vector<Block> blocks{{1, 0}, {2, 1}, {3, 3}, {3, 6}, {1, 9}};
  // Block ordering
  // [1, 0, 2, 4, 5]
  std::vector<int> block_ordering{{1, 0, 2, 4, 3}};

  // Expected ordering
  // [1, 2, 0, 3, 4, 5, 9, 6, 7, 8]
  std::vector<int> expected_scalar_ordering{{1, 2, 0, 3, 4, 5, 9, 6, 7, 8}};

  std::vector<int> scalar_ordering;
  BlockOrderingToScalarOrdering(blocks, block_ordering, &scalar_ordering);
  EXPECT_EQ(scalar_ordering.size(), expected_scalar_ordering.size());
  for (int i = 0; i < expected_scalar_ordering.size(); ++i) {
    EXPECT_EQ(scalar_ordering[i], expected_scalar_ordering[i]);
  }
}

static void FillBlock(const std::vector<Block>& row_blocks,
                      const std::vector<Block>& col_blocks,
                      const int row_block_id,
                      const int col_block_id,
                      std::vector<Eigen::Triplet<double>>* triplets) {
  for (int r = 0; r < row_blocks[row_block_id].size; ++r) {
    for (int c = 0; c < col_blocks[col_block_id].size; ++c) {
      triplets->push_back(
          Eigen::Triplet<double>(row_blocks[row_block_id].position + r,
                                 col_blocks[col_block_id].position + c,
                                 1.0));
    }
  }
}

TEST(_, ScalarMatrixToBlockMatrix) {
  // Block sparsity.
  //
  //     [1 2 3 2]
  // [1]  x   x
  // [2]    x   x
  // [2]  x x
  // num_nonzeros = 1 + 3 + 4 + 4 + 1 + 2 = 15

  std::vector<Block> col_blocks{{1, 0}, {2, 1}, {3, 3}, {2, 5}};
  const int num_cols = NumScalarEntries(col_blocks);

  std::vector<Block> row_blocks{{1, 0}, {2, 1}, {2, 3}};
  const int num_rows = NumScalarEntries(row_blocks);

  std::vector<Eigen::Triplet<double>> triplets;
  FillBlock(row_blocks, col_blocks, 0, 0, &triplets);
  FillBlock(row_blocks, col_blocks, 2, 0, &triplets);
  FillBlock(row_blocks, col_blocks, 1, 1, &triplets);
  FillBlock(row_blocks, col_blocks, 2, 1, &triplets);
  FillBlock(row_blocks, col_blocks, 0, 2, &triplets);
  FillBlock(row_blocks, col_blocks, 1, 3, &triplets);
  Eigen::SparseMatrix<double> sparse_matrix(num_rows, num_cols);
  sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());

  const std::vector<int> expected_compressed_block_rows{{0, 2, 1, 2, 0, 1}};
  const std::vector<int> expected_compressed_block_cols{{0, 2, 4, 5, 6}};

  std::vector<int> compressed_block_rows;
  std::vector<int> compressed_block_cols;
  CompressedColumnScalarMatrixToBlockMatrix(sparse_matrix.innerIndexPtr(),
                                            sparse_matrix.outerIndexPtr(),
                                            row_blocks,
                                            col_blocks,
                                            &compressed_block_rows,
                                            &compressed_block_cols);

  EXPECT_EQ(compressed_block_rows, expected_compressed_block_rows);
  EXPECT_EQ(compressed_block_cols, expected_compressed_block_cols);
}

class SolveUpperTriangularTest : public ::testing::Test {
 protected:
  const std::vector<int>& cols() const { return cols_; }
  const std::vector<int>& rows() const { return rows_; }
  const std::vector<double>& values() const { return values_; }

 private:
  const std::vector<int> cols_ = {0, 1, 2, 4, 7};
  const std::vector<int> rows_ = {0, 1, 1, 2, 0, 1, 3};
  const std::vector<double> values_ = {
      0.50754, 0.80483, 0.14120, 0.3, 0.77696, 0.41860, 0.88979};
};

TEST_F(SolveUpperTriangularTest, SolveInPlace) {
  double rhs_and_solution[] = {1.0, 1.0, 2.0, 2.0};
  const double expected[] = {-1.4706, -1.0962, 6.6667, 2.2477};

  SolveUpperTriangularInPlace<int>(cols().size() - 1,
                                   rows().data(),
                                   cols().data(),
                                   values().data(),
                                   rhs_and_solution);

  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(rhs_and_solution[i], expected[i], 1e-4) << i;
  }
}

TEST_F(SolveUpperTriangularTest, TransposeSolveInPlace) {
  double rhs_and_solution[] = {1.0, 1.0, 2.0, 2.0};
  double expected[] = {1.970288, 1.242498, 6.081864, -0.057255};

  SolveUpperTriangularTransposeInPlace<int>(cols().size() - 1,
                                            rows().data(),
                                            cols().data(),
                                            values().data(),
                                            rhs_and_solution);

  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(rhs_and_solution[i], expected[i], 1e-4) << i;
  }
}

TEST_F(SolveUpperTriangularTest, RTRSolveWithSparseRHS) {
  double solution[4];
  // clang-format off
  double expected[] = { 6.8420e+00,   1.0057e+00,  -1.4907e-16,  -1.9335e+00,
                        1.0057e+00,   2.2275e+00,  -1.9493e+00,  -6.5693e-01,
                        -1.4907e-16,  -1.9493e+00,   1.1111e+01,   9.7381e-17,
                        -1.9335e+00,  -6.5693e-01,   9.7381e-17,   1.2631e+00 };
  // clang-format on

  for (int i = 0; i < 4; ++i) {
    SolveRTRWithSparseRHS<int>(cols().size() - 1,
                               rows().data(),
                               cols().data(),
                               values().data(),
                               i,
                               solution);
    for (int j = 0; j < 4; ++j) {
      EXPECT_NEAR(solution[j], expected[4 * i + j], 1e-3) << i;
    }
  }
}

}  // namespace ceres::internal
