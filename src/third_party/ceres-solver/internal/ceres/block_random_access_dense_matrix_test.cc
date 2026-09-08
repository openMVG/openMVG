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

#include "ceres/block_random_access_dense_matrix.h"

#include <vector>

#include "ceres/internal/eigen.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(BlockRandomAccessDenseMatrix, GetCell) {
  ContextImpl context;
  constexpr int num_threads = 1;

  std::vector<Block> blocks;
  blocks.emplace_back(3, 0);
  blocks.emplace_back(4, 3);
  blocks.emplace_back(5, 7);
  constexpr int num_rows = 3 + 4 + 5;
  BlockRandomAccessDenseMatrix m(blocks, &context, num_threads);
  EXPECT_EQ(m.num_rows(), num_rows);
  EXPECT_EQ(m.num_cols(), num_rows);

  for (int i = 0; i < blocks.size(); ++i) {
    const int row_idx = blocks[i].position;
    for (int j = 0; j < blocks.size(); ++j) {
      const int col_idx = blocks[j].position;
      int row;
      int col;
      int row_stride;
      int col_stride;
      CellInfo* cell = m.GetCell(i, j, &row, &col, &row_stride, &col_stride);

      EXPECT_TRUE(cell != nullptr);
      EXPECT_EQ(row, row_idx);
      EXPECT_EQ(col, col_idx);
      EXPECT_EQ(row_stride, 3 + 4 + 5);
      EXPECT_EQ(col_stride, 3 + 4 + 5);
    }
  }
}

TEST(BlockRandomAccessDenseMatrix, WriteCell) {
  ContextImpl context;
  constexpr int num_threads = 1;

  std::vector<Block> blocks;
  blocks.emplace_back(3, 0);
  blocks.emplace_back(4, 3);
  blocks.emplace_back(5, 7);
  constexpr int num_rows = 3 + 4 + 5;

  BlockRandomAccessDenseMatrix m(blocks, &context, num_threads);

  // Fill the cell (i,j) with (i + 1) * (j + 1)
  for (int i = 0; i < blocks.size(); ++i) {
    for (int j = 0; j < blocks.size(); ++j) {
      int row;
      int col;
      int row_stride;
      int col_stride;
      CellInfo* cell = m.GetCell(i, j, &row, &col, &row_stride, &col_stride);
      MatrixRef(cell->values, row_stride, col_stride)
          .block(row, col, blocks[i].size, blocks[j].size) =
          (i + 1) * (j + 1) * Matrix::Ones(blocks[i].size, blocks[j].size);
    }
  }

  // Check the values in the array are correct by going over the
  // entries of each block manually.
  for (int i = 0; i < blocks.size(); ++i) {
    const int row_idx = blocks[i].position;
    for (int j = 0; j < blocks.size(); ++j) {
      const int col_idx = blocks[j].position;
      // Check the values of this block.
      for (int r = 0; r < blocks[i].size; ++r) {
        for (int c = 0; c < blocks[j].size; ++c) {
          int pos = row_idx * num_rows + col_idx;
          EXPECT_EQ(m.values()[pos], (i + 1) * (j + 1));
        }
      }
    }
  }
}

}  // namespace ceres::internal
