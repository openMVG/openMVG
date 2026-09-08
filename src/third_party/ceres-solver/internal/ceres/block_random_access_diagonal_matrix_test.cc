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

#include "ceres/block_random_access_diagonal_matrix.h"

#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Cholesky"
#include "ceres/internal/eigen.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

class BlockRandomAccessDiagonalMatrixTest : public ::testing::Test {
 public:
  void SetUp() override {
    std::vector<Block> blocks;
    blocks.emplace_back(3, 0);
    blocks.emplace_back(4, 3);
    blocks.emplace_back(5, 7);

    const int num_rows = 3 + 4 + 5;
    num_nonzeros_ = 3 * 3 + 4 * 4 + 5 * 5;

    m_ =
        std::make_unique<BlockRandomAccessDiagonalMatrix>(blocks, &context_, 1);

    EXPECT_EQ(m_->num_rows(), num_rows);
    EXPECT_EQ(m_->num_cols(), num_rows);

    for (int i = 0; i < blocks.size(); ++i) {
      const int row_block_id = i;
      int col_block_id;
      int row;
      int col;
      int row_stride;
      int col_stride;

      for (int j = 0; j < blocks.size(); ++j) {
        col_block_id = j;
        CellInfo* cell = m_->GetCell(
            row_block_id, col_block_id, &row, &col, &row_stride, &col_stride);
        // Off diagonal entries are not present.
        if (i != j) {
          EXPECT_TRUE(cell == nullptr);
          continue;
        }

        EXPECT_TRUE(cell != nullptr);
        EXPECT_EQ(row, 0);
        EXPECT_EQ(col, 0);
        EXPECT_EQ(row_stride, blocks[row_block_id].size);
        EXPECT_EQ(col_stride, blocks[col_block_id].size);

        // Write into the block
        MatrixRef(cell->values, row_stride, col_stride)
            .block(row,
                   col,
                   blocks[row_block_id].size,
                   blocks[col_block_id].size) =
            (row_block_id + 1) * (col_block_id + 1) *
                Matrix::Ones(blocks[row_block_id].size,
                             blocks[col_block_id].size) +
            Matrix::Identity(blocks[row_block_id].size,
                             blocks[row_block_id].size);
      }
    }
  }

 protected:
  ContextImpl context_;
  int num_nonzeros_;
  std::unique_ptr<BlockRandomAccessDiagonalMatrix> m_;
};

TEST_F(BlockRandomAccessDiagonalMatrixTest, MatrixContents) {
  auto* crsm = m_->matrix();
  EXPECT_EQ(crsm->num_nonzeros(), num_nonzeros_);

  Matrix dense;
  crsm->ToDenseMatrix(&dense);

  double kTolerance = 1e-14;

  // (0,0)
  EXPECT_NEAR(
      (dense.block(0, 0, 3, 3) - (Matrix::Ones(3, 3) + Matrix::Identity(3, 3)))
          .norm(),
      0.0,
      kTolerance);

  // (1,1)
  EXPECT_NEAR((dense.block(3, 3, 4, 4) -
               (2 * 2 * Matrix::Ones(4, 4) + Matrix::Identity(4, 4)))
                  .norm(),
              0.0,
              kTolerance);

  // (1,1)
  EXPECT_NEAR((dense.block(7, 7, 5, 5) -
               (3 * 3 * Matrix::Ones(5, 5) + Matrix::Identity(5, 5)))
                  .norm(),
              0.0,
              kTolerance);

  // There is nothing else in the matrix besides these four blocks.
  EXPECT_NEAR(
      dense.norm(),
      sqrt(6 * 1.0 + 3 * 4.0 + 12 * 16.0 + 4 * 25.0 + 20 * 81.0 + 5 * 100.0),
      kTolerance);
}

TEST_F(BlockRandomAccessDiagonalMatrixTest, RightMultiplyAndAccumulate) {
  double kTolerance = 1e-14;
  auto* crsm = m_->matrix();
  Matrix dense;
  crsm->ToDenseMatrix(&dense);
  Vector x = Vector::Random(dense.rows());
  Vector expected_y = dense * x;
  Vector actual_y = Vector::Zero(dense.rows());
  m_->RightMultiplyAndAccumulate(x.data(), actual_y.data());
  EXPECT_NEAR((expected_y - actual_y).norm(), 0, kTolerance);
}

TEST_F(BlockRandomAccessDiagonalMatrixTest, Invert) {
  double kTolerance = 1e-14;
  auto* crsm = m_->matrix();
  Matrix dense;
  crsm->ToDenseMatrix(&dense);
  Matrix expected_inverse =
      dense.llt().solve(Matrix::Identity(dense.rows(), dense.rows()));

  m_->Invert();
  crsm->ToDenseMatrix(&dense);

  EXPECT_NEAR((expected_inverse - dense).norm(), 0.0, kTolerance);
}

}  // namespace ceres::internal
