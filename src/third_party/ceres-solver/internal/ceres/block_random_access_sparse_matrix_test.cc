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

#include "ceres/block_random_access_sparse_matrix.h"

#include <limits>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "ceres/internal/eigen.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(BlockRandomAccessSparseMatrix, GetCell) {
  ContextImpl context;
  constexpr int num_threads = 1;
  std::vector<Block> blocks;
  blocks.emplace_back(3, 0);
  blocks.emplace_back(4, 3);
  blocks.emplace_back(5, 7);
  constexpr int num_rows = 3 + 4 + 5;

  std::set<std::pair<int, int>> block_pairs;
  int num_nonzeros = 0;
  block_pairs.emplace(0, 0);
  num_nonzeros += blocks[0].size * blocks[0].size;

  block_pairs.emplace(1, 1);
  num_nonzeros += blocks[1].size * blocks[1].size;

  block_pairs.emplace(1, 2);
  num_nonzeros += blocks[1].size * blocks[2].size;

  block_pairs.emplace(0, 2);
  num_nonzeros += blocks[2].size * blocks[0].size;

  BlockRandomAccessSparseMatrix m(blocks, block_pairs, &context, num_threads);
  EXPECT_EQ(m.num_rows(), num_rows);
  EXPECT_EQ(m.num_cols(), num_rows);

  for (const auto& block_pair : block_pairs) {
    const int row_block_id = block_pair.first;
    const int col_block_id = block_pair.second;
    int row;
    int col;
    int row_stride;
    int col_stride;
    CellInfo* cell = m.GetCell(
        row_block_id, col_block_id, &row, &col, &row_stride, &col_stride);
    EXPECT_TRUE(cell != nullptr);
    EXPECT_EQ(row, 0);
    EXPECT_EQ(col, 0);
    EXPECT_EQ(row_stride, blocks[row_block_id].size);
    EXPECT_EQ(col_stride, blocks[col_block_id].size);

    // Write into the block
    MatrixRef(cell->values, row_stride, col_stride)
        .block(row, col, blocks[row_block_id].size, blocks[col_block_id].size) =
        (row_block_id + 1) * (col_block_id + 1) *
        Matrix::Ones(blocks[row_block_id].size, blocks[col_block_id].size);
  }

  const BlockSparseMatrix* bsm = m.matrix();
  EXPECT_EQ(bsm->num_nonzeros(), num_nonzeros);

  Matrix dense;
  bsm->ToDenseMatrix(&dense);

  double kTolerance = 1e-14;

  // (0, 0)
  EXPECT_NEAR(
      (dense.block(0, 0, 3, 3) - Matrix::Ones(3, 3)).norm(), 0.0, kTolerance);
  // (1, 1)
  EXPECT_NEAR((dense.block(3, 3, 4, 4) - 2 * 2 * Matrix::Ones(4, 4)).norm(),
              0.0,
              kTolerance);
  // (1, 2)
  EXPECT_NEAR((dense.block(3, 3 + 4, 4, 5) - 2 * 3 * Matrix::Ones(4, 5)).norm(),
              0.0,
              kTolerance);
  // (0, 2)
  EXPECT_NEAR((dense.block(0, 3 + 4, 3, 5) - 3 * 1 * Matrix::Ones(3, 5)).norm(),
              0.0,
              kTolerance);

  // There is nothing else in the matrix besides these four blocks.
  EXPECT_NEAR(
      dense.norm(), sqrt(9. + 16. * 16. + 36. * 20. + 9. * 15.), kTolerance);

  Vector x = Vector::Ones(dense.rows());
  Vector actual_y = Vector::Zero(dense.rows());
  Vector expected_y = Vector::Zero(dense.rows());

  expected_y += dense.selfadjointView<Eigen::Upper>() * x;
  m.SymmetricRightMultiplyAndAccumulate(x.data(), actual_y.data());
  EXPECT_NEAR((expected_y - actual_y).norm(), 0.0, kTolerance)
      << "actual: " << actual_y.transpose() << "\n"
      << "expected: " << expected_y.transpose() << "matrix: \n " << dense;
}

// IntPairToInt64 is private, thus this fixture is needed to access and
// test it.
class BlockRandomAccessSparseMatrixTest : public ::testing::Test {
 public:
  void SetUp() final {
    std::vector<Block> blocks;
    blocks.emplace_back(1, 0);
    std::set<std::pair<int, int>> block_pairs;
    block_pairs.emplace(0, 0);
    m_ = std::make_unique<BlockRandomAccessSparseMatrix>(
        blocks, block_pairs, &context_, 1);
  }

  void CheckIntPairToInt64(int a, int b) {
    int64_t value = m_->IntPairToInt64(a, b);
    EXPECT_GT(value, 0) << "Overflow a = " << a << " b = " << b;
    EXPECT_GT(value, a) << "Overflow a = " << a << " b = " << b;
    EXPECT_GT(value, b) << "Overflow a = " << a << " b = " << b;
  }

  void CheckInt64ToIntPair() {
    uint64_t max_rows = m_->kRowShift;
    for (int row = max_rows - 10; row < max_rows; ++row) {
      for (int col = 0; col < 10; ++col) {
        int row_computed;
        int col_computed;
        m_->Int64ToIntPair(
            m_->IntPairToInt64(row, col), &row_computed, &col_computed);
        EXPECT_EQ(row, row_computed);
        EXPECT_EQ(col, col_computed);
      }
    }
  }

 private:
  ContextImpl context_;
  std::unique_ptr<BlockRandomAccessSparseMatrix> m_;
};

TEST_F(BlockRandomAccessSparseMatrixTest, IntPairToInt64Overflow) {
  CheckIntPairToInt64(std::numeric_limits<int32_t>::max(),
                      std::numeric_limits<int32_t>::max());
}

TEST_F(BlockRandomAccessSparseMatrixTest, Int64ToIntPair) {
  CheckInt64ToIntPair();
}

}  // namespace ceres::internal
