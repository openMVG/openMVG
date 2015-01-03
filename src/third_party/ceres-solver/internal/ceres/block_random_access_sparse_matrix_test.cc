// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2010, 2011, 2012 Google Inc. All rights reserved.
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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include <limits>
#include <vector>
#include "ceres/block_random_access_sparse_matrix.h"
#include "ceres/internal/eigen.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(BlockRandomAccessSparseMatrix, GetCell) {
  vector<int> blocks;
  blocks.push_back(3);
  blocks.push_back(4);
  blocks.push_back(5);
  const int num_rows = 3 + 4 + 5;

  set< pair<int, int> > block_pairs;
  int num_nonzeros = 0;
  block_pairs.insert(make_pair(0, 0));
  num_nonzeros += blocks[0] * blocks[0];

  block_pairs.insert(make_pair(1, 1));
  num_nonzeros += blocks[1] * blocks[1];

  block_pairs.insert(make_pair(1, 2));
  num_nonzeros += blocks[1] * blocks[2];

  block_pairs.insert(make_pair(0, 2));
  num_nonzeros += blocks[2] * blocks[0];

  BlockRandomAccessSparseMatrix m(blocks, block_pairs);
  EXPECT_EQ(m.num_rows(), num_rows);
  EXPECT_EQ(m.num_cols(), num_rows);

  for (set<pair<int, int> >::const_iterator it = block_pairs.begin();
       it != block_pairs.end();
       ++it) {
    const int row_block_id = it->first;
    const int col_block_id = it->second;
    int row;
    int col;
    int row_stride;
    int col_stride;
    CellInfo* cell =  m.GetCell(row_block_id, col_block_id,
                                &row, &col,
                                &row_stride, &col_stride);
    EXPECT_TRUE(cell != NULL);
    EXPECT_EQ(row, 0);
    EXPECT_EQ(col, 0);
    EXPECT_EQ(row_stride, blocks[row_block_id]);
    EXPECT_EQ(col_stride, blocks[col_block_id]);

    // Write into the block
    MatrixRef(cell->values, row_stride, col_stride).block(
        row, col, blocks[row_block_id], blocks[col_block_id]) =
        (row_block_id + 1) * (col_block_id +1) *
        Matrix::Ones(blocks[row_block_id], blocks[col_block_id]);
  }

  const TripletSparseMatrix* tsm = m.matrix();
  EXPECT_EQ(tsm->num_nonzeros(), num_nonzeros);
  EXPECT_EQ(tsm->max_num_nonzeros(), num_nonzeros);

  Matrix dense;
  tsm->ToDenseMatrix(&dense);

  double kTolerance = 1e-14;

  // (0, 0)
  EXPECT_NEAR((dense.block(0, 0, 3, 3) - Matrix::Ones(3, 3)).norm(),
              0.0,
              kTolerance);
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
  EXPECT_NEAR(dense.norm(), sqrt(9. + 16. * 16. + 36. * 20. + 9. * 15.),
              kTolerance);

  Vector x = Vector::Ones(dense.rows());
  Vector actual_y = Vector::Zero(dense.rows());
  Vector expected_y = Vector::Zero(dense.rows());

  expected_y += dense.selfadjointView<Eigen::Upper>() * x;
  m.SymmetricRightMultiply(x.data(), actual_y.data());
  EXPECT_NEAR((expected_y - actual_y).norm(), 0.0, kTolerance)
      << "actual: " << actual_y.transpose() << "\n"
      << "expected: " << expected_y.transpose()
      << "matrix: \n " << dense;
}

// IntPairToLong is private, thus this fixture is needed to access and
// test it.
class BlockRandomAccessSparseMatrixTest : public ::testing::Test {
 public:
  virtual void SetUp() {
    vector<int> blocks;
    blocks.push_back(1);
    set< pair<int, int> > block_pairs;
    block_pairs.insert(make_pair(0, 0));
    m_.reset(new BlockRandomAccessSparseMatrix(blocks, block_pairs));
  }

  void CheckIntPairToLong(int a, int b) {
    int64 value = m_->IntPairToLong(a, b);
    EXPECT_GT(value, 0) << "Overflow a = " << a << " b = " << b;
    EXPECT_GT(value, a) << "Overflow a = " << a << " b = " << b;
    EXPECT_GT(value, b) << "Overflow a = " << a << " b = " << b;
  }

  void CheckLongToIntPair() {
    uint64 max_rows =  m_->kMaxRowBlocks;
    for (int row = max_rows - 10; row < max_rows; ++row) {
      for (int col = 0; col < 10; ++col) {
        int row_computed;
        int col_computed;
        m_->LongToIntPair(m_->IntPairToLong(row, col),
                          &row_computed,
                          &col_computed);
        EXPECT_EQ(row, row_computed);
        EXPECT_EQ(col, col_computed);
      }
    }
  }

 private:
  scoped_ptr<BlockRandomAccessSparseMatrix> m_;
};

TEST_F(BlockRandomAccessSparseMatrixTest, IntPairToLongOverflow) {
  CheckIntPairToLong(numeric_limits<int>::max(), numeric_limits<int>::max());
}

TEST_F(BlockRandomAccessSparseMatrixTest, LongToIntPair) {
  CheckLongToIntPair();
}

}  // namespace internal
}  // namespace ceres
