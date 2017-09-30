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

#include "ceres/block_sparse_matrix.h"

#include <string>
#include "ceres/casts.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/triplet_sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

class BlockSparseMatrixTest : public ::testing::Test {
 protected :
  virtual void SetUp() {
    scoped_ptr<LinearLeastSquaresProblem> problem(
        CreateLinearLeastSquaresProblemFromId(2));
    CHECK_NOTNULL(problem.get());
    A_.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));

    problem.reset(CreateLinearLeastSquaresProblemFromId(1));
    CHECK_NOTNULL(problem.get());
    B_.reset(down_cast<TripletSparseMatrix*>(problem->A.release()));

    CHECK_EQ(A_->num_rows(), B_->num_rows());
    CHECK_EQ(A_->num_cols(), B_->num_cols());
    CHECK_EQ(A_->num_nonzeros(), B_->num_nonzeros());
  }

  scoped_ptr<BlockSparseMatrix> A_;
  scoped_ptr<TripletSparseMatrix> B_;
};

TEST_F(BlockSparseMatrixTest, SetZeroTest) {
  A_->SetZero();
  EXPECT_EQ(13, A_->num_nonzeros());
}

TEST_F(BlockSparseMatrixTest, RightMultiplyTest) {
  Vector y_a = Vector::Zero(A_->num_rows());
  Vector y_b = Vector::Zero(A_->num_rows());
  for (int i = 0; i < A_->num_cols(); ++i) {
    Vector x = Vector::Zero(A_->num_cols());
    x[i] = 1.0;
    A_->RightMultiply(x.data(), y_a.data());
    B_->RightMultiply(x.data(), y_b.data());
    EXPECT_LT((y_a - y_b).norm(), 1e-12);
  }
}

TEST_F(BlockSparseMatrixTest, LeftMultiplyTest) {
  Vector y_a = Vector::Zero(A_->num_cols());
  Vector y_b = Vector::Zero(A_->num_cols());
  for (int i = 0; i < A_->num_rows(); ++i) {
    Vector x = Vector::Zero(A_->num_rows());
    x[i] = 1.0;
    A_->LeftMultiply(x.data(), y_a.data());
    B_->LeftMultiply(x.data(), y_b.data());
    EXPECT_LT((y_a - y_b).norm(), 1e-12);
  }
}

TEST_F(BlockSparseMatrixTest, SquaredColumnNormTest) {
  Vector y_a = Vector::Zero(A_->num_cols());
  Vector y_b = Vector::Zero(A_->num_cols());
  A_->SquaredColumnNorm(y_a.data());
  B_->SquaredColumnNorm(y_b.data());
  EXPECT_LT((y_a - y_b).norm(), 1e-12);
}

TEST_F(BlockSparseMatrixTest, ToDenseMatrixTest) {
  Matrix m_a;
  Matrix m_b;
  A_->ToDenseMatrix(&m_a);
  B_->ToDenseMatrix(&m_b);
  EXPECT_LT((m_a - m_b).norm(), 1e-12);
}

TEST_F(BlockSparseMatrixTest, AppendRows) {
  scoped_ptr<LinearLeastSquaresProblem> problem(
      CreateLinearLeastSquaresProblemFromId(2));
  scoped_ptr<BlockSparseMatrix> m(
      down_cast<BlockSparseMatrix*>(problem->A.release()));
  A_->AppendRows(*m);
  EXPECT_EQ(A_->num_rows(), 2 * m->num_rows());
  EXPECT_EQ(A_->num_cols(), m->num_cols());

  problem.reset(CreateLinearLeastSquaresProblemFromId(1));
  scoped_ptr<TripletSparseMatrix> m2(
      down_cast<TripletSparseMatrix*>(problem->A.release()));
  B_->AppendRows(*m2);

  Vector y_a = Vector::Zero(A_->num_rows());
  Vector y_b = Vector::Zero(A_->num_rows());
  for (int i = 0; i < A_->num_cols(); ++i) {
    Vector x = Vector::Zero(A_->num_cols());
    x[i] = 1.0;
    y_a.setZero();
    y_b.setZero();

    A_->RightMultiply(x.data(), y_a.data());
    B_->RightMultiply(x.data(), y_b.data());
    EXPECT_LT((y_a - y_b).norm(), 1e-12);
  }
}

TEST_F(BlockSparseMatrixTest, AppendAndDeleteBlockDiagonalMatrix) {
  const std::vector<Block>& column_blocks = A_->block_structure()->cols;
  const int num_cols =
      column_blocks.back().size + column_blocks.back().position;
  Vector diagonal(num_cols);
  for (int i = 0; i < num_cols; ++i) {
    diagonal(i) = 2 * i * i + 1;
  }
  scoped_ptr<BlockSparseMatrix> appendage(
      BlockSparseMatrix::CreateDiagonalMatrix(diagonal.data(), column_blocks));

  A_->AppendRows(*appendage);
  Vector y_a, y_b;
  y_a.resize(A_->num_rows());
  y_b.resize(A_->num_rows());
  for (int i = 0; i < A_->num_cols(); ++i) {
    Vector x = Vector::Zero(A_->num_cols());
    x[i] = 1.0;
    y_a.setZero();
    y_b.setZero();

    A_->RightMultiply(x.data(), y_a.data());
    B_->RightMultiply(x.data(), y_b.data());
    EXPECT_LT((y_a.head(B_->num_rows()) - y_b.head(B_->num_rows())).norm(), 1e-12);
    Vector expected_tail = Vector::Zero(A_->num_cols());
    expected_tail(i) = diagonal(i);
    EXPECT_LT((y_a.tail(A_->num_cols()) - expected_tail).norm(), 1e-12);
  }


  A_->DeleteRowBlocks(column_blocks.size());
  EXPECT_EQ(A_->num_rows(), B_->num_rows());
  EXPECT_EQ(A_->num_cols(), B_->num_cols());

  y_a.resize(A_->num_rows());
  y_b.resize(A_->num_rows());
  for (int i = 0; i < A_->num_cols(); ++i) {
    Vector x = Vector::Zero(A_->num_cols());
    x[i] = 1.0;
    y_a.setZero();
    y_b.setZero();

    A_->RightMultiply(x.data(), y_a.data());
    B_->RightMultiply(x.data(), y_b.data());
    EXPECT_LT((y_a - y_b).norm(), 1e-12);
  }
}

TEST(BlockSparseMatrix, CreateDiagonalMatrix) {
  std::vector<Block> column_blocks;
  column_blocks.push_back(Block(2, 0));
  column_blocks.push_back(Block(1, 2));
  column_blocks.push_back(Block(3, 3));
  const int num_cols =
      column_blocks.back().size + column_blocks.back().position;
  Vector diagonal(num_cols);
  for (int i = 0; i < num_cols; ++i) {
    diagonal(i) = 2 * i * i + 1;
  }

  scoped_ptr<BlockSparseMatrix> m(
      BlockSparseMatrix::CreateDiagonalMatrix(diagonal.data(), column_blocks));
  const CompressedRowBlockStructure* bs = m->block_structure();
  EXPECT_EQ(bs->cols.size(), column_blocks.size());
  for (int i = 0; i < column_blocks.size(); ++i) {
    EXPECT_EQ(bs->cols[i].size, column_blocks[i].size);
    EXPECT_EQ(bs->cols[i].position, column_blocks[i].position);
  }
  EXPECT_EQ(m->num_rows(), m->num_cols());
  Vector x = Vector::Ones(num_cols);
  Vector y = Vector::Zero(num_cols);
  m->RightMultiply(x.data(), y.data());
  for (int i = 0; i < num_cols; ++i) {
    EXPECT_NEAR(y[i], diagonal[i], std::numeric_limits<double>::epsilon());
  }
}


}  // namespace internal
}  // namespace ceres
