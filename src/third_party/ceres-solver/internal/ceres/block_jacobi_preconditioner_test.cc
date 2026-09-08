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

#include "ceres/block_jacobi_preconditioner.h"

#include <memory>
#include <random>
#include <vector>

#include "Eigen/Dense"
#include "ceres/block_random_access_diagonal_matrix.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/linear_least_squares_problems.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(BlockSparseJacobiPreconditioner, _) {
  constexpr int kNumtrials = 10;
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_col_blocks = 3;
  options.min_col_block_size = 1;
  options.max_col_block_size = 3;

  options.num_row_blocks = 5;
  options.min_row_block_size = 1;
  options.max_row_block_size = 4;
  options.block_density = 0.25;
  std::mt19937 prng;

  Preconditioner::Options preconditioner_options;
  ContextImpl context;
  preconditioner_options.context = &context;

  for (int trial = 0; trial < kNumtrials; ++trial) {
    auto jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);
    Vector diagonal = Vector::Ones(jacobian->num_cols());
    Matrix dense_jacobian;
    jacobian->ToDenseMatrix(&dense_jacobian);
    Matrix hessian = dense_jacobian.transpose() * dense_jacobian;
    hessian.diagonal() += diagonal.array().square().matrix();

    BlockSparseJacobiPreconditioner pre(preconditioner_options, *jacobian);
    pre.Update(*jacobian, diagonal.data());

    // The const_cast is needed to be able to call GetCell.
    auto* m = const_cast<BlockRandomAccessDiagonalMatrix*>(&pre.matrix());
    EXPECT_EQ(m->num_rows(), jacobian->num_cols());
    EXPECT_EQ(m->num_cols(), jacobian->num_cols());

    const CompressedRowBlockStructure* bs = jacobian->block_structure();
    for (int i = 0; i < bs->cols.size(); ++i) {
      const int block_size = bs->cols[i].size;
      int r, c, row_stride, col_stride;
      CellInfo* cell_info = m->GetCell(i, i, &r, &c, &row_stride, &col_stride);
      Matrix actual_block_inverse =
          MatrixRef(cell_info->values, row_stride, col_stride)
              .block(r, c, block_size, block_size);
      Matrix expected_block = hessian.block(
          bs->cols[i].position, bs->cols[i].position, block_size, block_size);
      const double residual = (actual_block_inverse * expected_block -
                               Matrix::Identity(block_size, block_size))
                                  .norm();
      EXPECT_NEAR(residual, 0.0, 1e-12) << "Block: " << i;
    }
    options.num_col_blocks++;
    options.num_row_blocks++;
  }
}

TEST(CompressedRowSparseJacobiPreconditioner, _) {
  constexpr int kNumtrials = 10;
  CompressedRowSparseMatrix::RandomMatrixOptions options;
  options.num_col_blocks = 3;
  options.min_col_block_size = 1;
  options.max_col_block_size = 3;

  options.num_row_blocks = 5;
  options.min_row_block_size = 1;
  options.max_row_block_size = 4;
  options.block_density = 0.25;
  std::mt19937 prng;

  Preconditioner::Options preconditioner_options;
  ContextImpl context;
  preconditioner_options.context = &context;

  for (int trial = 0; trial < kNumtrials; ++trial) {
    auto jacobian =
        CompressedRowSparseMatrix::CreateRandomMatrix(options, prng);
    Vector diagonal = Vector::Ones(jacobian->num_cols());

    Matrix dense_jacobian;
    jacobian->ToDenseMatrix(&dense_jacobian);
    Matrix hessian = dense_jacobian.transpose() * dense_jacobian;
    hessian.diagonal() += diagonal.array().square().matrix();

    BlockCRSJacobiPreconditioner pre(preconditioner_options, *jacobian);
    pre.Update(*jacobian, diagonal.data());
    auto& m = pre.matrix();

    EXPECT_EQ(m.num_rows(), jacobian->num_cols());
    EXPECT_EQ(m.num_cols(), jacobian->num_cols());

    const auto& col_blocks = jacobian->col_blocks();
    for (int i = 0, col = 0; i < col_blocks.size(); ++i) {
      const int block_size = col_blocks[i].size;
      int idx = m.rows()[col];
      for (int j = 0; j < block_size; ++j) {
        EXPECT_EQ(m.rows()[col + j + 1] - m.rows()[col + j], block_size);
        for (int k = 0; k < block_size; ++k, ++idx) {
          EXPECT_EQ(m.cols()[idx], col + k);
        }
      }

      ConstMatrixRef actual_block_inverse(
          m.values() + m.rows()[col], block_size, block_size);
      Matrix expected_block = hessian.block(col, col, block_size, block_size);
      const double residual = (actual_block_inverse * expected_block -
                               Matrix::Identity(block_size, block_size))
                                  .norm();
      EXPECT_NEAR(residual, 0.0, 1e-12) << "Block: " << i;
      col += block_size;
    }
    options.num_col_blocks++;
    options.num_row_blocks++;
  }
}

}  // namespace ceres::internal
