// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2014 Google Inc. All rights reserved.
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

#include "ceres/block_jacobi_preconditioner.h"

#include <vector>
#include "ceres/block_random_access_diagonal_matrix.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/internal/scoped_ptr.h"
#include "gtest/gtest.h"
#include "Eigen/Dense"

namespace ceres {
namespace internal {


class BlockJacobiPreconditionerTest : public ::testing::Test {
 protected:
  void SetUpFromProblemId(int problem_id) {
    scoped_ptr<LinearLeastSquaresProblem> problem(
        CreateLinearLeastSquaresProblemFromId(problem_id));

    CHECK_NOTNULL(problem.get());
    A.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));
    D.reset(problem->D.release());

    Matrix dense_a;
    A->ToDenseMatrix(&dense_a);
    dense_ata = dense_a.transpose() * dense_a;
    dense_ata += VectorRef(D.get(), A->num_cols())
        .array().square().matrix().asDiagonal();
  }

  void VerifyDiagonalBlocks(const int problem_id) {
    SetUpFromProblemId(problem_id);

    BlockJacobiPreconditioner pre(*A);
    pre.Update(*A, D.get());
    BlockRandomAccessDiagonalMatrix* m =
        const_cast<BlockRandomAccessDiagonalMatrix*>(&pre.matrix());
    EXPECT_EQ(m->num_rows(), A->num_cols());
    EXPECT_EQ(m->num_cols(), A->num_cols());

    const CompressedRowBlockStructure* bs = A->block_structure();
    for (int i = 0; i < bs->cols.size(); ++i) {
      const int block_size = bs->cols[i].size;
      int r, c, row_stride, col_stride;
      CellInfo* cell_info = m->GetCell(i, i,
                                       &r, &c,
                                       &row_stride, &col_stride);
      MatrixRef m(cell_info->values, row_stride, col_stride);
      Matrix actual_block_inverse = m.block(r, c, block_size, block_size);
      Matrix expected_block = dense_ata.block(bs->cols[i].position,
                                              bs->cols[i].position,
                                              block_size,
                                              block_size);
      const double residual = (actual_block_inverse * expected_block -
                               Matrix::Identity(block_size, block_size)).norm();
      EXPECT_NEAR(residual, 0.0, 1e-12) << "Block: " << i;
    }
  }

  scoped_ptr<BlockSparseMatrix> A;
  scoped_array<double> D;
  Matrix dense_ata;
};

TEST_F(BlockJacobiPreconditionerTest, SmallProblem) {
  VerifyDiagonalBlocks(2);
}

TEST_F(BlockJacobiPreconditionerTest, LargeProblem) {
  VerifyDiagonalBlocks(3);
}

}  // namespace internal
}  // namespace ceres
