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

#include "ceres/compressed_row_sparse_matrix.h"

#include "ceres/casts.h"
#include "ceres/crs_matrix.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/triplet_sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

void CompareMatrices(const SparseMatrix* a, const SparseMatrix* b) {
  EXPECT_EQ(a->num_rows(), b->num_rows());
  EXPECT_EQ(a->num_cols(), b->num_cols());

  int num_rows = a->num_rows();
  int num_cols = a->num_cols();

  for (int i = 0; i < num_cols; ++i) {
    Vector x = Vector::Zero(num_cols);
    x(i) = 1.0;

    Vector y_a = Vector::Zero(num_rows);
    Vector y_b = Vector::Zero(num_rows);

    a->RightMultiply(x.data(), y_a.data());
    b->RightMultiply(x.data(), y_b.data());

    EXPECT_EQ((y_a - y_b).norm(), 0);
  }
}

class CompressedRowSparseMatrixTest : public ::testing::Test {
 protected :
  virtual void SetUp() {
    scoped_ptr<LinearLeastSquaresProblem> problem(
        CreateLinearLeastSquaresProblemFromId(1));

    CHECK_NOTNULL(problem.get());

    tsm.reset(down_cast<TripletSparseMatrix*>(problem->A.release()));
    crsm.reset(new CompressedRowSparseMatrix(*tsm));

    num_rows = tsm->num_rows();
    num_cols = tsm->num_cols();
  }

  int num_rows;
  int num_cols;

  scoped_ptr<TripletSparseMatrix> tsm;
  scoped_ptr<CompressedRowSparseMatrix> crsm;
};

TEST_F(CompressedRowSparseMatrixTest, RightMultiply) {
  CompareMatrices(tsm.get(), crsm.get());
}

TEST_F(CompressedRowSparseMatrixTest, LeftMultiply) {
  for (int i = 0; i < num_rows; ++i) {
    Vector a = Vector::Zero(num_rows);
    a(i) = 1.0;

    Vector b1 = Vector::Zero(num_cols);
    Vector b2 = Vector::Zero(num_cols);

    tsm->LeftMultiply(a.data(), b1.data());
    crsm->LeftMultiply(a.data(), b2.data());

    EXPECT_EQ((b1 - b2).norm(), 0);
  }
}

TEST_F(CompressedRowSparseMatrixTest, ColumnNorm) {
  Vector b1 = Vector::Zero(num_cols);
  Vector b2 = Vector::Zero(num_cols);

  tsm->SquaredColumnNorm(b1.data());
  crsm->SquaredColumnNorm(b2.data());

  EXPECT_EQ((b1 - b2).norm(), 0);
}

TEST_F(CompressedRowSparseMatrixTest, Scale) {
  Vector scale(num_cols);
  for (int i = 0; i < num_cols; ++i) {
    scale(i) = i + 1;
  }

  tsm->ScaleColumns(scale.data());
  crsm->ScaleColumns(scale.data());
  CompareMatrices(tsm.get(), crsm.get());
}

TEST_F(CompressedRowSparseMatrixTest, DeleteRows) {
  for (int i = 0; i < num_rows; ++i) {
    tsm->Resize(num_rows - i, num_cols);
    crsm->DeleteRows(crsm->num_rows() - tsm->num_rows());
    CompareMatrices(tsm.get(), crsm.get());
  }
}

TEST_F(CompressedRowSparseMatrixTest, AppendRows) {
  for (int i = 0; i < num_rows; ++i) {
    TripletSparseMatrix tsm_appendage(*tsm);
    tsm_appendage.Resize(i, num_cols);

    tsm->AppendRows(tsm_appendage);
    CompressedRowSparseMatrix crsm_appendage(tsm_appendage);
    crsm->AppendRows(crsm_appendage);

    CompareMatrices(tsm.get(), crsm.get());
  }
}

TEST_F(CompressedRowSparseMatrixTest, ToDenseMatrix) {
  Matrix tsm_dense;
  Matrix crsm_dense;

  tsm->ToDenseMatrix(&tsm_dense);
  crsm->ToDenseMatrix(&crsm_dense);

  EXPECT_EQ((tsm_dense - crsm_dense).norm(), 0.0);
}

TEST_F(CompressedRowSparseMatrixTest, ToCRSMatrix) {
  CRSMatrix crs_matrix;
  crsm->ToCRSMatrix(&crs_matrix);
  EXPECT_EQ(crsm->num_rows(), crs_matrix.num_rows);
  EXPECT_EQ(crsm->num_cols(), crs_matrix.num_cols);
  EXPECT_EQ(crsm->num_rows() + 1, crs_matrix.rows.size());
  EXPECT_EQ(crsm->num_nonzeros(), crs_matrix.cols.size());
  EXPECT_EQ(crsm->num_nonzeros(), crs_matrix.values.size());

  for (int i = 0; i < crsm->num_rows() + 1; ++i) {
    EXPECT_EQ(crsm->rows()[i], crs_matrix.rows[i]);
  }

  for (int i = 0; i < crsm->num_nonzeros(); ++i) {
    EXPECT_EQ(crsm->cols()[i], crs_matrix.cols[i]);
    EXPECT_EQ(crsm->values()[i], crs_matrix.values[i]);
  }
}

TEST(CompressedRowSparseMatrix, CreateBlockDiagonalMatrix) {
  vector<int> blocks;
  blocks.push_back(1);
  blocks.push_back(2);
  blocks.push_back(2);

  Vector diagonal(5);
  for (int i = 0; i < 5; ++i) {
    diagonal(i) = i + 1;
  }

  scoped_ptr<CompressedRowSparseMatrix> matrix(
      CompressedRowSparseMatrix::CreateBlockDiagonalMatrix(
          diagonal.data(), blocks));

  EXPECT_EQ(matrix->num_rows(), 5);
  EXPECT_EQ(matrix->num_cols(), 5);
  EXPECT_EQ(matrix->num_nonzeros(), 9);
  EXPECT_EQ(blocks, matrix->row_blocks());
  EXPECT_EQ(blocks, matrix->col_blocks());

  Vector x(5);
  Vector y(5);

  x.setOnes();
  y.setZero();
  matrix->RightMultiply(x.data(), y.data());
  for (int i = 0; i < diagonal.size(); ++i) {
    EXPECT_EQ(y[i], diagonal[i]);
  }

  y.setZero();
  matrix->LeftMultiply(x.data(), y.data());
  for (int i = 0; i < diagonal.size(); ++i) {
    EXPECT_EQ(y[i], diagonal[i]);
  }

  Matrix dense;
  matrix->ToDenseMatrix(&dense);
  EXPECT_EQ((dense.diagonal() - diagonal).norm(), 0.0);
}

class SolveLowerTriangularTest : public ::testing::Test {
 protected:
  void SetUp() {
    matrix_.reset(new CompressedRowSparseMatrix(4, 4, 7));
    int* rows = matrix_->mutable_rows();
    int* cols = matrix_->mutable_cols();
    double* values = matrix_->mutable_values();

    rows[0] = 0;
    cols[0] = 0;
    values[0] = 0.50754;

    rows[1] = 1;
    cols[1] = 1;
    values[1] = 0.80483;

    rows[2] = 2;
    cols[2] = 1;
    values[2] = 0.14120;
    cols[3] = 2;
    values[3] = 0.3;

    rows[3] = 4;
    cols[4] = 0;
    values[4] = 0.77696;
    cols[5] = 1;
    values[5] = 0.41860;
    cols[6] = 3;
    values[6] = 0.88979;

    rows[4] = 7;
  }

  scoped_ptr<CompressedRowSparseMatrix> matrix_;
};

TEST_F(SolveLowerTriangularTest, SolveInPlace) {
  double rhs_and_solution[] = {1.0, 1.0, 2.0, 2.0};
  double expected[] = {1.970288,  1.242498,  6.081864, -0.057255};
  matrix_->SolveLowerTriangularInPlace(rhs_and_solution);
  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(rhs_and_solution[i], expected[i], 1e-4) << i;
  }
}

TEST_F(SolveLowerTriangularTest, TransposeSolveInPlace) {
  double rhs_and_solution[] = {1.0, 1.0, 2.0, 2.0};
  const double expected[] = { -1.4706, -1.0962, 6.6667, 2.2477};

  matrix_->SolveLowerTriangularTransposeInPlace(rhs_and_solution);
  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(rhs_and_solution[i], expected[i], 1e-4) << i;
  }
}

TEST(CompressedRowSparseMatrix, Transpose) {
  //  0  1  0  2  3  0
  //  4  6  7  0  0  8
  //  9 10  0 11 12  0
  // 13  0 14 15  9  0
  //  0 16 17  0  0  0

  CompressedRowSparseMatrix matrix(5, 6, 30);
  int* rows = matrix.mutable_rows();
  int* cols = matrix.mutable_cols();
  double* values = matrix.mutable_values();

  rows[0] = 0;
  cols[0] = 1;
  cols[1] = 3;
  cols[2] = 4;

  rows[1] = 3;
  cols[3] = 0;
  cols[4] = 1;
  cols[5] = 2;
  cols[6] = 5;


  rows[2] = 7;
  cols[7] = 0;
  cols[8] = 1;
  cols[9] = 3;
  cols[10] = 4;

  rows[3] = 11;
  cols[11] = 0;
  cols[12] = 2;
  cols[13] = 3;
  cols[14] = 4;

  rows[4] = 15;
  cols[15] = 1;
  cols[16] = 2;
  rows[5] = 17;

  copy(values, values + 17, cols);

  scoped_ptr<CompressedRowSparseMatrix> transpose(matrix.Transpose());

  Matrix dense_matrix;
  matrix.ToDenseMatrix(&dense_matrix);

  Matrix dense_transpose;
  transpose->ToDenseMatrix(&dense_transpose);
  EXPECT_NEAR((dense_matrix - dense_transpose.transpose()).norm(), 0.0, 1e-14);
}

}  // namespace internal
}  // namespace ceres
