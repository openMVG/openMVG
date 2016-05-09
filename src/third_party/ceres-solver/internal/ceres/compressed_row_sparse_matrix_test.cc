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

#include "ceres/compressed_row_sparse_matrix.h"

#include <numeric>
#include "ceres/casts.h"
#include "ceres/crs_matrix.h"
#include "ceres/cxsparse.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/random.h"
#include "ceres/triplet_sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

using std::vector;

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

    vector<int>* row_blocks = crsm->mutable_row_blocks();
    row_blocks->resize(num_rows);
    std::fill(row_blocks->begin(), row_blocks->end(), 1);

    vector<int>* col_blocks = crsm->mutable_col_blocks();
    col_blocks->resize(num_cols);
    std::fill(col_blocks->begin(), col_blocks->end(), 1);
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
  // Clear the row and column blocks as these are purely scalar tests.
  crsm->mutable_row_blocks()->clear();
  crsm->mutable_col_blocks()->clear();
  for (int i = 0; i < num_rows; ++i) {
    tsm->Resize(num_rows - i, num_cols);
    crsm->DeleteRows(crsm->num_rows() - tsm->num_rows());
    CompareMatrices(tsm.get(), crsm.get());
  }
}

TEST_F(CompressedRowSparseMatrixTest, AppendRows) {
  // Clear the row and column blocks as these are purely scalar tests.
  crsm->mutable_row_blocks()->clear();
  crsm->mutable_col_blocks()->clear();

  for (int i = 0; i < num_rows; ++i) {
    TripletSparseMatrix tsm_appendage(*tsm);
    tsm_appendage.Resize(i, num_cols);

    tsm->AppendRows(tsm_appendage);
    CompressedRowSparseMatrix crsm_appendage(tsm_appendage);
    crsm->AppendRows(crsm_appendage);

    CompareMatrices(tsm.get(), crsm.get());
  }
}

TEST_F(CompressedRowSparseMatrixTest, AppendAndDeleteBlockDiagonalMatrix) {
  int num_diagonal_rows = crsm->num_cols();

  scoped_array<double> diagonal(new double[num_diagonal_rows]);
  for (int i = 0; i < num_diagonal_rows; ++i) {
    diagonal[i] = i;
  }

  vector<int> row_and_column_blocks;
  row_and_column_blocks.push_back(1);
  row_and_column_blocks.push_back(2);
  row_and_column_blocks.push_back(2);

  const vector<int> pre_row_blocks = crsm->row_blocks();
  const vector<int> pre_col_blocks = crsm->col_blocks();

  scoped_ptr<CompressedRowSparseMatrix> appendage(
      CompressedRowSparseMatrix::CreateBlockDiagonalMatrix(
          diagonal.get(), row_and_column_blocks));
  LOG(INFO) << appendage->row_blocks().size();

  crsm->AppendRows(*appendage);

  const vector<int> post_row_blocks = crsm->row_blocks();
  const vector<int> post_col_blocks = crsm->col_blocks();

  vector<int> expected_row_blocks = pre_row_blocks;
  expected_row_blocks.insert(expected_row_blocks.end(),
                             row_and_column_blocks.begin(),
                             row_and_column_blocks.end());

  vector<int> expected_col_blocks = pre_col_blocks;

  EXPECT_EQ(expected_row_blocks, crsm->row_blocks());
  EXPECT_EQ(expected_col_blocks, crsm->col_blocks());

  crsm->DeleteRows(num_diagonal_rows);
  EXPECT_EQ(crsm->row_blocks(), pre_row_blocks);
  EXPECT_EQ(crsm->col_blocks(), pre_col_blocks);
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

  // Block structure:
  //  A  A  A  A  B  B
  //  A  A  A  A  B  B
  //  A  A  A  A  B  B
  //  C  C  C  C  D  D
  //  C  C  C  C  D  D
  //  C  C  C  C  D  D

  CompressedRowSparseMatrix matrix(5, 6, 30);
  int* rows = matrix.mutable_rows();
  int* cols = matrix.mutable_cols();
  double* values = matrix.mutable_values();
  matrix.mutable_row_blocks()->push_back(3);
  matrix.mutable_row_blocks()->push_back(3);
  matrix.mutable_col_blocks()->push_back(4);
  matrix.mutable_col_blocks()->push_back(2);

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

  std::copy(values, values + 17, cols);

  scoped_ptr<CompressedRowSparseMatrix> transpose(matrix.Transpose());

  ASSERT_EQ(transpose->row_blocks().size(), matrix.col_blocks().size());
  for (int i = 0; i < transpose->row_blocks().size(); ++i) {
    EXPECT_EQ(transpose->row_blocks()[i], matrix.col_blocks()[i]);
  }

  ASSERT_EQ(transpose->col_blocks().size(), matrix.row_blocks().size());
  for (int i = 0; i < transpose->col_blocks().size(); ++i) {
    EXPECT_EQ(transpose->col_blocks()[i], matrix.row_blocks()[i]);
  }

  Matrix dense_matrix;
  matrix.ToDenseMatrix(&dense_matrix);

  Matrix dense_transpose;
  transpose->ToDenseMatrix(&dense_transpose);
  EXPECT_NEAR((dense_matrix - dense_transpose.transpose()).norm(), 0.0, 1e-14);
}

#ifndef CERES_NO_CXSPARSE

struct RandomMatrixOptions {
  int num_row_blocks;
  int min_row_block_size;
  int max_row_block_size;
  int num_col_blocks;
  int min_col_block_size;
  int max_col_block_size;
  double block_density;
};

CompressedRowSparseMatrix* CreateRandomCompressedRowSparseMatrix(
    const RandomMatrixOptions& options) {
  vector<int> row_blocks;
  for (int i = 0; i < options.num_row_blocks; ++i) {
    const int delta_block_size =
        Uniform(options.max_row_block_size - options.min_row_block_size);
    row_blocks.push_back(options.min_row_block_size + delta_block_size);
  }

  vector<int> col_blocks;
  for (int i = 0; i < options.num_col_blocks; ++i) {
    const int delta_block_size =
        Uniform(options.max_col_block_size - options.min_col_block_size);
    col_blocks.push_back(options.min_col_block_size + delta_block_size);
  }

  vector<int> rows;
  vector<int> cols;
  vector<double> values;

  while (values.size() == 0) {
    int row_block_begin = 0;
    for (int r = 0; r < options.num_row_blocks; ++r) {
      int col_block_begin = 0;
      for (int c = 0; c < options.num_col_blocks; ++c) {
        if (RandDouble() <= options.block_density) {
          for (int i = 0; i < row_blocks[r]; ++i) {
            for (int j = 0; j < col_blocks[c]; ++j) {
              rows.push_back(row_block_begin + i);
              cols.push_back(col_block_begin + j);
              values.push_back(RandNormal());
            }
          }
        }
        col_block_begin += col_blocks[c];
      }
      row_block_begin += row_blocks[r];
    }
  }

  const int num_rows = std::accumulate(row_blocks.begin(), row_blocks.end(), 0);
  const int num_cols = std::accumulate(col_blocks.begin(), col_blocks.end(), 0);
  const int num_nonzeros = values.size();

  TripletSparseMatrix tsm(num_rows, num_cols, num_nonzeros);
  std::copy(rows.begin(), rows.end(), tsm.mutable_rows());
  std::copy(cols.begin(), cols.end(), tsm.mutable_cols());
  std::copy(values.begin(), values.end(), tsm.mutable_values());
  tsm.set_num_nonzeros(num_nonzeros);
  CompressedRowSparseMatrix* matrix = new CompressedRowSparseMatrix(tsm);
  (*matrix->mutable_row_blocks())  = row_blocks;
  (*matrix->mutable_col_blocks())  = col_blocks;
  return matrix;
}

void ToDenseMatrix(const cs_di* matrix, Matrix* dense_matrix) {
  dense_matrix->resize(matrix->m, matrix->n);
  dense_matrix->setZero();

  for (int c = 0; c < matrix->n; ++c) {
    for (int idx = matrix->p[c]; idx < matrix->p[c + 1]; ++idx) {
      const int r = matrix->i[idx];
      (*dense_matrix)(r, c) = matrix->x[idx];
    }
  }
}

TEST(CompressedRowSparseMatrix, ComputeOuterProduct) {
  // "Randomly generated seed."
  SetRandomState(29823);
  int kMaxNumRowBlocks = 10;
  int kMaxNumColBlocks = 10;
  int kNumTrials = 10;

  CXSparse cxsparse;
  const double kTolerance = 1e-18;

  // Create a random matrix, compute its outer product using CXSParse
  // and ComputeOuterProduct. Convert both matrices to dense matrices
  // and compare their upper triangular parts. They should be within
  // kTolerance of each other.
  for (int num_row_blocks = 1;
       num_row_blocks < kMaxNumRowBlocks;
       ++num_row_blocks) {
    for (int num_col_blocks = 1;
         num_col_blocks < kMaxNumColBlocks;
         ++num_col_blocks) {
      for (int trial = 0; trial < kNumTrials; ++trial) {
        RandomMatrixOptions options;
        options.num_row_blocks = num_row_blocks;
        options.num_col_blocks = num_col_blocks;
        options.min_row_block_size = 1;
        options.max_row_block_size = 5;
        options.min_col_block_size = 1;
        options.max_col_block_size = 10;
        options.block_density = std::max(0.1, RandDouble());

        VLOG(2) << "num row blocks: " << options.num_row_blocks;
        VLOG(2) << "num col blocks: " << options.num_col_blocks;
        VLOG(2) << "min row block size: " << options.min_row_block_size;
        VLOG(2) << "max row block size: " << options.max_row_block_size;
        VLOG(2) << "min col block size: " << options.min_col_block_size;
        VLOG(2) << "max col block size: " << options.max_col_block_size;
        VLOG(2) << "block density: " << options.block_density;

        scoped_ptr<CompressedRowSparseMatrix> matrix(
            CreateRandomCompressedRowSparseMatrix(options));

        cs_di cs_matrix_transpose =
            cxsparse.CreateSparseMatrixTransposeView(matrix.get());
        cs_di* cs_matrix = cxsparse.TransposeMatrix(&cs_matrix_transpose);
        cs_di* expected_outer_product =
            cxsparse.MatrixMatrixMultiply(&cs_matrix_transpose, cs_matrix);

        vector<int> program;
        scoped_ptr<CompressedRowSparseMatrix> outer_product(
            CompressedRowSparseMatrix::CreateOuterProductMatrixAndProgram(
                *matrix, &program));
        CompressedRowSparseMatrix::ComputeOuterProduct(*matrix,
                                                       program,
                                                       outer_product.get());

        cs_di actual_outer_product =
            cxsparse.CreateSparseMatrixTransposeView(outer_product.get());

        ASSERT_EQ(actual_outer_product.m, actual_outer_product.n);
        ASSERT_EQ(expected_outer_product->m, expected_outer_product->n);
        ASSERT_EQ(actual_outer_product.m, expected_outer_product->m);

        Matrix actual_matrix;
        Matrix expected_matrix;

        ToDenseMatrix(expected_outer_product, &expected_matrix);
        expected_matrix.triangularView<Eigen::StrictlyLower>().setZero();

        ToDenseMatrix(&actual_outer_product, &actual_matrix);
        const double diff_norm =
            (actual_matrix - expected_matrix).norm() / expected_matrix.norm();
        ASSERT_NEAR(diff_norm, 0.0, kTolerance)
            << "expected: \n"
            << expected_matrix
            << "\nactual: \n"
            << actual_matrix;

        cxsparse.Free(cs_matrix);
        cxsparse.Free(expected_outer_product);
      }
    }
  }
}

#endif  // CERES_NO_CXSPARSE

}  // namespace internal
}  // namespace ceres
