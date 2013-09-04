// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2013 Google Inc. All rights reserved.
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

#include "ceres/incomplete_lq_factorization.h"

#include "Eigen/Dense"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/internal/scoped_ptr.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

void ExpectMatricesAreEqual(const CompressedRowSparseMatrix& expected,
                            const CompressedRowSparseMatrix& actual,
                            const double tolerance) {
  EXPECT_EQ(expected.num_rows(), actual.num_rows());
  EXPECT_EQ(expected.num_cols(), actual.num_cols());
  for (int i = 0; i < expected.num_rows(); ++i) {
    EXPECT_EQ(expected.rows()[i], actual.rows()[i]);
  }

  for (int i = 0; i < actual.num_nonzeros(); ++i) {
    EXPECT_EQ(expected.cols()[i], actual.cols()[i]);
    EXPECT_NEAR(expected.values()[i], actual.values()[i], tolerance);
  }
}

TEST(IncompleteQRFactorization, OneByOneMatrix) {
  CompressedRowSparseMatrix matrix(1, 1, 1);
  matrix.mutable_rows()[0] = 0;
  matrix.mutable_rows()[1] = 1;
  matrix.mutable_cols()[0] = 0;
  matrix.mutable_values()[0] = 2;

  scoped_ptr<CompressedRowSparseMatrix> l(
      IncompleteLQFactorization(matrix, 1, 0.0, 1, 0.0));
  ExpectMatricesAreEqual(matrix, *l, 1e-16);
}

TEST(IncompleteLQFactorization, CompleteFactorization) {
  double dense_matrix[] = {
    0.00000,  0.00000,  0.20522,  0.00000,  0.49077,  0.92835,  0.00000,  0.83825,  0.00000,  0.00000,  // NOLINT
    0.00000,  0.00000,  0.00000,  0.62491,  0.38144,  0.00000,  0.79394,  0.79178,  0.00000,  0.44382,  // NOLINT
    0.00000,  0.00000,  0.00000,  0.61517,  0.55941,  0.00000,  0.00000,  0.00000,  0.00000,  0.60664,  // NOLINT
    0.00000,  0.00000,  0.00000,  0.00000,  0.45031,  0.00000,  0.64132,  0.00000,  0.38832,  0.00000,  // NOLINT
    0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.57121,  0.00000,  0.01375,  0.70640,  0.00000,  // NOLINT
    0.00000,  0.00000,  0.07451,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  // NOLINT
    0.68095,  0.00000,  0.00000,  0.95473,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  // NOLINT
    0.00000,  0.00000,  0.00000,  0.00000,  0.59374,  0.00000,  0.00000,  0.00000,  0.49139,  0.00000,  // NOLINT
    0.91276,  0.96641,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.91797,  // NOLINT
    0.96828,  0.00000,  0.00000,  0.72583,  0.00000,  0.00000,  0.81459,  0.00000,  0.04560,  0.00000   // NOLINT
  };

  CompressedRowSparseMatrix matrix(10, 10, 100);
  int* rows = matrix.mutable_rows();
  int* cols = matrix.mutable_cols();
  double* values = matrix.mutable_values();

  int idx = 0;
  for (int i = 0; i < 10; ++i) {
    rows[i] = idx;
    for (int j = 0; j < 10; ++j) {
      const double v = dense_matrix[i * 10 + j];
      if (fabs(v) > 1e-6) {
        cols[idx] = j;
        values[idx] = v;
        ++idx;
      }
    }
  }
  rows[10] = idx;

  scoped_ptr<CompressedRowSparseMatrix> lmatrix(
      IncompleteLQFactorization(matrix, 10, 0.0, 10, 0.0));

  ConstMatrixRef mref(dense_matrix, 10, 10);

  // Use Cholesky factorization to compute the L matrix.
  Matrix expected_l_matrix  = (mref * mref.transpose()).llt().matrixL();
  Matrix actual_l_matrix;
  lmatrix->ToDenseMatrix(&actual_l_matrix);

  EXPECT_NEAR((expected_l_matrix * expected_l_matrix.transpose() -
               actual_l_matrix * actual_l_matrix.transpose()).norm(),
              0.0,
              1e-10)
      << "expected: \n" << expected_l_matrix
      << "\actual: \n" << actual_l_matrix;
}

TEST(IncompleteLQFactorization, DropEntriesAndAddRow) {
  // Allocate space and then make it a zero sized matrix.
  CompressedRowSparseMatrix matrix(10, 10, 100);
  matrix.set_num_rows(0);

  vector<pair<int, double> > scratch(10);

  Vector dense_vector(10);
  dense_vector(0) = 5;
  dense_vector(1) = 1;
  dense_vector(2) = 2;
  dense_vector(3) = 3;
  dense_vector(4) = 1;
  dense_vector(5) = 4;

  // Add a row with just one entry.
  DropEntriesAndAddRow(dense_vector, 1, 1, 0, &scratch, &matrix);
  EXPECT_EQ(matrix.num_rows(), 1);
  EXPECT_EQ(matrix.num_cols(), 10);
  EXPECT_EQ(matrix.num_nonzeros(), 1);
  EXPECT_EQ(matrix.values()[0], 5.0);
  EXPECT_EQ(matrix.cols()[0], 0);

  // Add a row with six entries
  DropEntriesAndAddRow(dense_vector, 6, 10, 0, &scratch, &matrix);
  EXPECT_EQ(matrix.num_rows(), 2);
  EXPECT_EQ(matrix.num_cols(), 10);
  EXPECT_EQ(matrix.num_nonzeros(), 7);
  for (int idx = matrix.rows()[1]; idx < matrix.rows()[2]; ++idx) {
    EXPECT_EQ(matrix.cols()[idx], idx - matrix.rows()[1]);
    EXPECT_EQ(matrix.values()[idx], dense_vector(idx - matrix.rows()[1]));
  }

  // Add the top 3 entries.
  DropEntriesAndAddRow(dense_vector, 6, 3, 0, &scratch, &matrix);
  EXPECT_EQ(matrix.num_rows(), 3);
  EXPECT_EQ(matrix.num_cols(), 10);
  EXPECT_EQ(matrix.num_nonzeros(), 10);

  EXPECT_EQ(matrix.cols()[matrix.rows()[2]], 0);
  EXPECT_EQ(matrix.cols()[matrix.rows()[2] + 1], 3);
  EXPECT_EQ(matrix.cols()[matrix.rows()[2] + 2], 5);

  EXPECT_EQ(matrix.values()[matrix.rows()[2]], 5);
  EXPECT_EQ(matrix.values()[matrix.rows()[2] + 1], 3);
  EXPECT_EQ(matrix.values()[matrix.rows()[2] + 2], 4);

  // Only keep entries greater than 1.0;
  DropEntriesAndAddRow(dense_vector, 6, 6, 0.2, &scratch, &matrix);
  EXPECT_EQ(matrix.num_rows(), 4);
  EXPECT_EQ(matrix.num_cols(), 10);
  EXPECT_EQ(matrix.num_nonzeros(), 14);

  EXPECT_EQ(matrix.cols()[matrix.rows()[3]], 0);
  EXPECT_EQ(matrix.cols()[matrix.rows()[3] + 1], 2);
  EXPECT_EQ(matrix.cols()[matrix.rows()[3] + 2], 3);
  EXPECT_EQ(matrix.cols()[matrix.rows()[3] + 3], 5);

  EXPECT_EQ(matrix.values()[matrix.rows()[3]], 5);
  EXPECT_EQ(matrix.values()[matrix.rows()[3] + 1], 2);
  EXPECT_EQ(matrix.values()[matrix.rows()[3] + 2], 3);
  EXPECT_EQ(matrix.values()[matrix.rows()[3] + 3], 4);

  // Only keep the top 2 entries greater than 1.0
  DropEntriesAndAddRow(dense_vector, 6, 2, 0.2, &scratch, &matrix);
  EXPECT_EQ(matrix.num_rows(), 5);
  EXPECT_EQ(matrix.num_cols(), 10);
  EXPECT_EQ(matrix.num_nonzeros(), 16);

  EXPECT_EQ(matrix.cols()[matrix.rows()[4]], 0);
  EXPECT_EQ(matrix.cols()[matrix.rows()[4] + 1], 5);

  EXPECT_EQ(matrix.values()[matrix.rows()[4]], 5);
  EXPECT_EQ(matrix.values()[matrix.rows()[4] + 1], 4);
}


}  // namespace internal
}  // namespace ceres
