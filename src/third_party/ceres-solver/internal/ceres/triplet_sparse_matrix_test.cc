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

#include "ceres/triplet_sparse_matrix.h"

#include <memory>

#include "ceres/crs_matrix.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(TripletSparseMatrix, DefaultConstructorReturnsEmptyObject) {
  TripletSparseMatrix m;
  EXPECT_EQ(m.num_rows(), 0);
  EXPECT_EQ(m.num_cols(), 0);
  EXPECT_EQ(m.num_nonzeros(), 0);
  EXPECT_EQ(m.max_num_nonzeros(), 0);
}

TEST(TripletSparseMatrix, SimpleConstructorAndBasicOperations) {
  // Build a matrix
  TripletSparseMatrix m(2, 5, 4);
  EXPECT_EQ(m.num_rows(), 2);
  EXPECT_EQ(m.num_cols(), 5);
  EXPECT_EQ(m.num_nonzeros(), 0);
  EXPECT_EQ(m.max_num_nonzeros(), 4);

  m.mutable_rows()[0] = 0;
  m.mutable_cols()[0] = 1;
  m.mutable_values()[0] = 2.5;

  m.mutable_rows()[1] = 1;
  m.mutable_cols()[1] = 4;
  m.mutable_values()[1] = 5.2;
  m.set_num_nonzeros(2);

  EXPECT_EQ(m.num_nonzeros(), 2);

  ASSERT_TRUE(m.AllTripletsWithinBounds());

  // We should never be able resize and lose data
  EXPECT_DEATH_IF_SUPPORTED(m.Reserve(1), "Reallocation will cause data loss");

  // We should be able to resize while preserving data
  m.Reserve(50);
  EXPECT_EQ(m.max_num_nonzeros(), 50);

  m.Reserve(3);
  EXPECT_EQ(m.max_num_nonzeros(), 50);  // The space is already reserved.

  EXPECT_EQ(m.rows()[0], 0);
  EXPECT_EQ(m.rows()[1], 1);

  EXPECT_EQ(m.cols()[0], 1);
  EXPECT_EQ(m.cols()[1], 4);

  EXPECT_DOUBLE_EQ(m.values()[0], 2.5);
  EXPECT_DOUBLE_EQ(m.values()[1], 5.2);

  // Bounds check should fail
  m.mutable_rows()[0] = 10;
  EXPECT_FALSE(m.AllTripletsWithinBounds());

  m.mutable_rows()[0] = 1;
  m.mutable_cols()[0] = 100;
  EXPECT_FALSE(m.AllTripletsWithinBounds());

  // Remove all data and then resize the data store
  m.SetZero();
  EXPECT_EQ(m.num_nonzeros(), 0);
  m.Reserve(1);
}

TEST(TripletSparseMatrix, CopyConstructor) {
  TripletSparseMatrix orig(2, 5, 4);
  orig.mutable_rows()[0] = 0;
  orig.mutable_cols()[0] = 1;
  orig.mutable_values()[0] = 2.5;

  orig.mutable_rows()[1] = 1;
  orig.mutable_cols()[1] = 4;
  orig.mutable_values()[1] = 5.2;
  orig.set_num_nonzeros(2);

  TripletSparseMatrix cpy(orig);

  EXPECT_EQ(cpy.num_rows(), 2);
  EXPECT_EQ(cpy.num_cols(), 5);
  ASSERT_EQ(cpy.num_nonzeros(), 2);
  EXPECT_EQ(cpy.max_num_nonzeros(), 4);

  EXPECT_EQ(cpy.rows()[0], 0);
  EXPECT_EQ(cpy.rows()[1], 1);

  EXPECT_EQ(cpy.cols()[0], 1);
  EXPECT_EQ(cpy.cols()[1], 4);

  EXPECT_DOUBLE_EQ(cpy.values()[0], 2.5);
  EXPECT_DOUBLE_EQ(cpy.values()[1], 5.2);
}

TEST(TripletSparseMatrix, AssignmentOperator) {
  TripletSparseMatrix orig(2, 5, 4);
  orig.mutable_rows()[0] = 0;
  orig.mutable_cols()[0] = 1;
  orig.mutable_values()[0] = 2.5;

  orig.mutable_rows()[1] = 1;
  orig.mutable_cols()[1] = 4;
  orig.mutable_values()[1] = 5.2;
  orig.set_num_nonzeros(2);

  TripletSparseMatrix cpy(3, 50, 40);
  cpy.mutable_rows()[0] = 0;
  cpy.mutable_cols()[0] = 10;
  cpy.mutable_values()[0] = 10.22;

  cpy.mutable_rows()[1] = 2;
  cpy.mutable_cols()[1] = 23;
  cpy.mutable_values()[1] = 34.45;

  cpy.mutable_rows()[0] = 0;
  cpy.mutable_cols()[0] = 10;
  cpy.mutable_values()[0] = 10.22;

  cpy.mutable_rows()[1] = 0;
  cpy.mutable_cols()[1] = 3;
  cpy.mutable_values()[1] = 4.4;
  cpy.set_num_nonzeros(3);

  cpy = orig;

  EXPECT_EQ(cpy.num_rows(), 2);
  EXPECT_EQ(cpy.num_cols(), 5);
  ASSERT_EQ(cpy.num_nonzeros(), 2);
  EXPECT_EQ(cpy.max_num_nonzeros(), 4);

  EXPECT_EQ(cpy.rows()[0], 0);
  EXPECT_EQ(cpy.rows()[1], 1);

  EXPECT_EQ(cpy.cols()[0], 1);
  EXPECT_EQ(cpy.cols()[1], 4);

  EXPECT_DOUBLE_EQ(cpy.values()[0], 2.5);
  EXPECT_DOUBLE_EQ(cpy.values()[1], 5.2);
}

TEST(TripletSparseMatrix, AssignmentOperatorSelfAssignment) {
  TripletSparseMatrix orig(2, 5, 4);
  orig.mutable_rows()[0] = 0;
  orig.mutable_cols()[0] = 1;
  orig.mutable_values()[0] = 2.5;

  orig.mutable_rows()[1] = 1;
  orig.mutable_cols()[1] = 4;
  orig.mutable_values()[1] = 5.2;
  orig.set_num_nonzeros(2);

  // Who's on earth gonna do this?
  orig = orig;

  EXPECT_EQ(orig.num_rows(), 2);
  EXPECT_EQ(orig.num_cols(), 5);
  ASSERT_EQ(orig.num_nonzeros(), 2);
  EXPECT_EQ(orig.max_num_nonzeros(), 4);

  EXPECT_EQ(orig.rows()[0], 0);
  EXPECT_EQ(orig.rows()[1], 1);

  EXPECT_EQ(orig.cols()[0], 1);
  EXPECT_EQ(orig.cols()[1], 4);

  EXPECT_DOUBLE_EQ(orig.values()[0], 2.5);
  EXPECT_DOUBLE_EQ(orig.values()[1], 5.2);
}

TEST(TripletSparseMatrix, AppendRows) {
  // Build one matrix.
  TripletSparseMatrix m(2, 5, 4);
  m.mutable_rows()[0] = 0;
  m.mutable_cols()[0] = 1;
  m.mutable_values()[0] = 2.5;

  m.mutable_rows()[1] = 1;
  m.mutable_cols()[1] = 4;
  m.mutable_values()[1] = 5.2;
  m.set_num_nonzeros(2);

  // Build another matrix.
  TripletSparseMatrix a(10, 5, 4);
  a.mutable_rows()[0] = 0;
  a.mutable_cols()[0] = 1;
  a.mutable_values()[0] = 3.5;

  a.mutable_rows()[1] = 1;
  a.mutable_cols()[1] = 4;
  a.mutable_values()[1] = 6.2;

  a.mutable_rows()[2] = 9;
  a.mutable_cols()[2] = 5;
  a.mutable_values()[2] = 1;
  a.set_num_nonzeros(3);

  // Glue the second matrix to the bottom of the first.
  m.AppendRows(a);

  EXPECT_EQ(m.num_rows(), 12);
  EXPECT_EQ(m.num_cols(), 5);
  ASSERT_EQ(m.num_nonzeros(), 5);

  EXPECT_EQ(m.values()[0], 2.5);
  EXPECT_EQ(m.values()[1], 5.2);
  EXPECT_EQ(m.values()[2], 3.5);
  EXPECT_EQ(m.values()[3], 6.2);
  EXPECT_EQ(m.values()[4], 1);

  EXPECT_EQ(m.rows()[0], 0);
  EXPECT_EQ(m.rows()[1], 1);
  EXPECT_EQ(m.rows()[2], 2);
  EXPECT_EQ(m.rows()[3], 3);
  EXPECT_EQ(m.rows()[4], 11);

  EXPECT_EQ(m.cols()[0], 1);
  EXPECT_EQ(m.cols()[1], 4);
  EXPECT_EQ(m.cols()[2], 1);
  EXPECT_EQ(m.cols()[3], 4);
  EXPECT_EQ(m.cols()[4], 5);
}

TEST(TripletSparseMatrix, AppendCols) {
  // Build one matrix.
  TripletSparseMatrix m(2, 5, 4);
  m.mutable_rows()[0] = 0;
  m.mutable_cols()[0] = 1;
  m.mutable_values()[0] = 2.5;

  m.mutable_rows()[1] = 1;
  m.mutable_cols()[1] = 4;
  m.mutable_values()[1] = 5.2;
  m.set_num_nonzeros(2);

  // Build another matrix.
  TripletSparseMatrix a(2, 15, 4);
  a.mutable_rows()[0] = 0;
  a.mutable_cols()[0] = 1;
  a.mutable_values()[0] = 3.5;

  a.mutable_rows()[1] = 1;
  a.mutable_cols()[1] = 4;
  a.mutable_values()[1] = 6.2;

  a.mutable_rows()[2] = 0;
  a.mutable_cols()[2] = 10;
  a.mutable_values()[2] = 1;
  a.set_num_nonzeros(3);

  // Glue the second matrix to the left of the first.
  m.AppendCols(a);

  EXPECT_EQ(m.num_rows(), 2);
  EXPECT_EQ(m.num_cols(), 20);
  ASSERT_EQ(m.num_nonzeros(), 5);

  EXPECT_EQ(m.values()[0], 2.5);
  EXPECT_EQ(m.values()[1], 5.2);
  EXPECT_EQ(m.values()[2], 3.5);
  EXPECT_EQ(m.values()[3], 6.2);
  EXPECT_EQ(m.values()[4], 1);

  EXPECT_EQ(m.rows()[0], 0);
  EXPECT_EQ(m.rows()[1], 1);
  EXPECT_EQ(m.rows()[2], 0);
  EXPECT_EQ(m.rows()[3], 1);
  EXPECT_EQ(m.rows()[4], 0);

  EXPECT_EQ(m.cols()[0], 1);
  EXPECT_EQ(m.cols()[1], 4);
  EXPECT_EQ(m.cols()[2], 6);
  EXPECT_EQ(m.cols()[3], 9);
  EXPECT_EQ(m.cols()[4], 15);
}

TEST(TripletSparseMatrix, CreateDiagonalMatrix) {
  std::unique_ptr<double[]> values(new double[10]);
  for (int i = 0; i < 10; ++i) values[i] = i;

  std::unique_ptr<TripletSparseMatrix> m(
      TripletSparseMatrix::CreateSparseDiagonalMatrix(values.get(), 10));
  EXPECT_EQ(m->num_rows(), 10);
  EXPECT_EQ(m->num_cols(), 10);
  ASSERT_EQ(m->num_nonzeros(), 10);
  for (int i = 0; i < 10; ++i) {
    EXPECT_EQ(m->rows()[i], i);
    EXPECT_EQ(m->cols()[i], i);
    EXPECT_EQ(m->values()[i], i);
  }
}

TEST(TripletSparseMatrix, Resize) {
  TripletSparseMatrix m(10, 20, 200);
  int nnz = 0;
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 20; ++j) {
      m.mutable_rows()[nnz] = i;
      m.mutable_cols()[nnz] = j;
      m.mutable_values()[nnz++] = i + j;
    }
  }
  m.set_num_nonzeros(nnz);
  m.Resize(5, 6);
  EXPECT_EQ(m.num_rows(), 5);
  EXPECT_EQ(m.num_cols(), 6);
  ASSERT_EQ(m.num_nonzeros(), 30);
  for (int i = 0; i < 30; ++i) {
    EXPECT_EQ(m.values()[i], m.rows()[i] + m.cols()[i]);
  }
}

TEST(TripletSparseMatrix, ToCRSMatrix) {
  // Test matrix:
  // [1, 2, 0, 5, 6, 0,
  //  3, 4, 0, 7, 8, 0,
  //  0, 0, 9, 0, 0, 0]
  TripletSparseMatrix m(3,
                        6,
                        {0, 0, 0, 0, 1, 1, 1, 1, 2},
                        {0, 1, 3, 4, 0, 1, 3, 4, 2},
                        {1, 2, 3, 4, 5, 6, 7, 8, 9});
  CRSMatrix m_crs;
  m.ToCRSMatrix(&m_crs);
  EXPECT_EQ(m_crs.num_rows, 3);
  EXPECT_EQ(m_crs.num_cols, 6);

  EXPECT_EQ(m_crs.rows.size(), 4);
  EXPECT_EQ(m_crs.rows[0], 0);
  EXPECT_EQ(m_crs.rows[1], 4);
  EXPECT_EQ(m_crs.rows[2], 8);
  EXPECT_EQ(m_crs.rows[3], 9);

  EXPECT_EQ(m_crs.cols.size(), 9);
  EXPECT_EQ(m_crs.cols[0], 0);
  EXPECT_EQ(m_crs.cols[1], 1);
  EXPECT_EQ(m_crs.cols[2], 3);
  EXPECT_EQ(m_crs.cols[3], 4);
  EXPECT_EQ(m_crs.cols[4], 0);
  EXPECT_EQ(m_crs.cols[5], 1);
  EXPECT_EQ(m_crs.cols[6], 3);
  EXPECT_EQ(m_crs.cols[7], 4);
  EXPECT_EQ(m_crs.cols[8], 2);

  EXPECT_EQ(m_crs.values.size(), 9);
  for (int i = 0; i < 9; ++i) {
    EXPECT_EQ(m_crs.values[i], i + 1);
  }
}

}  // namespace ceres::internal
