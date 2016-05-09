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
// Author: richie.stebbing@gmail.com (Richard Stebbing)

#include "ceres/dynamic_compressed_row_sparse_matrix.h"

#include "ceres/casts.h"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/triplet_sparse_matrix.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

using std::copy;
using std::vector;

class DynamicCompressedRowSparseMatrixTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    num_rows = 7;
    num_cols = 4;

    // The number of additional elements reserved when `Finalize` is called
    // should have no effect on the number of rows, columns or nonzeros.
    // Set this to some nonzero value to be sure.
    num_additional_elements = 13;

    expected_num_nonzeros = num_rows * num_cols - std::min(num_rows, num_cols);

    InitialiseDenseReference();
    InitialiseSparseMatrixReferences();

    dcrsm.reset(new DynamicCompressedRowSparseMatrix(num_rows,
                                                     num_cols,
                                                     0));
  }

  void Finalize() {
    dcrsm->Finalize(num_additional_elements);
  }

  void InitialiseDenseReference() {
    dense.resize(num_rows, num_cols);
    dense.setZero();
    int num_nonzeros = 0;
    for (int i = 0; i < (num_rows * num_cols); ++i) {
      const int r = i / num_cols, c = i % num_cols;
      if (r != c) {
        dense(r, c) = i + 1;
        ++num_nonzeros;
      }
    }
    ASSERT_EQ(num_nonzeros, expected_num_nonzeros);
  }

  void InitialiseSparseMatrixReferences() {
    vector<int> rows, cols;
    vector<double> values;
    for (int i = 0; i < (num_rows * num_cols); ++i) {
      const int r = i / num_cols, c = i % num_cols;
      if (r != c) {
        rows.push_back(r);
        cols.push_back(c);
        values.push_back(i + 1);
      }
    }
    ASSERT_EQ(values.size(), expected_num_nonzeros);

    tsm.reset(new TripletSparseMatrix(num_rows,
                                      num_cols,
                                      expected_num_nonzeros));
    copy(rows.begin(), rows.end(), tsm->mutable_rows());
    copy(cols.begin(), cols.end(), tsm->mutable_cols());
    copy(values.begin(), values.end(), tsm->mutable_values());
    tsm->set_num_nonzeros(values.size());

    Matrix dense_from_tsm;
    tsm->ToDenseMatrix(&dense_from_tsm);
    ASSERT_TRUE((dense.array() == dense_from_tsm.array()).all());

    crsm.reset(new CompressedRowSparseMatrix(*tsm));
    Matrix dense_from_crsm;
    crsm->ToDenseMatrix(&dense_from_crsm);
    ASSERT_TRUE((dense.array() == dense_from_crsm.array()).all());
  }

  void InsertNonZeroEntriesFromDenseReference() {
    for (int r = 0; r < num_rows; ++r) {
      for (int c = 0; c < num_cols; ++c) {
        const double& v = dense(r, c);
        if (v != 0.0) {
          dcrsm->InsertEntry(r, c, v);
        }
      }
    }
  }

  void ExpectEmpty() {
    EXPECT_EQ(dcrsm->num_rows(), num_rows);
    EXPECT_EQ(dcrsm->num_cols(), num_cols);
    EXPECT_EQ(dcrsm->num_nonzeros(), 0);

    Matrix dense_from_dcrsm;
    dcrsm->ToDenseMatrix(&dense_from_dcrsm);
    EXPECT_EQ(dense_from_dcrsm.rows(), num_rows);
    EXPECT_EQ(dense_from_dcrsm.cols(), num_cols);
    EXPECT_TRUE((dense_from_dcrsm.array() == 0.0).all());
  }

  void ExpectEqualToDenseReference() {
    Matrix dense_from_dcrsm;
    dcrsm->ToDenseMatrix(&dense_from_dcrsm);
    EXPECT_TRUE((dense.array() == dense_from_dcrsm.array()).all());
  }

  void ExpectEqualToCompressedRowSparseMatrixReference() {
    typedef Eigen::Map<const Eigen::VectorXi> ConstIntVectorRef;

    ConstIntVectorRef crsm_rows(crsm->rows(), crsm->num_rows() + 1);
    ConstIntVectorRef dcrsm_rows(dcrsm->rows(), dcrsm->num_rows() + 1);
    EXPECT_TRUE((crsm_rows.array() == dcrsm_rows.array()).all());

    ConstIntVectorRef crsm_cols(crsm->cols(), crsm->num_nonzeros());
    ConstIntVectorRef dcrsm_cols(dcrsm->cols(), dcrsm->num_nonzeros());
    EXPECT_TRUE((crsm_cols.array() == dcrsm_cols.array()).all());

    ConstVectorRef crsm_values(crsm->values(), crsm->num_nonzeros());
    ConstVectorRef dcrsm_values(dcrsm->values(), dcrsm->num_nonzeros());
    EXPECT_TRUE((crsm_values.array() == dcrsm_values.array()).all());
  }

  int num_rows;
  int num_cols;

  int num_additional_elements;

  int expected_num_nonzeros;

  Matrix dense;
  scoped_ptr<TripletSparseMatrix> tsm;
  scoped_ptr<CompressedRowSparseMatrix> crsm;

  scoped_ptr<DynamicCompressedRowSparseMatrix> dcrsm;
};

TEST_F(DynamicCompressedRowSparseMatrixTest, Initialization) {
  ExpectEmpty();

  Finalize();
  ExpectEmpty();
}

TEST_F(DynamicCompressedRowSparseMatrixTest, InsertEntryAndFinalize) {
  InsertNonZeroEntriesFromDenseReference();
  ExpectEmpty();

  Finalize();
  ExpectEqualToDenseReference();
  ExpectEqualToCompressedRowSparseMatrixReference();
}

TEST_F(DynamicCompressedRowSparseMatrixTest, ClearRows) {
  InsertNonZeroEntriesFromDenseReference();
  Finalize();
  ExpectEqualToDenseReference();
  ExpectEqualToCompressedRowSparseMatrixReference();

  dcrsm->ClearRows(0, 0);
  Finalize();
  ExpectEqualToDenseReference();
  ExpectEqualToCompressedRowSparseMatrixReference();

  dcrsm->ClearRows(0, num_rows);
  ExpectEqualToCompressedRowSparseMatrixReference();

  Finalize();
  ExpectEmpty();

  InsertNonZeroEntriesFromDenseReference();
  dcrsm->ClearRows(1, 2);
  Finalize();
  dense.block(1, 0, 2, num_cols).setZero();
  ExpectEqualToDenseReference();

  InitialiseDenseReference();
}

}  // namespace internal
}  // namespace ceres
