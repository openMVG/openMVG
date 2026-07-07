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

#include "ceres/subset_preconditioner.h"

#include <memory>
#include <random>

#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "ceres/block_sparse_matrix.h"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/inner_product_computer.h"
#include "ceres/internal/config.h"
#include "ceres/internal/eigen.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

namespace {

// TODO(sameeragarwal): Refactor the following two functions out of
// here and sparse_cholesky_test.cc into a more suitable place.
template <int UpLoType>
bool SolveLinearSystemUsingEigen(const Matrix& lhs,
                                 const Vector rhs,
                                 Vector* solution) {
  Eigen::LLT<Matrix, UpLoType> llt = lhs.selfadjointView<UpLoType>().llt();
  if (llt.info() != Eigen::Success) {
    return false;
  }
  *solution = llt.solve(rhs);
  return (llt.info() == Eigen::Success);
}

// Use Eigen's Dense Cholesky solver to compute the solution to a
// sparse linear system.
bool ComputeExpectedSolution(const CompressedRowSparseMatrix& lhs,
                             const Vector& rhs,
                             Vector* solution) {
  Matrix dense_triangular_lhs;
  lhs.ToDenseMatrix(&dense_triangular_lhs);
  if (lhs.storage_type() ==
      CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR) {
    Matrix full_lhs = dense_triangular_lhs.selfadjointView<Eigen::Upper>();
    return SolveLinearSystemUsingEigen<Eigen::Upper>(full_lhs, rhs, solution);
  }
  return SolveLinearSystemUsingEigen<Eigen::Lower>(
      dense_triangular_lhs, rhs, solution);
}

using Param = ::testing::tuple<SparseLinearAlgebraLibraryType, bool>;

std::string ParamInfoToString(testing::TestParamInfo<Param> info) {
  Param param = info.param;
  std::stringstream ss;
  ss << SparseLinearAlgebraLibraryTypeToString(::testing::get<0>(param)) << "_"
     << (::testing::get<1>(param) ? "Diagonal" : "NoDiagonal");
  return ss.str();
}

}  // namespace

class SubsetPreconditionerTest : public ::testing::TestWithParam<Param> {
 protected:
  void SetUp() final {
    BlockSparseMatrix::RandomMatrixOptions options;
    options.num_col_blocks = 4;
    options.min_col_block_size = 1;
    options.max_col_block_size = 4;
    options.num_row_blocks = 8;
    options.min_row_block_size = 1;
    options.max_row_block_size = 4;
    options.block_density = 0.9;

    m_ = BlockSparseMatrix::CreateRandomMatrix(options, prng_);
    start_row_block_ = m_->block_structure()->rows.size();

    // Ensure that the bottom part of the matrix has the same column
    // block structure.
    options.col_blocks = m_->block_structure()->cols;
    b_ = BlockSparseMatrix::CreateRandomMatrix(options, prng_);
    m_->AppendRows(*b_);

    // Create a Identity block diagonal matrix with the same column
    // block structure.
    diagonal_ = Vector::Ones(m_->num_cols());
    block_diagonal_ = BlockSparseMatrix::CreateDiagonalMatrix(
        diagonal_.data(), b_->block_structure()->cols);

    // Unconditionally add the block diagonal to the matrix b_,
    // because either it is either part of b_ to make it full rank, or
    // we pass the same diagonal matrix later as the parameter D. In
    // either case the preconditioner matrix is b_' b + D'D.
    b_->AppendRows(*block_diagonal_);
    inner_product_computer_ = InnerProductComputer::Create(
        *b_, CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
    inner_product_computer_->Compute();
  }

  std::unique_ptr<BlockSparseMatrix> m_;
  std::unique_ptr<BlockSparseMatrix> b_;
  std::unique_ptr<BlockSparseMatrix> block_diagonal_;
  std::unique_ptr<InnerProductComputer> inner_product_computer_;
  std::unique_ptr<Preconditioner> preconditioner_;
  Vector diagonal_;
  int start_row_block_;
  std::mt19937 prng_;
};

TEST_P(SubsetPreconditionerTest, foo) {
  Param param = GetParam();
  Preconditioner::Options options;
  options.subset_preconditioner_start_row_block = start_row_block_;
  options.sparse_linear_algebra_library_type = ::testing::get<0>(param);
  preconditioner_ = std::make_unique<SubsetPreconditioner>(options, *m_);

  const bool with_diagonal = ::testing::get<1>(param);
  if (!with_diagonal) {
    m_->AppendRows(*block_diagonal_);
  }

  EXPECT_TRUE(
      preconditioner_->Update(*m_, with_diagonal ? diagonal_.data() : nullptr));

  // Repeatedly apply the preconditioner to random vectors and check
  // that the preconditioned value is the same as one obtained by
  // solving the linear system directly.
  for (int i = 0; i < 5; ++i) {
    CompressedRowSparseMatrix* lhs = inner_product_computer_->mutable_result();
    Vector rhs = Vector::Random(lhs->num_rows());
    Vector expected(lhs->num_rows());
    EXPECT_TRUE(ComputeExpectedSolution(*lhs, rhs, &expected));

    Vector actual(lhs->num_rows());
    preconditioner_->RightMultiplyAndAccumulate(rhs.data(), actual.data());

    Matrix eigen_lhs;
    lhs->ToDenseMatrix(&eigen_lhs);
    EXPECT_NEAR((actual - expected).norm() / actual.norm(),
                0.0,
                std::numeric_limits<double>::epsilon() * 10)
        << "\n"
        << eigen_lhs << "\n"
        << expected.transpose() << "\n"
        << actual.transpose();
  }
}

#ifndef CERES_NO_SUITESPARSE
INSTANTIATE_TEST_SUITE_P(SubsetPreconditionerWithSuiteSparse,
                         SubsetPreconditionerTest,
                         ::testing::Combine(::testing::Values(SUITE_SPARSE),
                                            ::testing::Values(true, false)),
                         ParamInfoToString);
#endif

#ifndef CERES_NO_ACCELERATE_SPARSE
INSTANTIATE_TEST_SUITE_P(
    SubsetPreconditionerWithAccelerateSparse,
    SubsetPreconditionerTest,
    ::testing::Combine(::testing::Values(ACCELERATE_SPARSE),
                       ::testing::Values(true, false)),
    ParamInfoToString);
#endif

#ifdef CERES_USE_EIGEN_SPARSE
INSTANTIATE_TEST_SUITE_P(SubsetPreconditionerWithEigenSparse,
                         SubsetPreconditionerTest,
                         ::testing::Combine(::testing::Values(EIGEN_SPARSE),
                                            ::testing::Values(true, false)),
                         ParamInfoToString);
#endif

}  // namespace ceres::internal
