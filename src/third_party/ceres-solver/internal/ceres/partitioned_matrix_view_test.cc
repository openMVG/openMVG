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

#include "ceres/partitioned_matrix_view.h"

#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "ceres/block_structure.h"
#include "ceres/casts.h"
#include "ceres/internal/eigen.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

const double kEpsilon = 1e-14;

// Param = <problem_id, num_threads>
using Param = ::testing::tuple<int, int>;

static std::string ParamInfoToString(testing::TestParamInfo<Param> info) {
  Param param = info.param;
  std::stringstream ss;
  ss << ::testing::get<0>(param) << "_" << ::testing::get<1>(param);
  return ss.str();
}

class PartitionedMatrixViewTest : public ::testing::TestWithParam<Param> {
 protected:
  void SetUp() final {
    const int problem_id = ::testing::get<0>(GetParam());
    const int num_threads = ::testing::get<1>(GetParam());
    auto problem = CreateLinearLeastSquaresProblemFromId(problem_id);
    CHECK(problem != nullptr);
    A_ = std::move(problem->A);
    auto block_sparse = down_cast<BlockSparseMatrix*>(A_.get());

    options_.num_threads = num_threads;
    options_.context = &context_;
    options_.elimination_groups.push_back(problem->num_eliminate_blocks);
    pmv_ = PartitionedMatrixViewBase::Create(options_, *block_sparse);

    LinearSolver::Options options_single_threaded = options_;
    options_single_threaded.num_threads = 1;
    pmv_single_threaded_ =
        PartitionedMatrixViewBase::Create(options_, *block_sparse);

    EXPECT_EQ(pmv_->num_col_blocks_e(), problem->num_eliminate_blocks);
    EXPECT_EQ(pmv_->num_col_blocks_f(),
              block_sparse->block_structure()->cols.size() -
                  problem->num_eliminate_blocks);
    EXPECT_EQ(pmv_->num_cols(), A_->num_cols());
    EXPECT_EQ(pmv_->num_rows(), A_->num_rows());
  }

  double RandDouble() { return distribution_(prng_); }

  LinearSolver::Options options_;
  ContextImpl context_;
  std::unique_ptr<LinearLeastSquaresProblem> problem_;
  std::unique_ptr<SparseMatrix> A_;
  std::unique_ptr<PartitionedMatrixViewBase> pmv_;
  std::unique_ptr<PartitionedMatrixViewBase> pmv_single_threaded_;
  std::mt19937 prng_;
  std::uniform_real_distribution<double> distribution_ =
      std::uniform_real_distribution<double>(0.0, 1.0);
};

TEST_P(PartitionedMatrixViewTest, RightMultiplyAndAccumulateE) {
  Vector x1(pmv_->num_cols_e());
  Vector x2(pmv_->num_cols());
  x2.setZero();

  for (int i = 0; i < pmv_->num_cols_e(); ++i) {
    x1(i) = x2(i) = RandDouble();
  }

  Vector expected = Vector::Zero(pmv_->num_rows());
  A_->RightMultiplyAndAccumulate(x2.data(), expected.data());

  Vector actual = Vector::Zero(pmv_->num_rows());
  pmv_->RightMultiplyAndAccumulateE(x1.data(), actual.data());

  for (int i = 0; i < pmv_->num_rows(); ++i) {
    EXPECT_NEAR(actual(i), expected(i), kEpsilon);
  }
}

TEST_P(PartitionedMatrixViewTest, RightMultiplyAndAccumulateF) {
  Vector x1(pmv_->num_cols_f());
  Vector x2(pmv_->num_cols());
  x2.setZero();

  for (int i = 0; i < pmv_->num_cols_f(); ++i) {
    x1(i) = x2(i + pmv_->num_cols_e()) = RandDouble();
  }

  Vector actual = Vector::Zero(pmv_->num_rows());
  pmv_->RightMultiplyAndAccumulateF(x1.data(), actual.data());

  Vector expected = Vector::Zero(pmv_->num_rows());
  A_->RightMultiplyAndAccumulate(x2.data(), expected.data());

  for (int i = 0; i < pmv_->num_rows(); ++i) {
    EXPECT_NEAR(actual(i), expected(i), kEpsilon);
  }
}

TEST_P(PartitionedMatrixViewTest, LeftMultiplyAndAccumulate) {
  Vector x = Vector::Zero(pmv_->num_rows());
  for (int i = 0; i < pmv_->num_rows(); ++i) {
    x(i) = RandDouble();
  }
  Vector x_pre = x;

  Vector expected = Vector::Zero(pmv_->num_cols());
  Vector e_actual = Vector::Zero(pmv_->num_cols_e());
  Vector f_actual = Vector::Zero(pmv_->num_cols_f());

  A_->LeftMultiplyAndAccumulate(x.data(), expected.data());
  pmv_->LeftMultiplyAndAccumulateE(x.data(), e_actual.data());
  pmv_->LeftMultiplyAndAccumulateF(x.data(), f_actual.data());

  for (int i = 0; i < pmv_->num_cols(); ++i) {
    EXPECT_NEAR(expected(i),
                (i < pmv_->num_cols_e()) ? e_actual(i)
                                         : f_actual(i - pmv_->num_cols_e()),
                kEpsilon);
  }
}

TEST_P(PartitionedMatrixViewTest, BlockDiagonalFtF) {
  std::unique_ptr<BlockSparseMatrix> block_diagonal_ff(
      pmv_->CreateBlockDiagonalFtF());
  const auto bs_diagonal = block_diagonal_ff->block_structure();
  const int num_rows = pmv_->num_rows();
  const int num_cols_f = pmv_->num_cols_f();
  const int num_cols_e = pmv_->num_cols_e();
  const int num_col_blocks_f = pmv_->num_col_blocks_f();
  const int num_col_blocks_e = pmv_->num_col_blocks_e();

  CHECK_EQ(block_diagonal_ff->num_rows(), num_cols_f);
  CHECK_EQ(block_diagonal_ff->num_cols(), num_cols_f);

  EXPECT_EQ(bs_diagonal->cols.size(), num_col_blocks_f);
  EXPECT_EQ(bs_diagonal->rows.size(), num_col_blocks_f);

  Matrix EF;
  A_->ToDenseMatrix(&EF);
  const auto F = EF.topRightCorner(num_rows, num_cols_f);

  Matrix expected_FtF = F.transpose() * F;
  Matrix actual_FtF;
  block_diagonal_ff->ToDenseMatrix(&actual_FtF);

  // FtF might be not block-diagonal
  auto bs = down_cast<BlockSparseMatrix*>(A_.get())->block_structure();
  for (int i = 0; i < num_col_blocks_f; ++i) {
    const auto col_block_f = bs->cols[num_col_blocks_e + i];
    const int block_size = col_block_f.size;
    const int block_pos = col_block_f.position - num_cols_e;
    const auto cell_expected =
        expected_FtF.block(block_pos, block_pos, block_size, block_size);
    auto cell_actual =
        actual_FtF.block(block_pos, block_pos, block_size, block_size);
    cell_actual -= cell_expected;
    EXPECT_NEAR(cell_actual.norm(), 0., kEpsilon);
  }
  // There should be nothing remaining outside block-diagonal
  EXPECT_NEAR(actual_FtF.norm(), 0., kEpsilon);
}

TEST_P(PartitionedMatrixViewTest, BlockDiagonalEtE) {
  std::unique_ptr<BlockSparseMatrix> block_diagonal_ee(
      pmv_->CreateBlockDiagonalEtE());
  const CompressedRowBlockStructure* bs = block_diagonal_ee->block_structure();
  const int num_rows = pmv_->num_rows();
  const int num_cols_e = pmv_->num_cols_e();
  const int num_col_blocks_e = pmv_->num_col_blocks_e();

  CHECK_EQ(block_diagonal_ee->num_rows(), num_cols_e);
  CHECK_EQ(block_diagonal_ee->num_cols(), num_cols_e);

  EXPECT_EQ(bs->cols.size(), num_col_blocks_e);
  EXPECT_EQ(bs->rows.size(), num_col_blocks_e);

  Matrix EF;
  A_->ToDenseMatrix(&EF);
  const auto E = EF.topLeftCorner(num_rows, num_cols_e);

  Matrix expected_EtE = E.transpose() * E;
  Matrix actual_EtE;
  block_diagonal_ee->ToDenseMatrix(&actual_EtE);

  EXPECT_NEAR((expected_EtE - actual_EtE).norm(), 0., kEpsilon);
}

TEST_P(PartitionedMatrixViewTest, UpdateBlockDiagonalEtE) {
  std::unique_ptr<BlockSparseMatrix> block_diagonal_ete(
      pmv_->CreateBlockDiagonalEtE());
  const int num_cols = pmv_->num_cols_e();

  Matrix multi_threaded(num_cols, num_cols);
  pmv_->UpdateBlockDiagonalEtE(block_diagonal_ete.get());
  block_diagonal_ete->ToDenseMatrix(&multi_threaded);

  Matrix single_threaded(num_cols, num_cols);
  pmv_single_threaded_->UpdateBlockDiagonalEtE(block_diagonal_ete.get());
  block_diagonal_ete->ToDenseMatrix(&single_threaded);

  EXPECT_NEAR((multi_threaded - single_threaded).norm(), 0., kEpsilon);
}

TEST_P(PartitionedMatrixViewTest, UpdateBlockDiagonalFtF) {
  std::unique_ptr<BlockSparseMatrix> block_diagonal_ftf(
      pmv_->CreateBlockDiagonalFtF());
  const int num_cols = pmv_->num_cols_f();

  Matrix multi_threaded(num_cols, num_cols);
  pmv_->UpdateBlockDiagonalFtF(block_diagonal_ftf.get());
  block_diagonal_ftf->ToDenseMatrix(&multi_threaded);

  Matrix single_threaded(num_cols, num_cols);
  pmv_single_threaded_->UpdateBlockDiagonalFtF(block_diagonal_ftf.get());
  block_diagonal_ftf->ToDenseMatrix(&single_threaded);

  EXPECT_NEAR((multi_threaded - single_threaded).norm(), 0., kEpsilon);
}

INSTANTIATE_TEST_SUITE_P(
    ParallelProducts,
    PartitionedMatrixViewTest,
    ::testing::Combine(::testing::Values(2, 4, 6),
                       ::testing::Values(1, 2, 3, 4, 5, 6, 7, 8)),
    ParamInfoToString);

}  // namespace internal
}  // namespace ceres
