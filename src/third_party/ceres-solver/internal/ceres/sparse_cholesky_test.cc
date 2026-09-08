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

#include "ceres/sparse_cholesky.h"

#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "ceres/block_sparse_matrix.h"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/inner_product_computer.h"
#include "ceres/internal/config.h"
#include "ceres/internal/eigen.h"
#include "ceres/iterative_refiner.h"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

namespace {

std::unique_ptr<BlockSparseMatrix> CreateRandomFullRankMatrix(
    const int num_col_blocks,
    const int min_col_block_size,
    const int max_col_block_size,
    const double block_density,
    std::mt19937& prng) {
  // Create a random matrix
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_col_blocks = num_col_blocks;
  options.min_col_block_size = min_col_block_size;
  options.max_col_block_size = max_col_block_size;

  options.num_row_blocks = 2 * num_col_blocks;
  options.min_row_block_size = 1;
  options.max_row_block_size = max_col_block_size;
  options.block_density = block_density;
  auto random_matrix = BlockSparseMatrix::CreateRandomMatrix(options, prng);

  // Add a diagonal block sparse matrix to make it full rank.
  Vector diagonal = Vector::Ones(random_matrix->num_cols());
  auto block_diagonal = BlockSparseMatrix::CreateDiagonalMatrix(
      diagonal.data(), random_matrix->block_structure()->cols);
  random_matrix->AppendRows(*block_diagonal);
  return random_matrix;
}

bool ComputeExpectedSolution(const CompressedRowSparseMatrix& lhs,
                             const Vector& rhs,
                             Vector* solution) {
  Matrix eigen_lhs;
  lhs.ToDenseMatrix(&eigen_lhs);
  if (lhs.storage_type() ==
      CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR) {
    Matrix full_lhs = eigen_lhs.selfadjointView<Eigen::Upper>();
    Eigen::LLT<Matrix, Eigen::Upper> llt =
        eigen_lhs.selfadjointView<Eigen::Upper>().llt();
    if (llt.info() != Eigen::Success) {
      return false;
    }
    *solution = llt.solve(rhs);
    return (llt.info() == Eigen::Success);
  }

  Matrix full_lhs = eigen_lhs.selfadjointView<Eigen::Lower>();
  Eigen::LLT<Matrix, Eigen::Lower> llt =
      eigen_lhs.selfadjointView<Eigen::Lower>().llt();
  if (llt.info() != Eigen::Success) {
    return false;
  }
  *solution = llt.solve(rhs);
  return (llt.info() == Eigen::Success);
}

void SparseCholeskySolverUnitTest(
    const SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type,
    const OrderingType ordering_type,
    const bool use_block_structure,
    const int num_blocks,
    const int min_block_size,
    const int max_block_size,
    const double block_density,
    std::mt19937& prng) {
  LinearSolver::Options sparse_cholesky_options;
  sparse_cholesky_options.sparse_linear_algebra_library_type =
      sparse_linear_algebra_library_type;
  sparse_cholesky_options.ordering_type = ordering_type;
  auto sparse_cholesky = SparseCholesky::Create(sparse_cholesky_options);
  const CompressedRowSparseMatrix::StorageType storage_type =
      sparse_cholesky->StorageType();

  auto m = CreateRandomFullRankMatrix(
      num_blocks, min_block_size, max_block_size, block_density, prng);
  auto inner_product_computer = InnerProductComputer::Create(*m, storage_type);
  inner_product_computer->Compute();
  CompressedRowSparseMatrix* lhs = inner_product_computer->mutable_result();

  if (!use_block_structure) {
    lhs->mutable_row_blocks()->clear();
    lhs->mutable_col_blocks()->clear();
  }

  Vector rhs = Vector::Random(lhs->num_rows());
  Vector expected(lhs->num_rows());
  Vector actual(lhs->num_rows());

  EXPECT_TRUE(ComputeExpectedSolution(*lhs, rhs, &expected));
  std::string message;
  EXPECT_EQ(
      sparse_cholesky->FactorAndSolve(lhs, rhs.data(), actual.data(), &message),
      LinearSolverTerminationType::SUCCESS);
  Matrix eigen_lhs;
  lhs->ToDenseMatrix(&eigen_lhs);
  EXPECT_NEAR((actual - expected).norm() / actual.norm(),
              0.0,
              std::numeric_limits<double>::epsilon() * 20)
      << "\n"
      << eigen_lhs;
}

using Param =
    ::testing::tuple<SparseLinearAlgebraLibraryType, OrderingType, bool>;

std::string ParamInfoToString(testing::TestParamInfo<Param> info) {
  Param param = info.param;
  std::stringstream ss;
  ss << SparseLinearAlgebraLibraryTypeToString(::testing::get<0>(param)) << "_"
     << ::testing::get<1>(param) << "_"
     << (::testing::get<2>(param) ? "UseBlockStructure" : "NoBlockStructure");
  return ss.str();
}

}  // namespace

class SparseCholeskyTest : public ::testing::TestWithParam<Param> {};

TEST_P(SparseCholeskyTest, FactorAndSolve) {
  constexpr int kMinNumBlocks = 1;
  constexpr int kMaxNumBlocks = 10;
  constexpr int kNumTrials = 10;
  constexpr int kMinBlockSize = 1;
  constexpr int kMaxBlockSize = 5;

  Param param = GetParam();

  std::mt19937 prng;
  std::uniform_real_distribution<double> distribution(0.1, 1.0);

  for (int num_blocks = kMinNumBlocks; num_blocks < kMaxNumBlocks;
       ++num_blocks) {
    for (int trial = 0; trial < kNumTrials; ++trial) {
      const double block_density = distribution(prng);
      SparseCholeskySolverUnitTest(::testing::get<0>(param),
                                   ::testing::get<1>(param),
                                   ::testing::get<2>(param),
                                   num_blocks,
                                   kMinBlockSize,
                                   kMaxBlockSize,
                                   block_density,
                                   prng);
    }
  }
}

namespace {

#ifndef CERES_NO_SUITESPARSE
INSTANTIATE_TEST_SUITE_P(
    SuiteSparseCholesky,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(SUITE_SPARSE),
                       ::testing::Values(OrderingType::AMD,
                                         OrderingType::NATURAL),
                       ::testing::Values(true, false)),
    ParamInfoToString);
#endif

#if !defined(CERES_NO_SUITESPARSE) && !defined(CERES_NO_CHOLMOD_PARTITION)
INSTANTIATE_TEST_SUITE_P(
    SuiteSparseCholeskyMETIS,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(SUITE_SPARSE),
                       ::testing::Values(OrderingType::NESDIS),
                       ::testing::Values(true, false)),
    ParamInfoToString);
#endif  // !defined(CERES_NO_SUITESPARSE) &&
        // !defined(CERES_NO_CHOLMOD_PARTITION)

#ifndef CERES_NO_ACCELERATE_SPARSE
INSTANTIATE_TEST_SUITE_P(
    AccelerateSparseCholesky,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(ACCELERATE_SPARSE),
                       ::testing::Values(OrderingType::AMD,
                                         OrderingType::NESDIS,
                                         OrderingType::NATURAL),
                       ::testing::Values(true, false)),
    ParamInfoToString);

INSTANTIATE_TEST_SUITE_P(
    AccelerateSparseCholeskySingle,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(ACCELERATE_SPARSE),
                       ::testing::Values(OrderingType::AMD,
                                         OrderingType::NESDIS,
                                         OrderingType::NATURAL),
                       ::testing::Values(true, false)),
    ParamInfoToString);
#endif

#ifdef CERES_USE_EIGEN_SPARSE
INSTANTIATE_TEST_SUITE_P(
    EigenSparseCholesky,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(EIGEN_SPARSE),
                       ::testing::Values(OrderingType::AMD,
                                         OrderingType::NATURAL),
                       ::testing::Values(true, false)),
    ParamInfoToString);

INSTANTIATE_TEST_SUITE_P(
    EigenSparseCholeskySingle,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(EIGEN_SPARSE),
                       ::testing::Values(OrderingType::AMD,
                                         OrderingType::NATURAL),
                       ::testing::Values(true, false)),
    ParamInfoToString);
#endif  // CERES_USE_EIGEN_SPARSE

#if defined(CERES_USE_EIGEN_SPARSE) && !defined(CERES_NO_EIGEN_METIS)
INSTANTIATE_TEST_SUITE_P(
    EigenSparseCholeskyMETIS,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(EIGEN_SPARSE),
                       ::testing::Values(OrderingType::NESDIS),
                       ::testing::Values(true, false)),
    ParamInfoToString);

INSTANTIATE_TEST_SUITE_P(
    EigenSparseCholeskySingleMETIS,
    SparseCholeskyTest,
    ::testing::Combine(::testing::Values(EIGEN_SPARSE),
                       ::testing::Values(OrderingType::NESDIS),
                       ::testing::Values(true, false)),
    ParamInfoToString);
#endif  // defined(CERES_USE_EIGEN_SPARSE) && !defined(CERES_NO_EIGEN_METIS)

class MockSparseCholesky : public SparseCholesky {
 public:
  MOCK_CONST_METHOD0(StorageType, CompressedRowSparseMatrix::StorageType());
  MOCK_METHOD2(Factorize,
               LinearSolverTerminationType(CompressedRowSparseMatrix* lhs,
                                           std::string* message));
  MOCK_METHOD3(Solve,
               LinearSolverTerminationType(const double* rhs,
                                           double* solution,
                                           std::string* message));
};

class MockSparseIterativeRefiner : public SparseIterativeRefiner {
 public:
  MockSparseIterativeRefiner() : SparseIterativeRefiner(1) {}
  MOCK_METHOD4(Refine,
               void(const SparseMatrix& lhs,
                    const double* rhs,
                    SparseCholesky* sparse_cholesky,
                    double* solution));
};

using testing::_;
using testing::Return;

TEST(RefinedSparseCholesky, StorageType) {
  auto sparse_cholesky = std::make_unique<MockSparseCholesky>();
  auto iterative_refiner = std::make_unique<MockSparseIterativeRefiner>();
  EXPECT_CALL(*sparse_cholesky, StorageType())
      .Times(1)
      .WillRepeatedly(
          Return(CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR));
  EXPECT_CALL(*iterative_refiner, Refine(_, _, _, _)).Times(0);
  RefinedSparseCholesky refined_sparse_cholesky(std::move(sparse_cholesky),
                                                std::move(iterative_refiner));
  EXPECT_EQ(refined_sparse_cholesky.StorageType(),
            CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
};

TEST(RefinedSparseCholesky, Factorize) {
  auto* mock_sparse_cholesky = new MockSparseCholesky;
  auto* mock_iterative_refiner = new MockSparseIterativeRefiner;
  EXPECT_CALL(*mock_sparse_cholesky, Factorize(_, _))
      .Times(1)
      .WillRepeatedly(Return(LinearSolverTerminationType::SUCCESS));
  EXPECT_CALL(*mock_iterative_refiner, Refine(_, _, _, _)).Times(0);
  std::unique_ptr<SparseCholesky> sparse_cholesky(mock_sparse_cholesky);
  std::unique_ptr<SparseIterativeRefiner> iterative_refiner(
      mock_iterative_refiner);
  RefinedSparseCholesky refined_sparse_cholesky(std::move(sparse_cholesky),
                                                std::move(iterative_refiner));
  CompressedRowSparseMatrix m(1, 1, 1);
  std::string message;
  EXPECT_EQ(refined_sparse_cholesky.Factorize(&m, &message),
            LinearSolverTerminationType::SUCCESS);
};

TEST(RefinedSparseCholesky, FactorAndSolveWithUnsuccessfulFactorization) {
  auto* mock_sparse_cholesky = new MockSparseCholesky;
  auto* mock_iterative_refiner = new MockSparseIterativeRefiner;
  EXPECT_CALL(*mock_sparse_cholesky, Factorize(_, _))
      .Times(1)
      .WillRepeatedly(Return(LinearSolverTerminationType::FAILURE));
  EXPECT_CALL(*mock_sparse_cholesky, Solve(_, _, _)).Times(0);
  EXPECT_CALL(*mock_iterative_refiner, Refine(_, _, _, _)).Times(0);
  std::unique_ptr<SparseCholesky> sparse_cholesky(mock_sparse_cholesky);
  std::unique_ptr<SparseIterativeRefiner> iterative_refiner(
      mock_iterative_refiner);
  RefinedSparseCholesky refined_sparse_cholesky(std::move(sparse_cholesky),
                                                std::move(iterative_refiner));
  CompressedRowSparseMatrix m(1, 1, 1);
  std::string message;
  double rhs;
  double solution;
  EXPECT_EQ(
      refined_sparse_cholesky.FactorAndSolve(&m, &rhs, &solution, &message),
      LinearSolverTerminationType::FAILURE);
};

TEST(RefinedSparseCholesky, FactorAndSolveWithSuccess) {
  auto* mock_sparse_cholesky = new MockSparseCholesky;
  std::unique_ptr<MockSparseIterativeRefiner> mock_iterative_refiner(
      new MockSparseIterativeRefiner);
  EXPECT_CALL(*mock_sparse_cholesky, Factorize(_, _))
      .Times(1)
      .WillRepeatedly(Return(LinearSolverTerminationType::SUCCESS));
  EXPECT_CALL(*mock_sparse_cholesky, Solve(_, _, _))
      .Times(1)
      .WillRepeatedly(Return(LinearSolverTerminationType::SUCCESS));
  EXPECT_CALL(*mock_iterative_refiner, Refine(_, _, _, _)).Times(1);

  std::unique_ptr<SparseCholesky> sparse_cholesky(mock_sparse_cholesky);
  std::unique_ptr<SparseIterativeRefiner> iterative_refiner(
      std::move(mock_iterative_refiner));
  RefinedSparseCholesky refined_sparse_cholesky(std::move(sparse_cholesky),
                                                std::move(iterative_refiner));
  CompressedRowSparseMatrix m(1, 1, 1);
  std::string message;
  double rhs;
  double solution;
  EXPECT_EQ(
      refined_sparse_cholesky.FactorAndSolve(&m, &rhs, &solution, &message),
      LinearSolverTerminationType::SUCCESS);
};

}  // namespace

}  // namespace ceres::internal
