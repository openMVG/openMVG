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

#include "ceres/schur_eliminator.h"

#include <algorithm>
#include <memory>
#include <random>
#include <vector>

#include "Eigen/Dense"
#include "ceres/block_random_access_dense_matrix.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/block_structure.h"
#include "ceres/casts.h"
#include "ceres/context_impl.h"
#include "ceres/detect_structure.h"
#include "ceres/internal/eigen.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/test_util.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

// TODO(sameeragarwal): Reduce the size of these tests and redo the
// parameterization to be more efficient.

namespace ceres::internal {

class SchurEliminatorTest : public ::testing::Test {
 protected:
  void SetUpFromId(int id) {
    auto problem = CreateLinearLeastSquaresProblemFromId(id);
    CHECK(problem != nullptr);
    SetupHelper(problem.get());
  }

  void SetupHelper(LinearLeastSquaresProblem* problem) {
    A.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));
    b = std::move(problem->b);
    D = std::move(problem->D);

    num_eliminate_blocks = problem->num_eliminate_blocks;
    num_eliminate_cols = 0;
    const CompressedRowBlockStructure* bs = A->block_structure();

    for (int i = 0; i < num_eliminate_blocks; ++i) {
      num_eliminate_cols += bs->cols[i].size;
    }
  }

  // Compute the golden values for the reduced linear system and the
  // solution to the linear least squares problem using dense linear
  // algebra.
  void ComputeReferenceSolution(const Vector& D) {
    Matrix J;
    A->ToDenseMatrix(&J);
    VectorRef f(b.get(), J.rows());

    Matrix H = (D.cwiseProduct(D)).asDiagonal();
    H.noalias() += J.transpose() * J;

    const Vector g = J.transpose() * f;
    const int schur_size = J.cols() - num_eliminate_cols;

    lhs_expected.resize(schur_size, schur_size);
    lhs_expected.setZero();

    rhs_expected.resize(schur_size);
    rhs_expected.setZero();

    sol_expected.resize(J.cols());
    sol_expected.setZero();

    Matrix P = H.block(0, 0, num_eliminate_cols, num_eliminate_cols);
    Matrix Q = H.block(0, num_eliminate_cols, num_eliminate_cols, schur_size);
    Matrix R =
        H.block(num_eliminate_cols, num_eliminate_cols, schur_size, schur_size);
    int row = 0;
    const CompressedRowBlockStructure* bs = A->block_structure();
    for (int i = 0; i < num_eliminate_blocks; ++i) {
      const int block_size = bs->cols[i].size;
      P.block(row, row, block_size, block_size) =
          P.block(row, row, block_size, block_size)
              .llt()
              .solve(Matrix::Identity(block_size, block_size));
      row += block_size;
    }

    lhs_expected.triangularView<Eigen::Upper>() = R - Q.transpose() * P * Q;
    rhs_expected =
        g.tail(schur_size) - Q.transpose() * P * g.head(num_eliminate_cols);
    sol_expected = H.llt().solve(g);
  }

  void EliminateSolveAndCompare(const VectorRef& diagonal,
                                bool use_static_structure,
                                const double relative_tolerance) {
    const CompressedRowBlockStructure* bs = A->block_structure();
    const int num_col_blocks = bs->cols.size();
    auto blocks = Tail(bs->cols, num_col_blocks - num_eliminate_blocks);
    BlockRandomAccessDenseMatrix lhs(blocks, &context_, 1);

    const int num_cols = A->num_cols();
    const int schur_size = lhs.num_rows();

    Vector rhs(schur_size);

    LinearSolver::Options options;
    options.context = &context_;
    options.elimination_groups.push_back(num_eliminate_blocks);
    if (use_static_structure) {
      DetectStructure(*bs,
                      num_eliminate_blocks,
                      &options.row_block_size,
                      &options.e_block_size,
                      &options.f_block_size);
    }

    std::unique_ptr<SchurEliminatorBase> eliminator =
        SchurEliminatorBase::Create(options);
    const bool kFullRankETE = true;
    eliminator->Init(num_eliminate_blocks, kFullRankETE, A->block_structure());
    eliminator->Eliminate(
        BlockSparseMatrixData(*A), b.get(), diagonal.data(), &lhs, rhs.data());

    MatrixRef lhs_ref(lhs.mutable_values(), lhs.num_rows(), lhs.num_cols());
    Vector reduced_sol =
        lhs_ref.selfadjointView<Eigen::Upper>().llt().solve(rhs);

    // Solution to the linear least squares problem.
    Vector sol(num_cols);
    sol.setZero();
    sol.tail(schur_size) = reduced_sol;
    eliminator->BackSubstitute(BlockSparseMatrixData(*A),
                               b.get(),
                               diagonal.data(),
                               reduced_sol.data(),
                               sol.data());

    Matrix delta = (lhs_ref - lhs_expected).selfadjointView<Eigen::Upper>();
    double diff = delta.norm();
    EXPECT_NEAR(diff / lhs_expected.norm(), 0.0, relative_tolerance);
    EXPECT_NEAR((rhs - rhs_expected).norm() / rhs_expected.norm(),
                0.0,
                relative_tolerance);
    EXPECT_NEAR((sol - sol_expected).norm() / sol_expected.norm(),
                0.0,
                relative_tolerance);
  }

  ContextImpl context_;

  std::unique_ptr<BlockSparseMatrix> A;
  std::unique_ptr<double[]> b;
  std::unique_ptr<double[]> D;
  int num_eliminate_blocks;
  int num_eliminate_cols;

  Matrix lhs_expected;
  Vector rhs_expected;
  Vector sol_expected;
};

TEST_F(SchurEliminatorTest, ScalarProblemNoRegularization) {
  SetUpFromId(2);
  Vector zero(A->num_cols());
  zero.setZero();

  ComputeReferenceSolution(VectorRef(zero.data(), A->num_cols()));
  EliminateSolveAndCompare(VectorRef(zero.data(), A->num_cols()), true, 1e-14);
  EliminateSolveAndCompare(VectorRef(zero.data(), A->num_cols()), false, 1e-14);
}

TEST_F(SchurEliminatorTest, ScalarProblemWithRegularization) {
  SetUpFromId(2);
  ComputeReferenceSolution(VectorRef(D.get(), A->num_cols()));
  EliminateSolveAndCompare(VectorRef(D.get(), A->num_cols()), true, 1e-14);
  EliminateSolveAndCompare(VectorRef(D.get(), A->num_cols()), false, 1e-14);
}

TEST_F(SchurEliminatorTest, VaryingFBlockSizeWithStaticStructure) {
  SetUpFromId(4);
  ComputeReferenceSolution(VectorRef(D.get(), A->num_cols()));
  EliminateSolveAndCompare(VectorRef(D.get(), A->num_cols()), true, 1e-14);
}

TEST_F(SchurEliminatorTest, VaryingFBlockSizeWithoutStaticStructure) {
  SetUpFromId(4);
  ComputeReferenceSolution(VectorRef(D.get(), A->num_cols()));
  EliminateSolveAndCompare(VectorRef(D.get(), A->num_cols()), false, 1e-14);
}

TEST(SchurEliminatorForOneFBlock, MatchesSchurEliminator) {
  constexpr int kRowBlockSize = 2;
  constexpr int kEBlockSize = 3;
  constexpr int kFBlockSize = 6;
  constexpr int num_e_blocks = 5;

  ContextImpl context;

  auto* bs = new CompressedRowBlockStructure;
  bs->cols.resize(num_e_blocks + 1);
  int col_pos = 0;
  for (int i = 0; i < num_e_blocks; ++i) {
    bs->cols[i].position = col_pos;
    bs->cols[i].size = kEBlockSize;
    col_pos += kEBlockSize;
  }
  bs->cols.back().position = col_pos;
  bs->cols.back().size = kFBlockSize;

  bs->rows.resize(2 * num_e_blocks + 1);
  int row_pos = 0;
  int cell_pos = 0;
  for (int i = 0; i < num_e_blocks; ++i) {
    {
      auto& row = bs->rows[2 * i];
      row.block.position = row_pos;
      row.block.size = kRowBlockSize;
      row_pos += kRowBlockSize;
      auto& cells = row.cells;
      cells.resize(2);
      cells[0].block_id = i;
      cells[0].position = cell_pos;
      cell_pos += kRowBlockSize * kEBlockSize;
      cells[1].block_id = num_e_blocks;
      cells[1].position = cell_pos;
      cell_pos += kRowBlockSize * kFBlockSize;
    }
    {
      auto& row = bs->rows[2 * i + 1];
      row.block.position = row_pos;
      row.block.size = kRowBlockSize;
      row_pos += kRowBlockSize;
      auto& cells = row.cells;
      cells.resize(1);
      cells[0].block_id = i;
      cells[0].position = cell_pos;
      cell_pos += kRowBlockSize * kEBlockSize;
    }
  }

  {
    auto& row = bs->rows.back();
    row.block.position = row_pos;
    row.block.size = kEBlockSize;
    row_pos += kRowBlockSize;
    auto& cells = row.cells;
    cells.resize(1);
    cells[0].block_id = num_e_blocks;
    cells[0].position = cell_pos;
    cell_pos += kEBlockSize * kEBlockSize;
  }

  BlockSparseMatrix matrix(bs);
  double* values = matrix.mutable_values();
  std::mt19937 prng;
  std::normal_distribution<> standard_normal;
  std::generate_n(values, matrix.num_nonzeros(), [&prng, &standard_normal] {
    return standard_normal(prng);
  });

  Vector b(matrix.num_rows());
  b.setRandom();

  Vector diagonal(matrix.num_cols());
  diagonal.setOnes();

  std::vector<Block> blocks;
  blocks.emplace_back(kFBlockSize, 0);
  BlockRandomAccessDenseMatrix actual_lhs(blocks, &context, 1);
  BlockRandomAccessDenseMatrix expected_lhs(blocks, &context, 1);
  Vector actual_rhs(kFBlockSize);
  Vector expected_rhs(kFBlockSize);

  Vector f_sol(kFBlockSize);
  f_sol.setRandom();
  Vector actual_e_sol(num_e_blocks * kEBlockSize);
  actual_e_sol.setZero();
  Vector expected_e_sol(num_e_blocks * kEBlockSize);
  expected_e_sol.setZero();

  {
    LinearSolver::Options linear_solver_options;
    linear_solver_options.e_block_size = kEBlockSize;
    linear_solver_options.row_block_size = kRowBlockSize;
    linear_solver_options.f_block_size = kFBlockSize;
    linear_solver_options.context = &context;
    std::unique_ptr<SchurEliminatorBase> eliminator(
        SchurEliminatorBase::Create(linear_solver_options));
    eliminator->Init(num_e_blocks, true, matrix.block_structure());
    eliminator->Eliminate(BlockSparseMatrixData(matrix),
                          b.data(),
                          diagonal.data(),
                          &expected_lhs,
                          expected_rhs.data());
    eliminator->BackSubstitute(BlockSparseMatrixData(matrix),
                               b.data(),
                               diagonal.data(),
                               f_sol.data(),
                               actual_e_sol.data());
  }

  {
    SchurEliminatorForOneFBlock<2, 3, 6> eliminator;
    eliminator.Init(num_e_blocks, true, matrix.block_structure());
    eliminator.Eliminate(BlockSparseMatrixData(matrix),
                         b.data(),
                         diagonal.data(),
                         &actual_lhs,
                         actual_rhs.data());
    eliminator.BackSubstitute(BlockSparseMatrixData(matrix),
                              b.data(),
                              diagonal.data(),
                              f_sol.data(),
                              expected_e_sol.data());
  }
  ConstMatrixRef actual_lhsref(
      actual_lhs.values(), actual_lhs.num_cols(), actual_lhs.num_cols());
  ConstMatrixRef expected_lhsref(
      expected_lhs.values(), actual_lhs.num_cols(), actual_lhs.num_cols());

  EXPECT_NEAR((actual_lhsref - expected_lhsref).norm() / expected_lhsref.norm(),
              0.0,
              1e-12)
      << "expected: \n"
      << expected_lhsref << "\nactual: \n"
      << actual_lhsref;

  EXPECT_NEAR(
      (actual_rhs - expected_rhs).norm() / expected_rhs.norm(), 0.0, 1e-12)
      << "expected: \n"
      << expected_rhs << "\nactual: \n"
      << actual_rhs;

  EXPECT_NEAR((actual_e_sol - expected_e_sol).norm() / expected_e_sol.norm(),
              0.0,
              1e-12)
      << "expected: \n"
      << expected_e_sol << "\nactual: \n"
      << actual_e_sol;
}

}  // namespace ceres::internal
