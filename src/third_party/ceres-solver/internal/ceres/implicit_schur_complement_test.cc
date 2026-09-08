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

#include "ceres/implicit_schur_complement.h"

#include <cstddef>
#include <memory>

#include "Eigen/Dense"
#include "ceres/block_random_access_dense_matrix.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/casts.h"
#include "ceres/context_impl.h"
#include "ceres/internal/eigen.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/schur_eliminator.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

using testing::AssertionResult;

const double kEpsilon = 1e-14;

class ImplicitSchurComplementTest : public ::testing::Test {
 protected:
  void SetUp() final {
    auto problem = CreateLinearLeastSquaresProblemFromId(2);

    CHECK(problem != nullptr);
    A_.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));
    b_ = std::move(problem->b);
    D_ = std::move(problem->D);

    num_cols_ = A_->num_cols();
    num_rows_ = A_->num_rows();
    num_eliminate_blocks_ = problem->num_eliminate_blocks;
  }

  void ReducedLinearSystemAndSolution(double* D,
                                      Matrix* lhs,
                                      Vector* rhs,
                                      Vector* solution) {
    const CompressedRowBlockStructure* bs = A_->block_structure();
    const int num_col_blocks = bs->cols.size();
    auto blocks = Tail(bs->cols, num_col_blocks - num_eliminate_blocks_);
    BlockRandomAccessDenseMatrix blhs(blocks, &context_, 1);
    const int num_schur_rows = blhs.num_rows();

    LinearSolver::Options options;
    options.elimination_groups.push_back(num_eliminate_blocks_);
    options.type = DENSE_SCHUR;
    ContextImpl context;
    options.context = &context;

    std::unique_ptr<SchurEliminatorBase> eliminator =
        SchurEliminatorBase::Create(options);
    CHECK(eliminator != nullptr);
    const bool kFullRankETE = true;
    eliminator->Init(num_eliminate_blocks_, kFullRankETE, bs);

    lhs->resize(num_schur_rows, num_schur_rows);
    rhs->resize(num_schur_rows);

    eliminator->Eliminate(
        BlockSparseMatrixData(*A_), b_.get(), D, &blhs, rhs->data());

    MatrixRef lhs_ref(blhs.mutable_values(), num_schur_rows, num_schur_rows);

    // lhs_ref is an upper triangular matrix. Construct a full version
    // of lhs_ref in lhs by transposing lhs_ref, choosing the strictly
    // lower triangular part of the matrix and adding it to lhs_ref.
    *lhs = lhs_ref;
    lhs->triangularView<Eigen::StrictlyLower>() =
        lhs_ref.triangularView<Eigen::StrictlyUpper>().transpose();

    solution->resize(num_cols_);
    solution->setZero();
    VectorRef schur_solution(solution->data() + num_cols_ - num_schur_rows,
                             num_schur_rows);
    schur_solution = lhs->selfadjointView<Eigen::Upper>().llt().solve(*rhs);
    eliminator->BackSubstitute(BlockSparseMatrixData(*A_),
                               b_.get(),
                               D,
                               schur_solution.data(),
                               solution->data());
  }

  AssertionResult TestImplicitSchurComplement(double* D) {
    Matrix lhs;
    Vector rhs;
    Vector reference_solution;
    ReducedLinearSystemAndSolution(D, &lhs, &rhs, &reference_solution);

    LinearSolver::Options options;
    options.elimination_groups.push_back(num_eliminate_blocks_);
    options.preconditioner_type = JACOBI;
    ContextImpl context;
    options.context = &context;
    ImplicitSchurComplement isc(options);
    isc.Init(*A_, D, b_.get());

    const int num_f_cols = lhs.cols();
    const int num_e_cols = num_cols_ - num_f_cols;

    Matrix A_dense, E, F, DE, DF;
    A_->ToDenseMatrix(&A_dense);
    E = A_dense.leftCols(A_->num_cols() - num_f_cols);
    F = A_dense.rightCols(num_f_cols);
    if (D) {
      DE = VectorRef(D, num_e_cols).asDiagonal();
      DF = VectorRef(D + num_e_cols, num_f_cols).asDiagonal();
    } else {
      DE = Matrix::Zero(num_e_cols, num_e_cols);
      DF = Matrix::Zero(num_f_cols, num_f_cols);
    }

    // Z = (block_diagonal(F'F))^-1 F'E (E'E)^-1 E'F
    // Here, assuming that block_diagonal(F'F) == diagonal(F'F)
    Matrix Z_reference =
        (F.transpose() * F + DF).diagonal().asDiagonal().inverse() *
        F.transpose() * E * (E.transpose() * E + DE).inverse() * E.transpose() *
        F;

    for (int i = 0; i < num_f_cols; ++i) {
      Vector x(num_f_cols);
      x.setZero();
      x(i) = 1.0;

      Vector y(num_f_cols);
      y = lhs * x;

      Vector z(num_f_cols);
      isc.RightMultiplyAndAccumulate(x.data(), z.data());

      // The i^th column of the implicit schur complement is the same as
      // the explicit schur complement.
      if ((y - z).norm() > kEpsilon) {
        return testing::AssertionFailure()
               << "Explicit and Implicit SchurComplements differ in "
               << "column " << i << ". explicit: " << y.transpose()
               << " implicit: " << z.transpose();
      }

      y.setZero();
      y = Z_reference * x;
      z.setZero();
      isc.InversePowerSeriesOperatorRightMultiplyAccumulate(x.data(), z.data());

      // The i^th column of operator Z stored implicitly is the same as its
      // explicit version.
      if ((y - z).norm() > kEpsilon) {
        return testing::AssertionFailure()
               << "Explicit and Implicit operators used to approximate the "
                  "inversion of schur complement via power series expansion "
                  "differ in column "
               << i << ". explicit: " << y.transpose()
               << " implicit: " << z.transpose();
      }
    }

    // Compare the rhs of the reduced linear system
    if ((isc.rhs() - rhs).norm() > kEpsilon) {
      return testing::AssertionFailure()
             << "Explicit and Implicit SchurComplements differ in "
             << "rhs. explicit: " << rhs.transpose()
             << " implicit: " << isc.rhs().transpose();
    }

    // Reference solution to the f_block.
    const Vector reference_f_sol =
        lhs.selfadjointView<Eigen::Upper>().llt().solve(rhs);

    // Backsubstituted solution from the implicit schur solver using the
    // reference solution to the f_block.
    Vector sol(num_cols_);
    isc.BackSubstitute(reference_f_sol.data(), sol.data());
    if ((sol - reference_solution).norm() > kEpsilon) {
      return testing::AssertionFailure()
             << "Explicit and Implicit SchurComplements solutions differ. "
             << "explicit: " << reference_solution.transpose()
             << " implicit: " << sol.transpose();
    }

    return testing::AssertionSuccess();
  }

  ContextImpl context_;
  int num_rows_;
  int num_cols_;
  int num_eliminate_blocks_;

  std::unique_ptr<BlockSparseMatrix> A_;
  std::unique_ptr<double[]> b_;
  std::unique_ptr<double[]> D_;
};

// Verify that the Schur Complement matrix implied by the
// ImplicitSchurComplement class matches the one explicitly computed
// by the SchurComplement solver.
//
// We do this with and without regularization to check that the
// support for the LM diagonal is correct.
TEST_F(ImplicitSchurComplementTest, SchurMatrixValuesTest) {
  EXPECT_TRUE(TestImplicitSchurComplement(nullptr));
  EXPECT_TRUE(TestImplicitSchurComplement(D_.get()));
}

}  // namespace ceres::internal
