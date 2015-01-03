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

#include "ceres/casts.h"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"


namespace ceres {
namespace internal {

class UnsymmetricLinearSolverTest : public ::testing::Test {
 protected :
  virtual void SetUp() {
    scoped_ptr<LinearLeastSquaresProblem> problem(
        CreateLinearLeastSquaresProblemFromId(0));

    CHECK_NOTNULL(problem.get());
    A_.reset(down_cast<TripletSparseMatrix*>(problem->A.release()));
    b_.reset(problem->b.release());
    D_.reset(problem->D.release());
    sol_unregularized_.reset(problem->x.release());
    sol_regularized_.reset(problem->x_D.release());
  }

  void TestSolver(const LinearSolver::Options& options) {


    LinearSolver::PerSolveOptions per_solve_options;
    LinearSolver::Summary unregularized_solve_summary;
    LinearSolver::Summary regularized_solve_summary;
    Vector x_unregularized(A_->num_cols());
    Vector x_regularized(A_->num_cols());

    scoped_ptr<SparseMatrix> transformed_A;

    if (options.type == DENSE_QR ||
        options.type == DENSE_NORMAL_CHOLESKY) {
      transformed_A.reset(new DenseSparseMatrix(*A_));
    } else if (options.type == SPARSE_NORMAL_CHOLESKY) {
      CompressedRowSparseMatrix* crsm =  new CompressedRowSparseMatrix(*A_);
      // Add row/column blocks structure.
      for (int i = 0; i < A_->num_rows(); ++i) {
        crsm->mutable_row_blocks()->push_back(1);
      }

      for (int i = 0; i < A_->num_cols(); ++i) {
        crsm->mutable_col_blocks()->push_back(1);
      }
      transformed_A.reset(crsm);
    } else {
      LOG(FATAL) << "Unknown linear solver : " << options.type;
    }

    // Unregularized
    scoped_ptr<LinearSolver> solver(LinearSolver::Create(options));
    unregularized_solve_summary =
        solver->Solve(transformed_A.get(),
                      b_.get(),
                      per_solve_options,
                      x_unregularized.data());

    // Sparsity structure is changing, reset the solver.
    solver.reset(LinearSolver::Create(options));
    // Regularized solution
    per_solve_options.D = D_.get();
    regularized_solve_summary =
        solver->Solve(transformed_A.get(),
                      b_.get(),
                      per_solve_options,
                      x_regularized.data());

    EXPECT_EQ(unregularized_solve_summary.termination_type,
              LINEAR_SOLVER_SUCCESS);

    for (int i = 0; i < A_->num_cols(); ++i) {
      EXPECT_NEAR(sol_unregularized_[i], x_unregularized[i], 1e-8)
          << "\nExpected: "
          << ConstVectorRef(sol_unregularized_.get(), A_->num_cols()).transpose()
          << "\nActual: " << x_unregularized.transpose();
    }

    EXPECT_EQ(regularized_solve_summary.termination_type,
              LINEAR_SOLVER_SUCCESS);
    for (int i = 0; i < A_->num_cols(); ++i) {
      EXPECT_NEAR(sol_regularized_[i], x_regularized[i], 1e-8)
          << "\nExpected: "
          << ConstVectorRef(sol_regularized_.get(), A_->num_cols()).transpose()
          << "\nActual: " << x_regularized.transpose();
    }
  }

  scoped_ptr<TripletSparseMatrix> A_;
  scoped_array<double> b_;
  scoped_array<double> D_;
  scoped_array<double> sol_unregularized_;
  scoped_array<double> sol_regularized_;
};

TEST_F(UnsymmetricLinearSolverTest, EigenDenseQR) {
  LinearSolver::Options options;
  options.type = DENSE_QR;
  options.dense_linear_algebra_library_type = EIGEN;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest, EigenDenseNormalCholesky) {
  LinearSolver::Options options;
  options.dense_linear_algebra_library_type = EIGEN;
  options.type = DENSE_NORMAL_CHOLESKY;
  TestSolver(options);
}

#ifndef CERES_NO_LAPACK
TEST_F(UnsymmetricLinearSolverTest, LAPACKDenseQR) {
  LinearSolver::Options options;
  options.type = DENSE_QR;
  options.dense_linear_algebra_library_type = LAPACK;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest, LAPACKDenseNormalCholesky) {
  LinearSolver::Options options;
  options.dense_linear_algebra_library_type = LAPACK;
  options.type = DENSE_NORMAL_CHOLESKY;
  TestSolver(options);
}
#endif

#ifndef CERES_NO_SUITESPARSE
TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingSuiteSparsePreOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = false;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingSuiteSparsePostOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = true;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingSuiteSparseDynamicSparsity) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.dynamic_sparsity = true;
  TestSolver(options);
}
#endif

#ifndef CERES_NO_CXSPARSE
TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingCXSparsePreOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = false;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingCXSparsePostOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = true;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingCXSparseDynamicSparsity) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.dynamic_sparsity = true;
  TestSolver(options);
}
#endif

#ifdef CERES_USE_EIGEN_SPARSE
TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingEigenPreOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = false;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingEigenPostOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = true;
  TestSolver(options);
}

TEST_F(UnsymmetricLinearSolverTest,
       SparseNormalCholeskyUsingEigenDynamicSparsity) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.dynamic_sparsity = true;
  TestSolver(options);
}

#endif  // CERES_USE_EIGEN_SPARSE

}  // namespace internal
}  // namespace ceres
