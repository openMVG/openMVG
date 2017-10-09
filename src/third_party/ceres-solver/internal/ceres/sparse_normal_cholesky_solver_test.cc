// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2017 Google Inc. All rights reserved.
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

#include "ceres/block_sparse_matrix.h"
#include "ceres/casts.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

#include "Eigen/Cholesky"

namespace ceres {
namespace internal {

// TODO(sameeragarwal): These tests needs to be re-written, since
// SparseNormalCholeskySolver is a composition of two classes now,
// InnerProductComputer and SparseCholesky.
//
// So the test should exercise the composition, rather than the
// numerics of the solver, which are well covered by tests for those
// classes.
class SparseNormalCholeskySolverTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    scoped_ptr<LinearLeastSquaresProblem> problem(
        CreateLinearLeastSquaresProblemFromId(2));

    CHECK_NOTNULL(problem.get());
    A_.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));
    b_.reset(problem->b.release());
    D_.reset(problem->D.release());
  }

  void TestSolver(const LinearSolver::Options& options, double* D) {
    Matrix dense_A;
    A_->ToDenseMatrix(&dense_A);
    Matrix lhs = dense_A.transpose() * dense_A;
    if (D != NULL) {
      lhs += (ConstVectorRef(D, A_->num_cols()).array() *
              ConstVectorRef(D, A_->num_cols()).array())
                 .matrix()
                 .asDiagonal();
    }

    Vector rhs(A_->num_cols());
    rhs.setZero();
    A_->LeftMultiply(b_.get(), rhs.data());
    Vector expected_solution = lhs.llt().solve(rhs);

    scoped_ptr<LinearSolver> solver(LinearSolver::Create(options));
    LinearSolver::PerSolveOptions per_solve_options;
    per_solve_options.D = D;
    Vector actual_solution(A_->num_cols());
    LinearSolver::Summary summary;
    summary = solver->Solve(
        A_.get(), b_.get(), per_solve_options, actual_solution.data());

    EXPECT_EQ(summary.termination_type, LINEAR_SOLVER_SUCCESS);

    for (int i = 0; i < A_->num_cols(); ++i) {
      EXPECT_NEAR(expected_solution(i), actual_solution(i), 1e-8)
          << "\nExpected: " << expected_solution.transpose()
          << "\nActual: " << actual_solution.transpose();
    }
  }

  void TestSolver(const LinearSolver::Options& options) {
    TestSolver(options, NULL);
    TestSolver(options, D_.get());
  }

  scoped_ptr<BlockSparseMatrix> A_;
  scoped_array<double> b_;
  scoped_array<double> D_;
};

#ifndef CERES_NO_SUITESPARSE
TEST_F(SparseNormalCholeskySolverTest,
       SparseNormalCholeskyUsingSuiteSparsePreOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = false;
  TestSolver(options);
}

TEST_F(SparseNormalCholeskySolverTest,
       SparseNormalCholeskyUsingSuiteSparsePostOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = true;
  TestSolver(options);
}
#endif

#ifndef CERES_NO_CXSPARSE
TEST_F(SparseNormalCholeskySolverTest,
       SparseNormalCholeskyUsingCXSparsePreOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = false;
  TestSolver(options);
}

TEST_F(SparseNormalCholeskySolverTest,
       SparseNormalCholeskyUsingCXSparsePostOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = true;
  TestSolver(options);
}
#endif

#ifdef CERES_USE_EIGEN_SPARSE
TEST_F(SparseNormalCholeskySolverTest,
       SparseNormalCholeskyUsingEigenPreOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = false;
  TestSolver(options);
}

TEST_F(SparseNormalCholeskySolverTest,
       SparseNormalCholeskyUsingEigenPostOrdering) {
  LinearSolver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.type = SPARSE_NORMAL_CHOLESKY;
  options.use_postordering = true;
  TestSolver(options);
}
#endif  // CERES_USE_EIGEN_SPARSE

}  // namespace internal
}  // namespace ceres
