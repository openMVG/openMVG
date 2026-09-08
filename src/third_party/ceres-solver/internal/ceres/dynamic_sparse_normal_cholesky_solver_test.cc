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

#include <memory>

#include "Eigen/Cholesky"
#include "ceres/casts.h"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/context_impl.h"
#include "ceres/internal/config.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

// TODO(sameeragarwal): These tests needs to be re-written to be more
// thorough, they do not really test the dynamic nature of the
// sparsity.
class DynamicSparseNormalCholeskySolverTest : public ::testing::Test {
 protected:
  void SetUp() final {
    std::unique_ptr<LinearLeastSquaresProblem> problem =
        CreateLinearLeastSquaresProblemFromId(1);
    A_ = CompressedRowSparseMatrix::FromTripletSparseMatrix(
        *down_cast<TripletSparseMatrix*>(problem->A.get()));
    b_ = std::move(problem->b);
    D_ = std::move(problem->D);
  }

  void TestSolver(const LinearSolver::Options& options, double* D) {
    Matrix dense_A;
    A_->ToDenseMatrix(&dense_A);
    Matrix lhs = dense_A.transpose() * dense_A;
    if (D != nullptr) {
      lhs += (ConstVectorRef(D, A_->num_cols()).array() *
              ConstVectorRef(D, A_->num_cols()).array())
                 .matrix()
                 .asDiagonal();
    }

    Vector rhs(A_->num_cols());
    rhs.setZero();
    A_->LeftMultiplyAndAccumulate(b_.get(), rhs.data());
    Vector expected_solution = lhs.llt().solve(rhs);

    std::unique_ptr<LinearSolver> solver(LinearSolver::Create(options));
    LinearSolver::PerSolveOptions per_solve_options;
    per_solve_options.D = D;
    Vector actual_solution(A_->num_cols());
    LinearSolver::Summary summary;
    summary = solver->Solve(
        A_.get(), b_.get(), per_solve_options, actual_solution.data());

    EXPECT_EQ(summary.termination_type, LinearSolverTerminationType::SUCCESS);

    for (int i = 0; i < A_->num_cols(); ++i) {
      EXPECT_NEAR(expected_solution(i), actual_solution(i), 1e-8)
          << "\nExpected: " << expected_solution.transpose()
          << "\nActual: " << actual_solution.transpose();
    }
  }

  void TestSolver(
      const SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type,
      const OrderingType ordering_type) {
    LinearSolver::Options options;
    options.type = SPARSE_NORMAL_CHOLESKY;
    options.dynamic_sparsity = true;
    options.sparse_linear_algebra_library_type =
        sparse_linear_algebra_library_type;
    options.ordering_type = ordering_type;
    ContextImpl context;
    options.context = &context;
    TestSolver(options, nullptr);
    TestSolver(options, D_.get());
  }

  std::unique_ptr<CompressedRowSparseMatrix> A_;
  std::unique_ptr<double[]> b_;
  std::unique_ptr<double[]> D_;
};

#ifndef CERES_NO_SUITESPARSE
TEST_F(DynamicSparseNormalCholeskySolverTest, SuiteSparseAMD) {
  TestSolver(SUITE_SPARSE, OrderingType::AMD);
}

#ifndef CERES_NO_CHOLMOD_PARTITION
TEST_F(DynamicSparseNormalCholeskySolverTest, SuiteSparseNESDIS) {
  TestSolver(SUITE_SPARSE, OrderingType::NESDIS);
}
#endif
#endif

#ifdef CERES_USE_EIGEN_SPARSE
TEST_F(DynamicSparseNormalCholeskySolverTest, EigenAMD) {
  TestSolver(EIGEN_SPARSE, OrderingType::AMD);
}

#ifndef CERES_NO_EIGEN_METIS
TEST_F(DynamicSparseNormalCholeskySolverTest, EigenNESDIS) {
  TestSolver(EIGEN_SPARSE, OrderingType::NESDIS);
}
#endif
#endif  // CERES_USE_EIGEN_SPARSE

}  // namespace internal
}  // namespace ceres
