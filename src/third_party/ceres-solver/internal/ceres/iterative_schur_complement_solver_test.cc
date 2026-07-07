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
//
// TODO(sameeragarwal): Add support for larger, more complicated and
// poorly conditioned problems both for correctness testing as well as
// benchmarking.

#include "ceres/iterative_schur_complement_solver.h"

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

namespace ceres {
namespace internal {

using testing::AssertionResult;

const double kEpsilon = 1e-14;

class IterativeSchurComplementSolverTest : public ::testing::Test {
 protected:
  void SetUpProblem(int problem_id) {
    std::unique_ptr<LinearLeastSquaresProblem> problem =
        CreateLinearLeastSquaresProblemFromId(problem_id);

    CHECK(problem != nullptr);
    A_.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));
    b_ = std::move(problem->b);
    D_ = std::move(problem->D);

    num_cols_ = A_->num_cols();
    num_rows_ = A_->num_rows();
    num_eliminate_blocks_ = problem->num_eliminate_blocks;
  }

  AssertionResult TestSolver(double* D,
                             PreconditionerType preconditioner_type,
                             bool use_spse_initialization) {
    TripletSparseMatrix triplet_A(
        A_->num_rows(), A_->num_cols(), A_->num_nonzeros());
    A_->ToTripletSparseMatrix(&triplet_A);

    DenseSparseMatrix dense_A(triplet_A);

    LinearSolver::Options options;
    options.type = DENSE_QR;
    ContextImpl context;
    options.context = &context;
    std::unique_ptr<LinearSolver> qr(LinearSolver::Create(options));

    LinearSolver::PerSolveOptions per_solve_options;
    per_solve_options.D = D;
    Vector reference_solution(num_cols_);
    qr->Solve(&dense_A, b_.get(), per_solve_options, reference_solution.data());

    options.elimination_groups.push_back(num_eliminate_blocks_);
    options.elimination_groups.push_back(0);
    options.max_num_iterations = num_cols_;
    options.max_num_spse_iterations = 1;
    options.use_spse_initialization = use_spse_initialization;
    options.preconditioner_type = preconditioner_type;
    IterativeSchurComplementSolver isc(options);

    Vector isc_sol(num_cols_);
    per_solve_options.r_tolerance = 1e-12;
    isc.Solve(A_.get(), b_.get(), per_solve_options, isc_sol.data());
    double diff = (isc_sol - reference_solution).norm();
    if (diff < kEpsilon) {
      return testing::AssertionSuccess();
    } else {
      return testing::AssertionFailure()
             << "The reference solution differs from the ITERATIVE_SCHUR"
             << " solution by " << diff << " which is more than " << kEpsilon;
    }
  }

  int num_rows_;
  int num_cols_;
  int num_eliminate_blocks_;
  std::unique_ptr<BlockSparseMatrix> A_;
  std::unique_ptr<double[]> b_;
  std::unique_ptr<double[]> D_;
};

TEST_F(IterativeSchurComplementSolverTest, NormalProblemSchurJacobi) {
  SetUpProblem(2);
  EXPECT_TRUE(TestSolver(nullptr, SCHUR_JACOBI, false));
  EXPECT_TRUE(TestSolver(D_.get(), SCHUR_JACOBI, false));
}

TEST_F(IterativeSchurComplementSolverTest,
       NormalProblemSchurJacobiWithPowerSeriesExpansionInitialization) {
  SetUpProblem(2);
  EXPECT_TRUE(TestSolver(nullptr, SCHUR_JACOBI, true));
  EXPECT_TRUE(TestSolver(D_.get(), SCHUR_JACOBI, true));
}

TEST_F(IterativeSchurComplementSolverTest,
       NormalProblemPowerSeriesExpansionPreconditioner) {
  SetUpProblem(5);
  EXPECT_TRUE(TestSolver(nullptr, SCHUR_POWER_SERIES_EXPANSION, false));
  EXPECT_TRUE(TestSolver(D_.get(), SCHUR_POWER_SERIES_EXPANSION, false));
}

TEST_F(IterativeSchurComplementSolverTest, ProblemWithNoFBlocks) {
  SetUpProblem(3);
  EXPECT_TRUE(TestSolver(nullptr, SCHUR_JACOBI, false));
  EXPECT_TRUE(TestSolver(D_.get(), SCHUR_JACOBI, false));
}

}  // namespace internal
}  // namespace ceres
