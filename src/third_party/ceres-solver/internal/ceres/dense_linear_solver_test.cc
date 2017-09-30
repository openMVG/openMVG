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

#include "ceres/casts.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

typedef ::testing::
    tuple<LinearSolverType, DenseLinearAlgebraLibraryType, bool, int>
        Param;

std::string ParamInfoToString(testing::TestParamInfo<Param> info) {
  Param param = info.param;
  std::stringstream ss;
  ss << LinearSolverTypeToString(::testing::get<0>(param)) << "_"
     << DenseLinearAlgebraLibraryTypeToString(::testing::get<1>(param)) << "_"
     << (::testing::get<2>(param) ? "Regularized" : "Unregularized") << "_"
     << ::testing::get<3>(param);
  return ss.str();
}

class DenseLinearSolverTest : public ::testing::TestWithParam<Param> {};

TEST_P(DenseLinearSolverTest, _) {
  Param param = GetParam();
  const bool regularized = testing::get<2>(param);

  scoped_ptr<LinearLeastSquaresProblem> problem(
      CreateLinearLeastSquaresProblemFromId(testing::get<3>(param)));
  DenseSparseMatrix lhs(*down_cast<TripletSparseMatrix*>(problem->A.get()));

  const int num_cols = lhs.num_cols();
  const int num_rows = lhs.num_rows();

  Vector rhs = Vector::Zero(num_rows + num_cols);
  rhs.head(num_rows) = ConstVectorRef(problem->b.get(), num_rows);

  LinearSolver::Options options;
  options.type = ::testing::get<0>(param);
  options.dense_linear_algebra_library_type = ::testing::get<1>(param);
  scoped_ptr<LinearSolver> solver(LinearSolver::Create(options));

  LinearSolver::PerSolveOptions per_solve_options;
  if (regularized) {
    per_solve_options.D = problem->D.get();
  }

  Vector solution(num_cols);
  LinearSolver::Summary summary =
      solver->Solve(&lhs, rhs.data(), per_solve_options, solution.data());
  EXPECT_EQ(summary.termination_type, LINEAR_SOLVER_SUCCESS);

  // If solving for the regularized solution, add the diagonal to the
  // matrix. This makes subsequent computations simpler.
  if (testing::get<2>(param)) {
    lhs.AppendDiagonal(problem->D.get());
  };

  Vector tmp = Vector::Zero(num_rows + num_cols);
  lhs.RightMultiply(solution.data(), tmp.data());
  Vector actual_normal_rhs = Vector::Zero(num_cols);
  lhs.LeftMultiply(tmp.data(), actual_normal_rhs.data());

  Vector expected_normal_rhs = Vector::Zero(num_cols);
  lhs.LeftMultiply(rhs.data(), expected_normal_rhs.data());
  const double residual = (expected_normal_rhs - actual_normal_rhs).norm() /
                          expected_normal_rhs.norm();

  EXPECT_NEAR(residual, 0.0, 10 * std::numeric_limits<double>::epsilon());
}

// TODO(sameeragarwal): Should we move away from hard coded linear
// least squares problem to randomly generated ones?
#ifndef CERES_NO_LAPACK

INSTANTIATE_TEST_CASE_P(
    DenseLinearSolver,
    DenseLinearSolverTest,
    ::testing::Combine(::testing::Values(DENSE_QR, DENSE_NORMAL_CHOLESKY),
                       ::testing::Values(EIGEN, LAPACK),
                       ::testing::Values(true, false),
                       ::testing::Values(0, 1)),
    ParamInfoToString);

#else

INSTANTIATE_TEST_CASE_P(
    DenseLinearSolver,
    DenseLinearSolverTest,
    ::testing::Combine(::testing::Values(DENSE_QR, DENSE_NORMAL_CHOLESKY),
                       ::testing::Values(EIGEN),
                       ::testing::Values(true, false),
                       ::testing::Values(0, 1)),
    ParamInfoToString);

#endif

}  // namespace internal
}  // namespace ceres
