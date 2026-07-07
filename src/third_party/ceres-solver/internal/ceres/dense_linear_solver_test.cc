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

#include "ceres/casts.h"
#include "ceres/context_impl.h"
#include "ceres/internal/config.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

using Param = ::testing::
    tuple<LinearSolverType, DenseLinearAlgebraLibraryType, bool, int>;

static std::string ParamInfoToString(testing::TestParamInfo<Param> info) {
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

  std::unique_ptr<LinearLeastSquaresProblem> problem =
      CreateLinearLeastSquaresProblemFromId(testing::get<3>(param));
  DenseSparseMatrix lhs(*down_cast<TripletSparseMatrix*>(problem->A.get()));

  const int num_cols = lhs.num_cols();
  const int num_rows = lhs.num_rows();

  Vector rhs = Vector::Zero(num_rows + num_cols);
  rhs.head(num_rows) = ConstVectorRef(problem->b.get(), num_rows);

  LinearSolver::Options options;
  options.type = ::testing::get<0>(param);
  options.dense_linear_algebra_library_type = ::testing::get<1>(param);
  ContextImpl context;
  options.context = &context;
  std::unique_ptr<LinearSolver> solver(LinearSolver::Create(options));

  LinearSolver::PerSolveOptions per_solve_options;
  if (regularized) {
    per_solve_options.D = problem->D.get();
  }

  Vector solution(num_cols);
  LinearSolver::Summary summary =
      solver->Solve(&lhs, rhs.data(), per_solve_options, solution.data());
  EXPECT_EQ(summary.termination_type, LinearSolverTerminationType::SUCCESS);

  Vector normal_rhs = lhs.matrix().transpose() * rhs.head(num_rows);
  Matrix normal_lhs = lhs.matrix().transpose() * lhs.matrix();

  if (regularized) {
    ConstVectorRef diagonal(problem->D.get(), num_cols);
    normal_lhs += diagonal.array().square().matrix().asDiagonal();
  }

  Vector actual_normal_rhs = normal_lhs * solution;

  const double normalized_residual =
      (normal_rhs - actual_normal_rhs).norm() / normal_rhs.norm();

  EXPECT_NEAR(
      normalized_residual, 0.0, 10 * std::numeric_limits<double>::epsilon())
      << "\nexpected: " << normal_rhs.transpose()
      << "\nactual: " << actual_normal_rhs.transpose();
}

namespace {

// TODO(sameeragarwal): Should we move away from hard coded linear
// least squares problem to randomly generated ones?
#ifndef CERES_NO_LAPACK

INSTANTIATE_TEST_SUITE_P(
    DenseLinearSolver,
    DenseLinearSolverTest,
    ::testing::Combine(::testing::Values(DENSE_QR, DENSE_NORMAL_CHOLESKY),
                       ::testing::Values(EIGEN, LAPACK),
                       ::testing::Values(true, false),
                       ::testing::Values(0, 1)),
    ParamInfoToString);

#else

INSTANTIATE_TEST_SUITE_P(
    DenseLinearSolver,
    DenseLinearSolverTest,
    ::testing::Combine(::testing::Values(DENSE_QR, DENSE_NORMAL_CHOLESKY),
                       ::testing::Values(EIGEN),
                       ::testing::Values(true, false),
                       ::testing::Values(0, 1)),
    ParamInfoToString);

#endif
}  // namespace
}  // namespace ceres::internal
