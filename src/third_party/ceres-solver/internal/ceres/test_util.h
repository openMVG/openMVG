// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
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
// Author: keir@google.com (Keir Mierle)

#include <string>
#include "ceres/internal/port.h"
#include "ceres/problem.h"
#include "ceres/solver.h"
#include "ceres/stringprintf.h"
#include "gtest/gtest.h"


#ifndef CERES_INTERNAL_TEST_UTIL_H_
#define CERES_INTERNAL_TEST_UTIL_H_

namespace ceres {
namespace internal {

// Expects that x and y have a relative difference of no more than
// max_abs_relative_difference. If either x or y is zero, then the relative
// difference is interpreted as an absolute difference.
bool ExpectClose(double x, double y, double max_abs_relative_difference);

// Expects that for all i = 1,.., n - 1
//
//   |p[i] - q[i]| / max(|p[i]|, |q[i]|) < tolerance
void ExpectArraysClose(int n,
                       const double* p,
                       const double* q,
                       double tolerance);

// Expects that for all i = 1,.., n - 1
//
//   |p[i] / max_norm_p - q[i] / max_norm_q| < tolerance
//
// where max_norm_p and max_norm_q are the max norms of the arrays p
// and q respectively.
void ExpectArraysCloseUptoScale(int n,
                                const double* p,
                                const double* q,
                                double tolerance);

// Construct a fully qualified path for the test file depending on the
// local build/testing environment.
std::string TestFileAbsolutePath(const std::string& filename);

// Struct used for configuring the solver. Used by end-to-end tests.
struct SolverConfig {
  SolverConfig(
      LinearSolverType linear_solver_type,
      SparseLinearAlgebraLibraryType
      sparse_linear_algebra_library_type = NO_SPARSE,
      bool use_automatic_ordering = true,
      PreconditionerType preconditioner_type = IDENTITY,
      int num_threads = 1)
      : linear_solver_type(linear_solver_type),
        sparse_linear_algebra_library_type(sparse_linear_algebra_library_type),
        use_automatic_ordering(use_automatic_ordering),
        preconditioner_type(preconditioner_type),
        num_threads(num_threads) {
  }

  std::string ToString() const {
    return StringPrintf(
        "(%s, %s, %s, %s, %d)",
        LinearSolverTypeToString(linear_solver_type),
        SparseLinearAlgebraLibraryTypeToString(
            sparse_linear_algebra_library_type),
        use_automatic_ordering ? "AUTOMATIC" : "USER",
        PreconditionerTypeToString(preconditioner_type),
        num_threads);
  }

  void UpdateOptions(Solver::Options* options) const {
    options->linear_solver_type = linear_solver_type;
    options->sparse_linear_algebra_library_type =
        sparse_linear_algebra_library_type;
    options->preconditioner_type = preconditioner_type;
    options->num_threads = num_threads;
    options->num_linear_solver_threads = num_threads;

    if (use_automatic_ordering) {
      options->linear_solver_ordering.reset();
    }
  }

  LinearSolverType linear_solver_type;
  SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type;
  bool use_automatic_ordering;
  PreconditionerType preconditioner_type;
  int num_threads;
};

SolverConfig ThreadedSolverConfig(
    LinearSolverType linear_solver_type,
    SparseLinearAlgebraLibraryType
    sparse_linear_algebra_library_type = NO_SPARSE,
    bool use_automatic_ordering = true,
    PreconditionerType preconditioner_type = IDENTITY);

// A templated test fixture, that is used for testing Ceres end to end
// by computing a solution to the problem for a given solver
// configuration and comparing it to a reference solver configuration.
//
// It is assumed that the SystemTestProblem has an Solver::Options
// struct that contains the reference Solver configuration.
template <class SystemTestProblem>
class SystemTest : public::testing::Test {
 protected:
  virtual void SetUp() {
    SystemTestProblem system_test_problem;
    SolveAndEvaluateFinalResiduals(
        *system_test_problem.mutable_solver_options(),
        system_test_problem.mutable_problem(),
        &expected_final_residuals_);
  }

  void RunSolverForConfigAndExpectResidualsMatch(const SolverConfig& config) {
    LOG(INFO) << "Running solver configuration: "
              << config.ToString();

    SystemTestProblem system_test_problem;
    config.UpdateOptions(system_test_problem.mutable_solver_options());
    std::vector<double> final_residuals;
    SolveAndEvaluateFinalResiduals(
        *system_test_problem.mutable_solver_options(),
        system_test_problem.mutable_problem(),
        &final_residuals);

    // We compare solutions by comparing their residual vectors. We do
    // not compare parameter vectors because it is much more brittle
    // and error prone to do so, since the same problem can have
    // nearly the same residuals at two completely different positions
    // in parameter space.
    CHECK_EQ(expected_final_residuals_.size(), final_residuals.size());
    for (int i = 0; i < final_residuals.size(); ++i) {
      EXPECT_NEAR(final_residuals[i],
                  expected_final_residuals_[i],
                  SystemTestProblem::kResidualTolerance)
          << "Not close enough residual:" << i;
    }
  }

  void SolveAndEvaluateFinalResiduals(const Solver::Options& options,
                                      Problem* problem,
                                      std::vector<double>* final_residuals) {
    Solver::Summary summary;
    Solve(options, problem, &summary);
    CHECK_NE(summary.termination_type, ceres::FAILURE);
    problem->Evaluate(Problem::EvaluateOptions(),
                      NULL,
                      final_residuals,
                      NULL,
                      NULL);
  }

  std::vector<double> expected_final_residuals_;
};

}  // namespace internal
}  // namespace ceres

#endif  // CERES_INTERNAL_TEST_UTIL_H_
