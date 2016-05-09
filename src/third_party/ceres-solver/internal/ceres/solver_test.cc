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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/solver.h"

#include <limits>
#include <cmath>
#include <vector>
#include "gtest/gtest.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/autodiff_cost_function.h"
#include "ceres/sized_cost_function.h"
#include "ceres/problem.h"
#include "ceres/problem_impl.h"

namespace ceres {
namespace internal {

using std::string;

TEST(SolverOptions, DefaultTrustRegionOptionsAreValid) {
  Solver::Options options;
  options.minimizer_type = TRUST_REGION;
  string error;
  EXPECT_TRUE(options.IsValid(&error)) << error;
}

TEST(SolverOptions, DefaultLineSearchOptionsAreValid) {
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  string error;
  EXPECT_TRUE(options.IsValid(&error)) << error;
}

struct QuadraticCostFunctor {
  template <typename T> bool operator()(const T* const x,
                                        T* residual) const {
    residual[0] = T(5.0) - *x;
    return true;
  }

  static CostFunction* Create() {
    return new AutoDiffCostFunction<QuadraticCostFunctor, 1, 1>(
        new QuadraticCostFunctor);
  }
};

struct RememberingCallback : public IterationCallback {
  explicit RememberingCallback(double *x) : calls(0), x(x) {}
  virtual ~RememberingCallback() {}
  virtual CallbackReturnType operator()(const IterationSummary& summary) {
    x_values.push_back(*x);
    return SOLVER_CONTINUE;
  }
  int calls;
  double *x;
  std::vector<double> x_values;
};

TEST(Solver, UpdateStateEveryIterationOption) {
  double x = 50.0;
  const double original_x = x;

  scoped_ptr<CostFunction> cost_function(QuadraticCostFunctor::Create());
  Problem::Options problem_options;
  problem_options.cost_function_ownership = DO_NOT_TAKE_OWNERSHIP;
  Problem problem(problem_options);
  problem.AddResidualBlock(cost_function.get(), NULL, &x);

  Solver::Options options;
  options.linear_solver_type = DENSE_QR;

  RememberingCallback callback(&x);
  options.callbacks.push_back(&callback);

  Solver::Summary summary;

  int num_iterations;

  // First try: no updating.
  Solve(options, &problem, &summary);
  num_iterations = summary.num_successful_steps +
                   summary.num_unsuccessful_steps;
  EXPECT_GT(num_iterations, 1);
  for (int i = 0; i < callback.x_values.size(); ++i) {
    EXPECT_EQ(50.0, callback.x_values[i]);
  }

  // Second try: with updating
  x = 50.0;
  options.update_state_every_iteration = true;
  callback.x_values.clear();
  Solve(options, &problem, &summary);
  num_iterations = summary.num_successful_steps +
                   summary.num_unsuccessful_steps;
  EXPECT_GT(num_iterations, 1);
  EXPECT_EQ(original_x, callback.x_values[0]);
  EXPECT_NE(original_x, callback.x_values[1]);
}

// The parameters must be in separate blocks so that they can be individually
// set constant or not.
struct Quadratic4DCostFunction {
  template <typename T> bool operator()(const T* const x,
                                        const T* const y,
                                        const T* const z,
                                        const T* const w,
                                        T* residual) const {
    // A 4-dimension axis-aligned quadratic.
    residual[0] = T(10.0) - *x +
                  T(20.0) - *y +
                  T(30.0) - *z +
                  T(40.0) - *w;
    return true;
  }

  static CostFunction* Create() {
    return new AutoDiffCostFunction<Quadratic4DCostFunction, 1, 1, 1, 1, 1>(
        new Quadratic4DCostFunction);
  }
};

// A cost function that simply returns its argument.
class UnaryIdentityCostFunction : public SizedCostFunction<1, 1> {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    residuals[0] = parameters[0][0];
    if (jacobians != NULL && jacobians[0] != NULL) {
      jacobians[0][0] = 1.0;
    }
    return true;
  }
};

TEST(Solver, TrustRegionProblemHasNoParameterBlocks) {
  Problem problem;
  Solver::Options options;
  options.minimizer_type = TRUST_REGION;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.message,
            "Function tolerance reached. "
            "No non-constant parameter blocks found.");
}

TEST(Solver, LineSearchProblemHasNoParameterBlocks) {
  Problem problem;
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.message,
            "Function tolerance reached. "
            "No non-constant parameter blocks found.");
}

TEST(Solver, TrustRegionProblemHasZeroResiduals) {
  Problem problem;
  double x = 1;
  problem.AddParameterBlock(&x, 1);
  Solver::Options options;
  options.minimizer_type = TRUST_REGION;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.message,
            "Function tolerance reached. "
            "No non-constant parameter blocks found.");
}

TEST(Solver, LineSearchProblemHasZeroResiduals) {
  Problem problem;
  double x = 1;
  problem.AddParameterBlock(&x, 1);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.message,
            "Function tolerance reached. "
            "No non-constant parameter blocks found.");
}

TEST(Solver, TrustRegionProblemIsConstant) {
  Problem problem;
  double x = 1;
  problem.AddResidualBlock(new UnaryIdentityCostFunction, NULL, &x);
  problem.SetParameterBlockConstant(&x);
  Solver::Options options;
  options.minimizer_type = TRUST_REGION;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.initial_cost, 1.0 / 2.0);
  EXPECT_EQ(summary.final_cost, 1.0 / 2.0);
}

TEST(Solver, LineSearchProblemIsConstant) {
  Problem problem;
  double x = 1;
  problem.AddResidualBlock(new UnaryIdentityCostFunction, NULL, &x);
  problem.SetParameterBlockConstant(&x);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.initial_cost, 1.0 / 2.0);
  EXPECT_EQ(summary.final_cost, 1.0 / 2.0);
}

#if defined(CERES_NO_SUITESPARSE)
TEST(Solver, SparseNormalCholeskyNoSuiteSparse) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, SparseSchurNoSuiteSparse) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.linear_solver_type = SPARSE_SCHUR;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}
#endif

#if defined(CERES_NO_CXSPARSE)
TEST(Solver, SparseNormalCholeskyNoCXSparse) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, SparseSchurNoCXSparse) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = CX_SPARSE;
  options.linear_solver_type = SPARSE_SCHUR;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}
#endif

#if !defined(CERES_USE_EIGEN_SPARSE)
TEST(Solver, SparseNormalCholeskyNoEigenSparse) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, SparseSchurNoEigenSparse) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.linear_solver_type = SPARSE_SCHUR;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}
#endif

TEST(Solver, SparseNormalCholeskyNoSparseLibrary) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, SparseSchurNoSparseLibrary) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  options.linear_solver_type = SPARSE_SCHUR;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeSchurWithClusterJacobiPerconditionerNoSparseLibrary) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  options.linear_solver_type = ITERATIVE_SCHUR;
  // Requires SuiteSparse.
  options.preconditioner_type = CLUSTER_JACOBI;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeSchurWithClusterTridiagonalPerconditionerNoSparseLibrary) {
  Solver::Options options;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  options.linear_solver_type = ITERATIVE_SCHUR;
  // Requires SuiteSparse.
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  string message;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeLinearSolverForDogleg) {
  Solver::Options options;
  options.trust_region_strategy_type = DOGLEG;
  string message;
  options.linear_solver_type = ITERATIVE_SCHUR;
  EXPECT_FALSE(options.IsValid(&message));

  options.linear_solver_type = CGNR;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, LinearSolverTypeNormalOperation) {
  Solver::Options options;
  options.linear_solver_type = DENSE_QR;

  string message;
  EXPECT_TRUE(options.IsValid(&message));

  options.linear_solver_type = DENSE_NORMAL_CHOLESKY;
  EXPECT_TRUE(options.IsValid(&message));

  options.linear_solver_type = DENSE_SCHUR;
  EXPECT_TRUE(options.IsValid(&message));

  options.linear_solver_type = SPARSE_SCHUR;
#if defined(CERES_NO_SUITESPARSE) &&            \
    defined(CERES_NO_CXSPARSE) &&               \
   !defined(CERES_USE_EIGEN_SPARSE)
  EXPECT_FALSE(options.IsValid(&message));
#else
  EXPECT_TRUE(options.IsValid(&message));
#endif

  options.linear_solver_type = ITERATIVE_SCHUR;
  EXPECT_TRUE(options.IsValid(&message));
}

template<int kNumResiduals, int N1 = 0, int N2 = 0, int N3 = 0>
class DummyCostFunction : public SizedCostFunction<kNumResiduals, N1, N2, N3> {
 public:
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const {
    for (int i = 0; i < kNumResiduals; ++i) {
      residuals[i] = kNumResiduals * kNumResiduals + i;
    }

    return true;
  }
};

TEST(Solver, FixedCostForConstantProblem) {
  double x = 1.0;
  Problem problem;
  problem.AddResidualBlock(new DummyCostFunction<2, 1>(), NULL, &x);
  problem.SetParameterBlockConstant(&x);
  const double expected_cost = 41.0 / 2.0;  // 1/2 * ((4 + 0)^2 + (4 + 1)^2)
  Solver::Options options;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_TRUE(summary.IsSolutionUsable());
  EXPECT_EQ(summary.fixed_cost, expected_cost);
  EXPECT_EQ(summary.initial_cost, expected_cost);
  EXPECT_EQ(summary.final_cost, expected_cost);
  EXPECT_EQ(summary.iterations.size(), 0);
}

}  // namespace internal
}  // namespace ceres
