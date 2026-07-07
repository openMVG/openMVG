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

#include "ceres/solver.h"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "ceres/autodiff_cost_function.h"
#include "ceres/evaluation_callback.h"
#include "ceres/manifold.h"
#include "ceres/problem.h"
#include "ceres/problem_impl.h"
#include "ceres/sized_cost_function.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(SolverOptions, DefaultTrustRegionOptionsAreValid) {
  Solver::Options options;
  options.minimizer_type = TRUST_REGION;
  std::string error;
  EXPECT_TRUE(options.IsValid(&error)) << error;
}

TEST(SolverOptions, DefaultLineSearchOptionsAreValid) {
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  std::string error;
  EXPECT_TRUE(options.IsValid(&error)) << error;
}

struct QuadraticCostFunctor {
  template <typename T>
  bool operator()(const T* const x, T* residual) const {
    residual[0] = T(5.0) - *x;
    return true;
  }

  static CostFunction* Create() {
    return new AutoDiffCostFunction<QuadraticCostFunctor, 1, 1>(
        new QuadraticCostFunctor);
  }
};

struct RememberingCallback : public IterationCallback {
  explicit RememberingCallback(double* x) : calls(0), x(x) {}
  CallbackReturnType operator()(const IterationSummary& summary) final {
    x_values.push_back(*x);
    return SOLVER_CONTINUE;
  }
  int calls;
  double* x;
  std::vector<double> x_values;
};

struct NoOpEvaluationCallback : EvaluationCallback {
  void PrepareForEvaluation(bool evaluate_jacobians,
                            bool new_evaluation_point) final {
    (void)evaluate_jacobians;
    (void)new_evaluation_point;
  }
};

TEST(Solver, UpdateStateEveryIterationOptionNoEvaluationCallback) {
  double x = 50.0;
  const double original_x = x;

  Problem::Options problem_options;
  Problem problem(problem_options);
  problem.AddResidualBlock(QuadraticCostFunctor::Create(), nullptr, &x);

  Solver::Options options;
  options.linear_solver_type = DENSE_QR;

  RememberingCallback callback(&x);
  options.callbacks.push_back(&callback);

  Solver::Summary summary;

  int num_iterations;

  // First: update_state_every_iteration=false, evaluation_callback=nullptr.
  Solve(options, &problem, &summary);
  num_iterations =
      summary.num_successful_steps + summary.num_unsuccessful_steps;
  EXPECT_GT(num_iterations, 1);
  for (double value : callback.x_values) {
    EXPECT_EQ(50.0, value);
  }

  // Second: update_state_every_iteration=true, evaluation_callback=nullptr.
  x = 50.0;
  options.update_state_every_iteration = true;
  callback.x_values.clear();
  Solve(options, &problem, &summary);
  num_iterations =
      summary.num_successful_steps + summary.num_unsuccessful_steps;
  EXPECT_GT(num_iterations, 1);
  EXPECT_EQ(original_x, callback.x_values[0]);
  EXPECT_NE(original_x, callback.x_values[1]);
}

TEST(Solver, UpdateStateEveryIterationOptionWithEvaluationCallback) {
  double x = 50.0;
  const double original_x = x;

  Problem::Options problem_options;
  NoOpEvaluationCallback evaluation_callback;
  problem_options.evaluation_callback = &evaluation_callback;

  Problem problem(problem_options);
  problem.AddResidualBlock(QuadraticCostFunctor::Create(), nullptr, &x);

  Solver::Options options;
  options.linear_solver_type = DENSE_QR;
  RememberingCallback callback(&x);
  options.callbacks.push_back(&callback);

  Solver::Summary summary;

  int num_iterations;

  // First: update_state_every_iteration=true, evaluation_callback=!nullptr.
  x = 50.0;
  options.update_state_every_iteration = true;
  callback.x_values.clear();
  Solve(options, &problem, &summary);
  num_iterations =
      summary.num_successful_steps + summary.num_unsuccessful_steps;
  EXPECT_GT(num_iterations, 1);
  EXPECT_EQ(original_x, callback.x_values[0]);
  EXPECT_NE(original_x, callback.x_values[1]);

  // Second: update_state_every_iteration=false, evaluation_callback=!nullptr.
  x = 50.0;
  options.update_state_every_iteration = false;
  callback.x_values.clear();
  Solve(options, &problem, &summary);
  num_iterations =
      summary.num_successful_steps + summary.num_unsuccessful_steps;
  EXPECT_GT(num_iterations, 1);
  EXPECT_EQ(original_x, callback.x_values[0]);
  EXPECT_NE(original_x, callback.x_values[1]);
}

TEST(Solver, CantMixEvaluationCallbackWithInnerIterations) {
  double x = 50.0;
  double y = 60.0;

  Problem::Options problem_options;
  NoOpEvaluationCallback evaluation_callback;
  problem_options.evaluation_callback = &evaluation_callback;

  Problem problem(problem_options);
  problem.AddResidualBlock(QuadraticCostFunctor::Create(), nullptr, &x);
  problem.AddResidualBlock(QuadraticCostFunctor::Create(), nullptr, &y);

  Solver::Options options;
  options.use_inner_iterations = true;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, FAILURE);

  options.use_inner_iterations = false;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
}

// The parameters must be in separate blocks so that they can be individually
// set constant or not.
struct Quadratic4DCostFunction {
  template <typename T>
  bool operator()(const T* const x,
                  const T* const y,
                  const T* const z,
                  const T* const w,
                  T* residual) const {
    // A 4-dimension axis-aligned quadratic.
    residual[0] = T(10.0) - *x + T(20.0) - *y + T(30.0) - *z + T(40.0) - *w;
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
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    residuals[0] = parameters[0][0];
    if (jacobians != nullptr && jacobians[0] != nullptr) {
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
  problem.AddResidualBlock(new UnaryIdentityCostFunction, nullptr, &x);
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
  problem.AddResidualBlock(new UnaryIdentityCostFunction, nullptr, &x);
  problem.SetParameterBlockConstant(&x);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_EQ(summary.initial_cost, 1.0 / 2.0);
  EXPECT_EQ(summary.final_cost, 1.0 / 2.0);
}

template <int kNumResiduals, int... Ns>
class DummyCostFunction : public SizedCostFunction<kNumResiduals, Ns...> {
 public:
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    for (int i = 0; i < kNumResiduals; ++i) {
      residuals[i] = kNumResiduals * kNumResiduals + i;
    }

    return true;
  }
};

TEST(Solver, FixedCostForConstantProblem) {
  double x = 1.0;
  Problem problem;
  problem.AddResidualBlock(new DummyCostFunction<2, 1>(), nullptr, &x);
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

struct LinearCostFunction {
  template <typename T>
  bool operator()(const T* x, const T* y, T* residual) const {
    residual[0] = T(10.0) - *x;
    residual[1] = T(5.0) - *y;
    return true;
  }
  static CostFunction* Create() {
    return new AutoDiffCostFunction<LinearCostFunction, 2, 1, 1>(
        new LinearCostFunction);
  }
};

TEST(Solver, ZeroSizedManifoldHoldsParameterBlockConstant) {
  double x = 0.0;
  double y = 1.0;
  Problem problem;
  problem.AddResidualBlock(LinearCostFunction::Create(), nullptr, &x, &y);
  problem.SetManifold(&y, new SubsetManifold(1, {0}));
  EXPECT_TRUE(problem.IsParameterBlockConstant(&y));

  Solver::Options options;
  options.function_tolerance = 0.0;
  options.gradient_tolerance = 0.0;
  options.parameter_tolerance = 0.0;
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  EXPECT_EQ(summary.termination_type, CONVERGENCE);
  EXPECT_NEAR(x, 10.0, 1e-7);
  EXPECT_EQ(y, 1.0);
}

TEST(Solver, DenseNormalCholeskyOptions) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = DENSE_NORMAL_CHOLESKY;
  EXPECT_TRUE(options.IsValid(&message));

  options.dense_linear_algebra_library_type = EIGEN;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.use_mixed_precision_solves = true;
  EXPECT_TRUE(options.IsValid(&message));

  if (IsDenseLinearAlgebraLibraryTypeAvailable(LAPACK)) {
    options.use_mixed_precision_solves = false;
    options.dense_linear_algebra_library_type = LAPACK;

    EXPECT_TRUE(options.IsValid(&message));
    options.use_mixed_precision_solves = true;
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    options.use_mixed_precision_solves = false;
    options.dense_linear_algebra_library_type = LAPACK;
    EXPECT_FALSE(options.IsValid(&message));
  }
}

TEST(Solver, DenseQrOptions) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = DENSE_QR;

  options.use_mixed_precision_solves = false;
  options.dense_linear_algebra_library_type = EIGEN;
  EXPECT_TRUE(options.IsValid(&message));

  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  if (IsDenseLinearAlgebraLibraryTypeAvailable(LAPACK)) {
    options.use_mixed_precision_solves = false;
    options.dense_linear_algebra_library_type = LAPACK;
    EXPECT_TRUE(options.IsValid(&message));
    options.use_mixed_precision_solves = true;
    EXPECT_FALSE(options.IsValid(&message));
  } else {
    options.use_mixed_precision_solves = false;
    options.dense_linear_algebra_library_type = LAPACK;
    EXPECT_FALSE(options.IsValid(&message));
  }
}

TEST(Solver, SparseNormalCholeskyOptionsNoSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, SparseNormalCholeskyOptionsEigenSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.linear_solver_ordering_type = AMD;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    EXPECT_FALSE(options.IsValid(&message));
  }

  if (IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_TRUE(options.IsValid(&message));
  }

#ifndef CERES_NO_EIGEN_METIS
  options.linear_solver_ordering_type = NESDIS;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_TRUE(options.IsValid(&message));
  }
#else
  options.linear_solver_ordering_type = NESDIS;
  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  EXPECT_FALSE(options.IsValid(&message));
#endif
}

TEST(Solver, SparseNormalCholeskyOptionsSuiteSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.linear_solver_ordering_type = AMD;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    EXPECT_FALSE(options.IsValid(&message));
  }

  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

#ifndef CERES_NO_CHOLMOD_PARTITION
  options.linear_solver_ordering_type = NESDIS;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }
#else
  options.linear_solver_ordering_type = NESDIS;
  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  EXPECT_FALSE(options.IsValid(&message));
#endif
}

TEST(Solver, SparseNormalCholeskyOptionsAccelerateSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.sparse_linear_algebra_library_type = ACCELERATE_SPARSE;
  options.linear_solver_ordering_type = AMD;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    EXPECT_FALSE(options.IsValid(&message));
  }

  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

  options.linear_solver_ordering_type = NESDIS;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }
}

TEST(Solver, DenseSchurOptions) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = DENSE_SCHUR;
  options.dense_linear_algebra_library_type = EIGEN;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.use_mixed_precision_solves = true;
  options.dynamic_sparsity = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.use_mixed_precision_solves = true;
  options.dynamic_sparsity = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.dense_linear_algebra_library_type = LAPACK;
  if (IsDenseLinearAlgebraLibraryTypeAvailable(
          options.dense_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }
}

TEST(Solver, SparseSchurOptionsNoSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_SCHUR;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, SparseSchurOptionsEigenSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_SCHUR;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.linear_solver_ordering_type = AMD;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    EXPECT_FALSE(options.IsValid(&message));
  }

  if (IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

#ifndef CERES_NO_EIGEN_METIS
  options.linear_solver_ordering_type = NESDIS;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(EIGEN_SPARSE)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }
#else
  options.linear_solver_ordering_type = NESDIS;
  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  EXPECT_FALSE(options.IsValid(&message));
#endif
}

TEST(Solver, SparseSchurOptionsSuiteSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_SCHUR;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.linear_solver_ordering_type = AMD;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    EXPECT_FALSE(options.IsValid(&message));
  }

  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

#ifndef CERES_NO_CHOLMOD_PARTITION
  options.linear_solver_ordering_type = NESDIS;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }
#else
  options.linear_solver_ordering_type = NESDIS;
  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  EXPECT_FALSE(options.IsValid(&message));
#endif
}

TEST(Solver, SparseSchurOptionsAccelerateSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = SPARSE_SCHUR;
  options.sparse_linear_algebra_library_type = ACCELERATE_SPARSE;
  options.linear_solver_ordering_type = AMD;

  options.use_mixed_precision_solves = false;
  options.dynamic_sparsity = false;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    EXPECT_TRUE(options.IsValid(&message));
  } else {
    EXPECT_FALSE(options.IsValid(&message));
  }

  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

  options.linear_solver_ordering_type = NESDIS;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.use_mixed_precision_solves = false;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));

    options.use_mixed_precision_solves = true;
    options.dynamic_sparsity = true;
    EXPECT_FALSE(options.IsValid(&message));
  }
}

TEST(Solver, CgnrOptionsIdentityPreconditioner) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = CGNR;
  options.preconditioner_type = IDENTITY;
  options.sparse_linear_algebra_library_type = NO_SPARSE;

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = ACCELERATE_SPARSE;
  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = CUDA_SPARSE;
  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(CUDA_SPARSE));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, CgnrOptionsJacobiPreconditioner) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = CGNR;
  options.preconditioner_type = JACOBI;
  options.sparse_linear_algebra_library_type = NO_SPARSE;

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = SUITE_SPARSE;

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = ACCELERATE_SPARSE;

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_TRUE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = CUDA_SPARSE;

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(CUDA_SPARSE));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, CgnrOptionsSubsetPreconditioner) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = CGNR;
  options.preconditioner_type = SUBSET;

  options.sparse_linear_algebra_library_type = NO_SPARSE;
  EXPECT_FALSE(options.IsValid(&message));
  options.residual_blocks_for_subset_preconditioner.insert(nullptr);
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));

  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.dynamic_sparsity = false;
    options.use_mixed_precision_solves = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.dynamic_sparsity = true;
    options.use_mixed_precision_solves = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.dynamic_sparsity = false;
    options.use_mixed_precision_solves = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.dynamic_sparsity = false;
    options.use_mixed_precision_solves = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.dynamic_sparsity = true;
    options.use_mixed_precision_solves = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.dynamic_sparsity = false;
    options.use_mixed_precision_solves = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

  options.sparse_linear_algebra_library_type = ACCELERATE_SPARSE;
  if (IsSparseLinearAlgebraLibraryTypeAvailable(
          options.sparse_linear_algebra_library_type)) {
    options.dynamic_sparsity = false;
    options.use_mixed_precision_solves = false;
    EXPECT_TRUE(options.IsValid(&message));

    options.dynamic_sparsity = true;
    options.use_mixed_precision_solves = false;
    EXPECT_FALSE(options.IsValid(&message));

    options.dynamic_sparsity = false;
    options.use_mixed_precision_solves = true;
    EXPECT_FALSE(options.IsValid(&message));
  }

  options.sparse_linear_algebra_library_type = CUDA_SPARSE;
  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = true;
  options.use_mixed_precision_solves = false;
  EXPECT_FALSE(options.IsValid(&message));

  options.dynamic_sparsity = false;
  options.use_mixed_precision_solves = true;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, CgnrOptionsSchurPreconditioners) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = CGNR;
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeSchurOptionsNoSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = ITERATIVE_SCHUR;
  options.sparse_linear_algebra_library_type = NO_SPARSE;
  options.preconditioner_type = IDENTITY;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = SUBSET;
  EXPECT_FALSE(options.IsValid(&message));

  options.use_explicit_schur_complement = true;
  options.preconditioner_type = IDENTITY;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeSchurOptionsEigenSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = ITERATIVE_SCHUR;
  options.sparse_linear_algebra_library_type = EIGEN_SPARSE;
  options.preconditioner_type = IDENTITY;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(
                options.sparse_linear_algebra_library_type));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(
                options.sparse_linear_algebra_library_type));
  options.preconditioner_type = SUBSET;
  EXPECT_FALSE(options.IsValid(&message));

  options.use_explicit_schur_complement = true;
  options.preconditioner_type = IDENTITY;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeSchurOptionsSuiteSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = ITERATIVE_SCHUR;
  options.sparse_linear_algebra_library_type = SUITE_SPARSE;
  options.preconditioner_type = IDENTITY;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(
                options.sparse_linear_algebra_library_type));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(
                options.sparse_linear_algebra_library_type));
  options.preconditioner_type = SUBSET;
  EXPECT_FALSE(options.IsValid(&message));

  options.use_explicit_schur_complement = true;
  options.preconditioner_type = IDENTITY;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_FALSE(options.IsValid(&message));
}

TEST(Solver, IterativeSchurOptionsAccelerateSparse) {
  std::string message;
  Solver::Options options;
  options.linear_solver_type = ITERATIVE_SCHUR;
  options.sparse_linear_algebra_library_type = ACCELERATE_SPARSE;
  options.preconditioner_type = IDENTITY;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(
                options.sparse_linear_algebra_library_type));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_EQ(options.IsValid(&message),
            IsSparseLinearAlgebraLibraryTypeAvailable(
                options.sparse_linear_algebra_library_type));
  options.preconditioner_type = SUBSET;
  EXPECT_FALSE(options.IsValid(&message));

  options.use_explicit_schur_complement = true;
  options.preconditioner_type = IDENTITY;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = SCHUR_JACOBI;
  EXPECT_TRUE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_JACOBI;
  EXPECT_FALSE(options.IsValid(&message));
  options.preconditioner_type = CLUSTER_TRIDIAGONAL;
  EXPECT_FALSE(options.IsValid(&message));
}

class LargeCostCostFunction : public SizedCostFunction<1, 1> {
 public:
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    residuals[0] = 1e300;
    if (jacobians && jacobians[0]) {
      jacobians[0][0] = 1.0;
    }
    return true;
  }
};

TEST(Solver, LargeCostProblem) {
  double x = 1;
  Problem problem;
  problem.AddResidualBlock(new LargeCostCostFunction, nullptr, &x);
  Solver::Options options;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  LOG(INFO) << summary.FullReport();
  EXPECT_EQ(summary.termination_type, FAILURE);
}

}  // namespace ceres::internal
