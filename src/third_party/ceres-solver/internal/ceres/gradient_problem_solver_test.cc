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
// Author: strandmark@google.com (Petter Strandmark)

#include "ceres/gradient_problem_solver.h"

#include "ceres/gradient_problem.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// Rosenbrock function; see http://en.wikipedia.org/wiki/Rosenbrock_function .
class Rosenbrock : public ceres::FirstOrderFunction {
 public:
  bool Evaluate(const double* parameters,
                double* cost,
                double* gradient) const final {
    const double x = parameters[0];
    const double y = parameters[1];

    cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
    if (gradient != nullptr) {
      gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
      gradient[1] = 200.0 * (y - x * x);
    }
    return true;
  }

  int NumParameters() const final { return 2; }
};

TEST(GradientProblemSolver, SolvesRosenbrockWithDefaultOptions) {
  const double expected_tolerance = 1e-9;
  double parameters[2] = {-1.2, 0.0};

  ceres::GradientProblemSolver::Options options;
  ceres::GradientProblemSolver::Summary summary;
  ceres::GradientProblem problem(new Rosenbrock());
  ceres::Solve(options, problem, parameters, &summary);

  EXPECT_EQ(CONVERGENCE, summary.termination_type);
  EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
  EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
}

class QuadraticFunction : public ceres::FirstOrderFunction {
  bool Evaluate(const double* parameters,
                double* cost,
                double* gradient) const final {
    const double x = parameters[0];
    *cost = 0.5 * (5.0 - x) * (5.0 - x);
    if (gradient != nullptr) {
      gradient[0] = x - 5.0;
    }

    return true;
  }
  int NumParameters() const final { return 1; }
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

TEST(Solver, UpdateStateEveryIterationOption) {
  double x = 50.0;
  const double original_x = x;

  ceres::GradientProblem problem(new QuadraticFunction);
  ceres::GradientProblemSolver::Options options;
  RememberingCallback callback(&x);
  options.callbacks.push_back(&callback);
  ceres::GradientProblemSolver::Summary summary;

  int num_iterations;

  // First try: no updating.
  ceres::Solve(options, problem, &x, &summary);
  num_iterations = summary.iterations.size() - 1;
  EXPECT_GT(num_iterations, 1);
  for (double value : callback.x_values) {
    EXPECT_EQ(50.0, value);
  }

  // Second try: with updating
  x = 50.0;
  options.update_state_every_iteration = true;
  callback.x_values.clear();
  ceres::Solve(options, problem, &x, &summary);
  num_iterations = summary.iterations.size() - 1;
  EXPECT_GT(num_iterations, 1);
  EXPECT_EQ(original_x, callback.x_values[0]);
  EXPECT_NE(original_x, callback.x_values[1]);
}

}  // namespace ceres::internal
