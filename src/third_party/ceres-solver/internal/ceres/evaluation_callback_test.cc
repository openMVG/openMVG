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
// Author: mierle@gmail.com (Keir Mierle)

#include "ceres/evaluation_callback.h"

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "ceres/autodiff_cost_function.h"
#include "ceres/problem.h"
#include "ceres/problem_impl.h"
#include "ceres/sized_cost_function.h"
#include "ceres/solver.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// Use an inline hash function to avoid portability wrangling. Algorithm from
// Daniel Bernstein, known as the "djb2" hash.
template <typename T>
uint64_t Djb2Hash(const T* data, const int size) {
  uint64_t hash = 5381;
  const auto* data_as_bytes = reinterpret_cast<const uint8_t*>(data);
  for (int i = 0; i < sizeof(*data) * size; ++i) {
    hash = hash * 33 + data_as_bytes[i];
  }
  return hash;
}

const double kUninitialized = 0;

// Generally multiple inheritance is a terrible idea, but in this (test)
// case it makes for a relatively elegant test implementation.
struct WigglyBowlCostFunctionAndEvaluationCallback : SizedCostFunction<2, 2>,
                                                     EvaluationCallback {
  explicit WigglyBowlCostFunctionAndEvaluationCallback(double* parameter)
      : EvaluationCallback(),
        user_parameter_block(parameter),
        prepare_num_calls(0),
        prepare_requested_jacobians(false),
        prepare_new_evaluation_point(false),
        prepare_parameter_hash(kUninitialized),
        evaluate_num_calls(0),
        evaluate_last_parameter_hash(kUninitialized) {}

  // Evaluation callback interface. This checks that all the preconditions are
  // met at the point that Ceres calls into it.
  void PrepareForEvaluation(bool evaluate_jacobians,
                            bool new_evaluation_point) final {
    // At this point, the incoming parameters are implicitly pushed by Ceres
    // into the user parameter blocks; in contrast to in Evaluate().
    uint64_t incoming_parameter_hash = Djb2Hash(user_parameter_block, 2);

    // Check: Prepare() & Evaluate() come in pairs, in that order. Before this
    // call, the number of calls excluding this one should match.
    EXPECT_EQ(prepare_num_calls, evaluate_num_calls);

    // Check: new_evaluation_point indicates that the parameter has changed.
    if (new_evaluation_point) {
      // If it's a new evaluation point, then the parameter should have
      // changed. Technically, it's not required that it must change but
      // in practice it does, and that helps with testing.
      EXPECT_NE(evaluate_last_parameter_hash, incoming_parameter_hash);
      EXPECT_NE(prepare_parameter_hash, incoming_parameter_hash);
    } else {
      // If this is the same evaluation point as last time, ensure that
      // the parameters match both from the previous evaluate, the
      // previous prepare, and the current prepare.
      EXPECT_EQ(evaluate_last_parameter_hash, prepare_parameter_hash);
      EXPECT_EQ(evaluate_last_parameter_hash, incoming_parameter_hash);
    }

    // Save details for to check at the next call to Evaluate().
    prepare_num_calls++;
    prepare_requested_jacobians = evaluate_jacobians;
    prepare_new_evaluation_point = new_evaluation_point;
    prepare_parameter_hash = incoming_parameter_hash;
  }

  // Cost function interface. This checks that preconditions that were
  // set as part of the PrepareForEvaluation() call are met in this one.
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    // Cost function implementation of the "Wiggly Bowl" function:
    //
    //   1/2 * [(y - a*sin(x))^2 + x^2],
    //
    // expressed as a Ceres cost function with two residuals:
    //
    //   r[0] = y - a*sin(x)
    //   r[1] = x.
    //
    // This is harder to optimize than the Rosenbrock function because the
    // minimizer has to navigate a sine-shaped valley while descending the 1D
    // parabola formed along the y axis. Note that the "a" needs to be more
    // than 5 to get a strong enough wiggle effect in the cost surface to
    // trigger failed iterations in the optimizer.
    const double a = 10.0;
    double x = (*parameters)[0];
    double y = (*parameters)[1];
    residuals[0] = y - a * sin(x);
    residuals[1] = x;
    if (jacobians != nullptr) {
      (*jacobians)[2 * 0 + 0] = -a * cos(x);  // df1/dx
      (*jacobians)[2 * 0 + 1] = 1.0;          // df1/dy
      (*jacobians)[2 * 1 + 0] = 1.0;          // df2/dx
      (*jacobians)[2 * 1 + 1] = 0.0;          // df2/dy
    }

    uint64_t incoming_parameter_hash = Djb2Hash(*parameters, 2);

    // Check: PrepareForEvaluation() & Evaluate() come in pairs, in that order.
    EXPECT_EQ(prepare_num_calls, evaluate_num_calls + 1);

    // Check: if new_evaluation_point indicates that the parameter has
    // changed, it has changed; otherwise it is the same.
    if (prepare_new_evaluation_point) {
      EXPECT_NE(evaluate_last_parameter_hash, incoming_parameter_hash);
    } else {
      EXPECT_NE(evaluate_last_parameter_hash, kUninitialized);
      EXPECT_EQ(evaluate_last_parameter_hash, incoming_parameter_hash);
    }

    // Check: Parameter matches value in in parameter blocks during prepare.
    EXPECT_EQ(prepare_parameter_hash, incoming_parameter_hash);

    // Check: jacobians are requested if they were in PrepareForEvaluation().
    EXPECT_EQ(prepare_requested_jacobians, jacobians != nullptr);

    evaluate_num_calls++;
    evaluate_last_parameter_hash = incoming_parameter_hash;
    return true;
  }

  // Pointer to the parameter block associated with this cost function.
  // Contents should get set by Ceres before calls to PrepareForEvaluation()
  // and Evaluate().
  double* user_parameter_block;

  // Track state: PrepareForEvaluation().
  //
  // These track details from the PrepareForEvaluation() call (hence the
  // "prepare_" prefix), which are checked for consistency in Evaluate().
  int prepare_num_calls;
  bool prepare_requested_jacobians;
  bool prepare_new_evaluation_point;
  uint64_t prepare_parameter_hash;

  // Track state: Evaluate().
  //
  // These track details from the Evaluate() call (hence the "evaluate_"
  // prefix), which are then checked for consistency in the calls to
  // PrepareForEvaluation(). Mutable is reasonable for this case.
  mutable int evaluate_num_calls;
  mutable uint64_t evaluate_last_parameter_hash;
};

TEST(EvaluationCallback, WithTrustRegionMinimizer) {
  double parameters[2] = {50.0, 50.0};
  const uint64_t original_parameters_hash = Djb2Hash(parameters, 2);

  WigglyBowlCostFunctionAndEvaluationCallback cost_function(parameters);
  Problem::Options problem_options;
  problem_options.evaluation_callback = &cost_function;
  problem_options.cost_function_ownership = DO_NOT_TAKE_OWNERSHIP;
  Problem problem(problem_options);
  problem.AddResidualBlock(&cost_function, nullptr, parameters);

  Solver::Options options;
  options.linear_solver_type = DENSE_QR;
  options.max_num_iterations = 50;

  // Run the solve. Checking is done inside the cost function / callback.
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Ensure that this was a hard cost function (not all steps succeed).
  EXPECT_GT(summary.num_successful_steps, 10);
  EXPECT_GT(summary.num_unsuccessful_steps, 10);

  // Ensure PrepareForEvaluation() is called the appropriate number of times.
  EXPECT_EQ(
      cost_function.prepare_num_calls,
      // Unsuccessful steps are evaluated only once (no jacobians).
      summary.num_unsuccessful_steps +
          // Successful steps are evaluated twice: with and without jacobians.
          2 * summary.num_successful_steps
          // Final iteration doesn't re-evaluate the jacobian.
          // Note: This may be sensitive to tweaks to the TR algorithm; if
          // this becomes too brittle, remove this EXPECT_EQ() entirely.
          - 1);

  // Ensure the callback calls ran a reasonable number of times.
  EXPECT_GT(cost_function.prepare_num_calls, 0);
  EXPECT_GT(cost_function.evaluate_num_calls, 0);
  EXPECT_EQ(cost_function.prepare_num_calls, cost_function.evaluate_num_calls);

  // Ensure that the parameters did actually change.
  EXPECT_NE(Djb2Hash(parameters, 2), original_parameters_hash);
}

// r = 1 - x
struct LinearResidual {
  template <typename T>
  bool operator()(const T* x, T* residuals) const {
    residuals[0] = 1.0 - x[0];
    return true;
  }

  static CostFunction* Create() {
    return new AutoDiffCostFunction<LinearResidual, 1, 1>(new LinearResidual);
  };
};

// Increments a counter everytime PrepareForEvaluation is called.
class IncrementingEvaluationCallback : public EvaluationCallback {
 public:
  void PrepareForEvaluation(bool evaluate_jacobians,
                            bool new_evaluation_point) final {
    (void)evaluate_jacobians;
    (void)new_evaluation_point;
    counter_ += 1.0;
  }

  double counter() const { return counter_; }

 private:
  double counter_ = -1;
};

// r = IncrementingEvaluationCallback::counter - x
struct EvaluationCallbackResidual {
  explicit EvaluationCallbackResidual(
      const IncrementingEvaluationCallback& callback)
      : callback(callback) {}

  template <typename T>
  bool operator()(const T* x, T* residuals) const {
    residuals[0] = callback.counter() - x[0];
    return true;
  }

  const IncrementingEvaluationCallback& callback;

  static CostFunction* Create(IncrementingEvaluationCallback& callback) {
    return new AutoDiffCostFunction<EvaluationCallbackResidual, 1, 1>(
        new EvaluationCallbackResidual(callback));
  };
};

// The following test, constructs a problem with residual blocks all
// of whose parameters are constant, so they are evaluated once
// outside the Minimizer to compute Solver::Summary::fixed_cost.
//
// The cost function for this residual block depends on the
// IncrementingEvaluationCallback::counter_, by checking the value of
// the fixed cost, we can check if the IncrementingEvaluationCallback
// was called.
TEST(EvaluationCallback, EvaluationCallbackIsCalledBeforeFixedCostIsEvaluated) {
  double x = 1;
  double y = 2;
  std::unique_ptr<IncrementingEvaluationCallback> callback(
      new IncrementingEvaluationCallback);
  Problem::Options problem_options;
  problem_options.evaluation_callback = callback.get();
  Problem problem(problem_options);
  problem.AddResidualBlock(LinearResidual::Create(), nullptr, &x);
  problem.AddResidualBlock(
      EvaluationCallbackResidual::Create(*callback), nullptr, &y);
  problem.SetParameterBlockConstant(&y);

  Solver::Options options;
  options.linear_solver_type = DENSE_QR;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  EXPECT_EQ(summary.fixed_cost, 2.0);
  EXPECT_EQ(summary.final_cost, summary.fixed_cost);
  EXPECT_GT(callback->counter(), 0);
}

static void WithLineSearchMinimizerImpl(
    LineSearchType line_search,
    LineSearchDirectionType line_search_direction,
    LineSearchInterpolationType line_search_interpolation) {
  double parameters[2] = {50.0, 50.0};
  const uint64_t original_parameters_hash = Djb2Hash(parameters, 2);

  WigglyBowlCostFunctionAndEvaluationCallback cost_function(parameters);
  Problem::Options problem_options;
  problem_options.evaluation_callback = &cost_function;
  problem_options.cost_function_ownership = DO_NOT_TAKE_OWNERSHIP;
  Problem problem(problem_options);
  problem.AddResidualBlock(&cost_function, nullptr, parameters);

  Solver::Options options;
  options.linear_solver_type = DENSE_QR;
  options.max_num_iterations = 50;
  options.minimizer_type = ceres::LINE_SEARCH;

  options.line_search_type = line_search;
  options.line_search_direction_type = line_search_direction;
  options.line_search_interpolation_type = line_search_interpolation;

  // Run the solve. Checking is done inside the cost function / callback.
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  // Ensure the callback calls ran a reasonable number of times.
  EXPECT_GT(summary.num_line_search_steps, 10);
  EXPECT_GT(cost_function.prepare_num_calls, 30);
  EXPECT_EQ(cost_function.prepare_num_calls, cost_function.evaluate_num_calls);

  // Ensure that the parameters did actually change.
  EXPECT_NE(Djb2Hash(parameters, 2), original_parameters_hash);
}

// Note: These tests omit combinations of Wolfe line search with bisection.
// Due to an implementation quirk in Wolfe line search with bisection, there
// are calls to re-evaluate an existing point with new_point = true. That
// causes the (overly) strict tests to break, since they check the new_point
// preconditions in an if-and-only-if way. Strictly speaking, if new_point =
// true, the interface does not *require* that the point has changed; only that
// if new_point = false, the same point is reused.
//
// Since the strict checking is useful to verify that there aren't missed
// optimizations, omit tests of the Wolfe with bisection cases.

// Wolfe with L-BFGS.
TEST(EvaluationCallback, WithLineSearchMinimizerWolfeLbfgsCubic) {
  WithLineSearchMinimizerImpl(WOLFE, LBFGS, CUBIC);
}
TEST(EvaluationCallback, WithLineSearchMinimizerWolfeLbfgsQuadratic) {
  WithLineSearchMinimizerImpl(WOLFE, LBFGS, QUADRATIC);
}

// Wolfe with full BFGS.
TEST(EvaluationCallback, WithLineSearchMinimizerWolfeBfgsCubic) {
  WithLineSearchMinimizerImpl(WOLFE, BFGS, CUBIC);
}

TEST(EvaluationCallback, WithLineSearchMinimizerWolfeBfgsQuadratic) {
  WithLineSearchMinimizerImpl(WOLFE, BFGS, QUADRATIC);
}

// Armijo with nonlinear conjugate gradient.
TEST(EvaluationCallback, WithLineSearchMinimizerArmijoCubic) {
  WithLineSearchMinimizerImpl(ARMIJO, NONLINEAR_CONJUGATE_GRADIENT, CUBIC);
}

TEST(EvaluationCallback, WithLineSearchMinimizerArmijoBisection) {
  WithLineSearchMinimizerImpl(ARMIJO, NONLINEAR_CONJUGATE_GRADIENT, BISECTION);
}

TEST(EvaluationCallback, WithLineSearchMinimizerArmijoQuadratic) {
  WithLineSearchMinimizerImpl(ARMIJO, NONLINEAR_CONJUGATE_GRADIENT, QUADRATIC);
}

}  // namespace ceres::internal
