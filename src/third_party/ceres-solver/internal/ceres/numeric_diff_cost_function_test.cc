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
//         tbennun@gmail.com (Tal Ben-Nun)

#include "ceres/numeric_diff_cost_function.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "ceres/internal/macros.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/array_utils.h"
#include "ceres/numeric_diff_test_utils.h"
#include "ceres/test_util.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(NumericDiffCostFunction, EasyCaseFunctorCentralDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor,
                                  CENTRAL,
                                  3,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new EasyFunctor));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, CENTRAL);
}

TEST(NumericDiffCostFunction, EasyCaseFunctorForwardDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor,
                                  FORWARD,
                                  3,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new EasyFunctor));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, FORWARD);
}

TEST(NumericDiffCostFunction, EasyCaseFunctorRidders) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor,
                                  RIDDERS,
                                  3,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new EasyFunctor));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, RIDDERS);
}

TEST(NumericDiffCostFunction, EasyCaseCostFunctionCentralDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyCostFunction,
                                  CENTRAL,
                                  3,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new EasyCostFunction, TAKE_OWNERSHIP));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, CENTRAL);
}

TEST(NumericDiffCostFunction, EasyCaseCostFunctionForwardDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyCostFunction,
                                  FORWARD,
                                  3,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new EasyCostFunction, TAKE_OWNERSHIP));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, FORWARD);
}

TEST(NumericDiffCostFunction, EasyCaseCostFunctionRidders) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyCostFunction,
                                  RIDDERS,
                                  3,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new EasyCostFunction, TAKE_OWNERSHIP));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, RIDDERS);
}

TEST(NumericDiffCostFunction,
     TranscendentalCaseFunctorCentralDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<TranscendentalFunctor,
                                  CENTRAL,
                                  2,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new TranscendentalFunctor));
  TranscendentalFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, CENTRAL);
}

TEST(NumericDiffCostFunction,
     TranscendentalCaseFunctorForwardDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<TranscendentalFunctor,
                                  FORWARD,
                                  2,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new TranscendentalFunctor));
  TranscendentalFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, FORWARD);
}

TEST(NumericDiffCostFunction,
     TranscendentalCaseFunctorRidders) {
  NumericDiffOptions options;

  // Using a smaller initial step size to overcome oscillatory function
  // behavior.
  options.ridders_relative_initial_step_size = 1e-3;

  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<TranscendentalFunctor,
                                  RIDDERS,
                                  2,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new TranscendentalFunctor, TAKE_OWNERSHIP, 2, options));
  TranscendentalFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, RIDDERS);
}

TEST(NumericDiffCostFunction,
     TranscendentalCaseCostFunctionCentralDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<TranscendentalCostFunction,
                                  CENTRAL,
                                  2,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new TranscendentalCostFunction, TAKE_OWNERSHIP));
  TranscendentalFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, CENTRAL);
}

TEST(NumericDiffCostFunction,
     TranscendentalCaseCostFunctionForwardDifferences) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<TranscendentalCostFunction,
                                  FORWARD,
                                  2,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new TranscendentalCostFunction, TAKE_OWNERSHIP));
  TranscendentalFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, FORWARD);
}

TEST(NumericDiffCostFunction,
     TranscendentalCaseCostFunctionRidders) {
  NumericDiffOptions options;

  // Using a smaller initial step size to overcome oscillatory function
  // behavior.
  options.ridders_relative_initial_step_size = 1e-3;

  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<TranscendentalCostFunction,
                                  RIDDERS,
                                  2,  /* number of residuals */
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
          new TranscendentalCostFunction, TAKE_OWNERSHIP, 2, options));
  TranscendentalFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, RIDDERS);
}

template<int num_rows, int num_cols>
class SizeTestingCostFunction : public SizedCostFunction<num_rows, num_cols> {
 public:
  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    return true;
  }
};

// As described in
// http://forum.kde.org/viewtopic.php?f=74&t=98536#p210774
// Eigen3 has restrictions on the Row/Column major storage of vectors,
// depending on their dimensions. This test ensures that the correct
// templates are instantiated for various shapes of the Jacobian
// matrix.
TEST(NumericDiffCostFunction, EigenRowMajorColMajorTest) {
  scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<SizeTestingCostFunction<1,1>,  CENTRAL, 1, 1>(
          new SizeTestingCostFunction<1,1>, ceres::TAKE_OWNERSHIP));

  cost_function.reset(
      new NumericDiffCostFunction<SizeTestingCostFunction<2,1>,  CENTRAL, 2, 1>(
          new SizeTestingCostFunction<2,1>, ceres::TAKE_OWNERSHIP));

  cost_function.reset(
      new NumericDiffCostFunction<SizeTestingCostFunction<1,2>,  CENTRAL, 1, 2>(
          new SizeTestingCostFunction<1,2>, ceres::TAKE_OWNERSHIP));

  cost_function.reset(
      new NumericDiffCostFunction<SizeTestingCostFunction<2,2>,  CENTRAL, 2, 2>(
          new SizeTestingCostFunction<2,2>, ceres::TAKE_OWNERSHIP));

  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor, CENTRAL, ceres::DYNAMIC, 1, 1>(
          new EasyFunctor, TAKE_OWNERSHIP, 1));

  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor, CENTRAL, ceres::DYNAMIC, 1, 1>(
          new EasyFunctor, TAKE_OWNERSHIP, 2));

  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor, CENTRAL, ceres::DYNAMIC, 1, 2>(
          new EasyFunctor, TAKE_OWNERSHIP, 1));

  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor, CENTRAL, ceres::DYNAMIC, 1, 2>(
          new EasyFunctor, TAKE_OWNERSHIP, 2));

  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor, CENTRAL, ceres::DYNAMIC, 2, 1>(
          new EasyFunctor, TAKE_OWNERSHIP, 1));

  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor, CENTRAL, ceres::DYNAMIC, 2, 1>(
          new EasyFunctor, TAKE_OWNERSHIP, 2));
}

TEST(NumericDiffCostFunction,
     EasyCaseFunctorCentralDifferencesAndDynamicNumResiduals) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<EasyFunctor,
                                  CENTRAL,
                                  ceres::DYNAMIC,
                                  5,  /* size of x1 */
                                  5   /* size of x2 */>(
                                      new EasyFunctor, TAKE_OWNERSHIP, 3));
  EasyFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function, CENTRAL);
}

TEST(NumericDiffCostFunction, ExponentialFunctorRidders) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<ExponentialFunctor,
                                  RIDDERS,
                                  1,  /* number of residuals */
                                  1   /* size of x1 */>(
             new ExponentialFunctor));
  ExponentialFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function);
}

TEST(NumericDiffCostFunction, ExponentialCostFunctionRidders) {
  internal::scoped_ptr<CostFunction> cost_function;
  cost_function.reset(
      new NumericDiffCostFunction<ExponentialCostFunction,
                                  RIDDERS,
                                  1,  /* number of residuals */
                                  1   /* size of x1 */>(
             new ExponentialCostFunction));
  ExponentialFunctor functor;
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function);
}

TEST(NumericDiffCostFunction, RandomizedFunctorRidders) {
  internal::scoped_ptr<CostFunction> cost_function;
  NumericDiffOptions options;
  // Larger initial step size is chosen to produce robust results in the
  // presence of random noise.
  options.ridders_relative_initial_step_size = 10.0;

  cost_function.reset(
      new NumericDiffCostFunction<RandomizedFunctor,
                                  RIDDERS,
                                  1,  /* number of residuals */
                                  1   /* size of x1 */>(
             new RandomizedFunctor(kNoiseFactor, kRandomSeed), TAKE_OWNERSHIP,
             1, options));
  RandomizedFunctor functor (kNoiseFactor, kRandomSeed);
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function);
}

TEST(NumericDiffCostFunction, RandomizedCostFunctionRidders) {
  internal::scoped_ptr<CostFunction> cost_function;
  NumericDiffOptions options;
  // Larger initial step size is chosen to produce robust results in the
  // presence of random noise.
  options.ridders_relative_initial_step_size = 10.0;

  cost_function.reset(
      new NumericDiffCostFunction<RandomizedCostFunction,
                                  RIDDERS,
                                  1,  /* number of residuals */
                                  1   /* size of x1 */>(
             new RandomizedCostFunction(kNoiseFactor, kRandomSeed),
             TAKE_OWNERSHIP, 1, options));
  RandomizedFunctor functor (kNoiseFactor, kRandomSeed);
  functor.ExpectCostFunctionEvaluationIsNearlyCorrect(*cost_function);
}

struct OnlyFillsOneOutputFunctor {
  bool operator()(const double* x, double* output) const {
    output[0] = x[0];
    return true;
  }
};

TEST(NumericDiffCostFunction, PartiallyFilledResidualShouldFailEvaluation) {
  double parameter = 1.0;
  double jacobian[2];
  double residuals[2];
  double* parameters[] = {&parameter};
  double* jacobians[] = {jacobian};

  scoped_ptr<CostFunction> cost_function(
      new NumericDiffCostFunction<OnlyFillsOneOutputFunctor, CENTRAL, 2, 1>(
          new OnlyFillsOneOutputFunctor));
  InvalidateArray(2, jacobian);
  InvalidateArray(2, residuals);
  EXPECT_TRUE(cost_function->Evaluate(parameters, residuals, jacobians));
  EXPECT_FALSE(IsArrayValid(2, residuals));
  InvalidateArray(2, residuals);
  EXPECT_TRUE(cost_function->Evaluate(parameters, residuals, NULL));
  // We are only testing residuals here, because the Jacobians are
  // computed using finite differencing from the residuals, so unless
  // we introduce a validation step after every evaluation of
  // residuals inside NumericDiffCostFunction, there is no way of
  // ensuring that the Jacobian array is invalid.
  EXPECT_FALSE(IsArrayValid(2, residuals));
}

}  // namespace internal
}  // namespace ceres
