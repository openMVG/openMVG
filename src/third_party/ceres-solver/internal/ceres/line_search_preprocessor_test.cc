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

#include <map>

#include "ceres/problem_impl.h"
#include "ceres/sized_cost_function.h"
#include "ceres/solver.h"
#include "ceres/line_search_preprocessor.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(LineSearchPreprocessor, ZeroProblem) {
  ProblemImpl problem;
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  LineSearchPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(LineSearchPreprocessor, ProblemWithInvalidParameterBlock) {
  ProblemImpl problem;
  double x = std::numeric_limits<double>::quiet_NaN();
  problem.AddParameterBlock(&x, 1);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  LineSearchPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(LineSearchPreprocessor, ParameterBlockHasBounds) {
  ProblemImpl problem;
  double x = 1.0;
  problem.AddParameterBlock(&x, 1);
  problem.SetParameterUpperBound(&x, 0, 1.0);
  problem.SetParameterLowerBound(&x, 0, 2.0);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  LineSearchPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

class FailingCostFunction : public SizedCostFunction<1, 1> {
 public:
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const {
    return false;
  }
};

TEST(LineSearchPreprocessor, RemoveParameterBlocksFailed) {
  ProblemImpl problem;
  double x = 3.0;
  problem.AddResidualBlock(new FailingCostFunction, NULL, &x);
  problem.SetParameterBlockConstant(&x);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  LineSearchPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(LineSearchPreprocessor, RemoveParameterBlocksSucceeds) {
  ProblemImpl problem;
  double x = 3.0;
  problem.AddParameterBlock(&x, 1);
  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;
  LineSearchPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem, &pp));
}

template<int kNumResiduals, int N1 = 0, int N2 = 0, int N3 = 0>
class DummyCostFunction : public SizedCostFunction<kNumResiduals, N1, N2, N3> {
 public:
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const {
    return true;
  }
};

TEST(LineSearchPreprocessor, NormalOperation) {
  ProblemImpl problem;
  double x = 1.0;
  double y = 1.0;
  double z = 1.0;
  problem.AddResidualBlock(new DummyCostFunction<1, 1, 1>, NULL, &x, &y);
  problem.AddResidualBlock(new DummyCostFunction<1, 1, 1>, NULL, &y, &z);

  Solver::Options options;
  options.minimizer_type = LINE_SEARCH;

  LineSearchPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem, &pp));
  EXPECT_EQ(pp.evaluator_options.linear_solver_type, CGNR);
  EXPECT_TRUE(pp.evaluator.get() != NULL);
}

}  // namespace internal
}  // namespace ceres
