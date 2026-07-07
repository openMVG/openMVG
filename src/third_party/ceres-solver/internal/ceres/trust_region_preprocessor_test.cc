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

#include "ceres/trust_region_preprocessor.h"

#include <array>
#include <map>

#include "ceres/internal/config.h"
#include "ceres/ordered_groups.h"
#include "ceres/problem_impl.h"
#include "ceres/sized_cost_function.h"
#include "ceres/solver.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(TrustRegionPreprocessor, ZeroProblem) {
  ProblemImpl problem;
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(TrustRegionPreprocessor, ProblemWithInvalidParameterBlock) {
  ProblemImpl problem;
  double x = std::numeric_limits<double>::quiet_NaN();
  problem.AddParameterBlock(&x, 1);
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(TrustRegionPreprocessor, ParameterBlockBoundsAreInvalid) {
  ProblemImpl problem;
  double x = 1.0;
  problem.AddParameterBlock(&x, 1);
  problem.SetParameterUpperBound(&x, 0, 1.0);
  problem.SetParameterLowerBound(&x, 0, 2.0);
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(TrustRegionPreprocessor, ParameterBlockIsInfeasible) {
  ProblemImpl problem;
  double x = 3.0;
  problem.AddParameterBlock(&x, 1);
  problem.SetParameterUpperBound(&x, 0, 1.0);
  problem.SetParameterLowerBound(&x, 0, 2.0);
  problem.SetParameterBlockConstant(&x);
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

class FailingCostFunction : public SizedCostFunction<1, 1> {
 public:
  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    return false;
  }
};

TEST(TrustRegionPreprocessor, RemoveParameterBlocksFailed) {
  ProblemImpl problem;
  double x = 3.0;
  problem.AddResidualBlock(new FailingCostFunction, nullptr, &x);
  problem.SetParameterBlockConstant(&x);
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem, &pp));
}

TEST(TrustRegionPreprocessor, RemoveParameterBlocksSucceeds) {
  ProblemImpl problem;
  double x = 3.0;
  problem.AddParameterBlock(&x, 1);
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem, &pp));
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

    if (jacobians == nullptr) {
      return true;
    }

    std::array<int, sizeof...(Ns)> N{Ns...};
    for (size_t i = 0; i < N.size(); ++i) {
      if (jacobians[i] != nullptr) {
        MatrixRef j(jacobians[i], kNumResiduals, N[i]);
        j.setOnes();
        j *= kNumResiduals * N[i];
      }
    }

    return true;
  }
};

class LinearSolverAndEvaluatorCreationTest : public ::testing::Test {
 public:
  void SetUp() final {
    x_ = 1.0;
    y_ = 1.0;
    z_ = 1.0;
    problem_.AddResidualBlock(
        new DummyCostFunction<1, 1, 1>, nullptr, &x_, &y_);
    problem_.AddResidualBlock(
        new DummyCostFunction<1, 1, 1>, nullptr, &y_, &z_);
  }

  void PreprocessForGivenLinearSolverAndVerify(
      const LinearSolverType linear_solver_type) {
    Solver::Options options;
    options.linear_solver_type = linear_solver_type;
    TrustRegionPreprocessor preprocessor;
    PreprocessedProblem pp;
    EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
    EXPECT_EQ(pp.options.linear_solver_type, linear_solver_type);
    EXPECT_EQ(pp.linear_solver_options.type, linear_solver_type);
    EXPECT_EQ(pp.evaluator_options.linear_solver_type, linear_solver_type);
    EXPECT_TRUE(pp.linear_solver.get() != nullptr);
    EXPECT_TRUE(pp.evaluator.get() != nullptr);
  }

 protected:
  ProblemImpl problem_;
  double x_;
  double y_;
  double z_;
};

TEST_F(LinearSolverAndEvaluatorCreationTest, DenseQR) {
  PreprocessForGivenLinearSolverAndVerify(DENSE_QR);
}

TEST_F(LinearSolverAndEvaluatorCreationTest, DenseNormalCholesky) {
  PreprocessForGivenLinearSolverAndVerify(DENSE_NORMAL_CHOLESKY);
}

TEST_F(LinearSolverAndEvaluatorCreationTest, DenseSchur) {
  PreprocessForGivenLinearSolverAndVerify(DENSE_SCHUR);
}

#if !defined(CERES_NO_SPARSE)
TEST_F(LinearSolverAndEvaluatorCreationTest, SparseNormalCholesky) {
  PreprocessForGivenLinearSolverAndVerify(SPARSE_NORMAL_CHOLESKY);
}
#endif

#if !defined(CERES_NO_SPARSE)
TEST_F(LinearSolverAndEvaluatorCreationTest, SparseSchur) {
  PreprocessForGivenLinearSolverAndVerify(SPARSE_SCHUR);
}
#endif

TEST_F(LinearSolverAndEvaluatorCreationTest, CGNR) {
  PreprocessForGivenLinearSolverAndVerify(CGNR);
}

TEST_F(LinearSolverAndEvaluatorCreationTest, IterativeSchur) {
  PreprocessForGivenLinearSolverAndVerify(ITERATIVE_SCHUR);
}

TEST_F(LinearSolverAndEvaluatorCreationTest, MinimizerIsAwareOfBounds) {
  problem_.SetParameterLowerBound(&x_, 0, 0.0);
  Solver::Options options;
  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
  EXPECT_EQ(pp.options.linear_solver_type, options.linear_solver_type);
  EXPECT_EQ(pp.linear_solver_options.type, options.linear_solver_type);
  EXPECT_EQ(pp.evaluator_options.linear_solver_type,
            options.linear_solver_type);
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
  EXPECT_TRUE(pp.minimizer_options.is_constrained);
}

TEST_F(LinearSolverAndEvaluatorCreationTest, SchurTypeSolverWithBadOrdering) {
  Solver::Options options;
  options.linear_solver_type = DENSE_SCHUR;
  options.linear_solver_ordering = std::make_shared<ParameterBlockOrdering>();
  options.linear_solver_ordering->AddElementToGroup(&x_, 0);
  options.linear_solver_ordering->AddElementToGroup(&y_, 0);
  options.linear_solver_ordering->AddElementToGroup(&z_, 1);

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem_, &pp));
}

TEST_F(LinearSolverAndEvaluatorCreationTest, SchurTypeSolverWithGoodOrdering) {
  Solver::Options options;
  options.linear_solver_type = DENSE_SCHUR;
  options.linear_solver_ordering = std::make_shared<ParameterBlockOrdering>();
  options.linear_solver_ordering->AddElementToGroup(&x_, 0);
  options.linear_solver_ordering->AddElementToGroup(&z_, 0);
  options.linear_solver_ordering->AddElementToGroup(&y_, 1);

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
  EXPECT_EQ(pp.options.linear_solver_type, DENSE_SCHUR);
  EXPECT_EQ(pp.linear_solver_options.type, DENSE_SCHUR);
  EXPECT_EQ(pp.evaluator_options.linear_solver_type, DENSE_SCHUR);
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
}

TEST_F(LinearSolverAndEvaluatorCreationTest,
       SchurTypeSolverWithEmptyFirstEliminationGroup) {
  problem_.SetParameterBlockConstant(&x_);
  problem_.SetParameterBlockConstant(&z_);

  Solver::Options options;
  options.linear_solver_type = DENSE_SCHUR;
  options.linear_solver_ordering = std::make_shared<ParameterBlockOrdering>();
  options.linear_solver_ordering->AddElementToGroup(&x_, 0);
  options.linear_solver_ordering->AddElementToGroup(&z_, 0);
  options.linear_solver_ordering->AddElementToGroup(&y_, 1);

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
  EXPECT_EQ(pp.options.linear_solver_type, DENSE_QR);
  EXPECT_EQ(pp.linear_solver_options.type, DENSE_QR);
  EXPECT_EQ(pp.evaluator_options.linear_solver_type, DENSE_QR);
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
}

TEST_F(LinearSolverAndEvaluatorCreationTest,
       SchurTypeSolverWithEmptySecondEliminationGroup) {
  problem_.SetParameterBlockConstant(&y_);

  Solver::Options options;
  options.linear_solver_type = DENSE_SCHUR;
  options.linear_solver_ordering = std::make_shared<ParameterBlockOrdering>();
  options.linear_solver_ordering->AddElementToGroup(&x_, 0);
  options.linear_solver_ordering->AddElementToGroup(&z_, 0);
  options.linear_solver_ordering->AddElementToGroup(&y_, 1);

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
  EXPECT_EQ(pp.options.linear_solver_type, DENSE_SCHUR);
  EXPECT_EQ(pp.linear_solver_options.type, DENSE_SCHUR);
  EXPECT_EQ(pp.evaluator_options.linear_solver_type, DENSE_SCHUR);
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
}

TEST(TrustRegionPreprocessorTest, InnerIterationsWithOneParameterBlock) {
  ProblemImpl problem;
  double x = 1.0;
  problem.AddResidualBlock(new DummyCostFunction<1, 1>, nullptr, &x);

  Solver::Options options;
  options.use_inner_iterations = true;

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem, &pp));
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
  EXPECT_TRUE(pp.inner_iteration_minimizer.get() == nullptr);
}

TEST_F(LinearSolverAndEvaluatorCreationTest,
       InnerIterationsWithTwoParameterBlocks) {
  Solver::Options options;
  options.use_inner_iterations = true;

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
  EXPECT_TRUE(pp.inner_iteration_minimizer.get() != nullptr);
}

TEST_F(LinearSolverAndEvaluatorCreationTest, InvalidInnerIterationsOrdering) {
  Solver::Options options;
  options.use_inner_iterations = true;
  options.inner_iteration_ordering = std::make_shared<ParameterBlockOrdering>();
  options.inner_iteration_ordering->AddElementToGroup(&x_, 0);
  options.inner_iteration_ordering->AddElementToGroup(&z_, 0);
  options.inner_iteration_ordering->AddElementToGroup(&y_, 0);

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_FALSE(preprocessor.Preprocess(options, &problem_, &pp));
}

TEST_F(LinearSolverAndEvaluatorCreationTest, ValidInnerIterationsOrdering) {
  Solver::Options options;
  options.use_inner_iterations = true;
  options.inner_iteration_ordering = std::make_shared<ParameterBlockOrdering>();
  options.inner_iteration_ordering->AddElementToGroup(&x_, 0);
  options.inner_iteration_ordering->AddElementToGroup(&z_, 0);
  options.inner_iteration_ordering->AddElementToGroup(&y_, 1);

  TrustRegionPreprocessor preprocessor;
  PreprocessedProblem pp;
  EXPECT_TRUE(preprocessor.Preprocess(options, &problem_, &pp));
  EXPECT_TRUE(pp.linear_solver.get() != nullptr);
  EXPECT_TRUE(pp.evaluator.get() != nullptr);
  EXPECT_TRUE(pp.inner_iteration_minimizer.get() != nullptr);
}

}  // namespace internal
}  // namespace ceres
