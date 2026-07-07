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
//         keir@google.com (Keir Mierle)

#include "ceres/problem.h"

#include <memory>
#include <string>
#include <vector>

#include "ceres/autodiff_cost_function.h"
#include "ceres/casts.h"
#include "ceres/cost_function.h"
#include "ceres/crs_matrix.h"
#include "ceres/evaluator_test_utils.h"
#include "ceres/internal/eigen.h"
#include "ceres/loss_function.h"
#include "ceres/map_util.h"
#include "ceres/parameter_block.h"
#include "ceres/problem_impl.h"
#include "ceres/program.h"
#include "ceres/sized_cost_function.h"
#include "ceres/sparse_matrix.h"
#include "ceres/types.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// The following three classes are for the purposes of defining
// function signatures. They have dummy Evaluate functions.

// Trivial cost function that accepts a single argument.
class UnaryCostFunction : public CostFunction {
 public:
  UnaryCostFunction(int num_residuals, int32_t parameter_block_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block_size);
  }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 1;
    }
    return true;
  }
};

// Trivial cost function that accepts two arguments.
class BinaryCostFunction : public CostFunction {
 public:
  BinaryCostFunction(int num_residuals,
                     int32_t parameter_block1_size,
                     int32_t parameter_block2_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block1_size);
    mutable_parameter_block_sizes()->push_back(parameter_block2_size);
  }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 2;
    }
    return true;
  }
};

// Trivial cost function that accepts three arguments.
class TernaryCostFunction : public CostFunction {
 public:
  TernaryCostFunction(int num_residuals,
                      int32_t parameter_block1_size,
                      int32_t parameter_block2_size,
                      int32_t parameter_block3_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block1_size);
    mutable_parameter_block_sizes()->push_back(parameter_block2_size);
    mutable_parameter_block_sizes()->push_back(parameter_block3_size);
  }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 3;
    }
    return true;
  }
};

TEST(Problem, MoveConstructor) {
  Problem src;
  double x;
  src.AddParameterBlock(&x, 1);
  Problem dst(std::move(src));
  EXPECT_TRUE(dst.HasParameterBlock(&x));
}

TEST(Problem, MoveAssignment) {
  Problem src;
  double x;
  src.AddParameterBlock(&x, 1);
  Problem dst;
  dst = std::move(src);
  EXPECT_TRUE(dst.HasParameterBlock(&x));
}

TEST(Problem, AddResidualWithNullCostFunctionDies) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);
  problem.AddParameterBlock(z, 5);

  EXPECT_DEATH_IF_SUPPORTED(problem.AddResidualBlock(nullptr, nullptr, x),
                            "cost_function != nullptr");
}

TEST(Problem, AddResidualWithIncorrectNumberOfParameterBlocksDies) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);
  problem.AddParameterBlock(z, 5);

  // UnaryCostFunction takes only one parameter, but two are passed.
  EXPECT_DEATH_IF_SUPPORTED(
      problem.AddResidualBlock(new UnaryCostFunction(2, 3), nullptr, x, y),
      "num_parameter_blocks");
}

TEST(Problem, AddResidualWithDifferentSizesOnTheSameVariableDies) {
  double x[3];

  Problem problem;
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), nullptr, x);
  EXPECT_DEATH_IF_SUPPORTED(
      problem.AddResidualBlock(
          new UnaryCostFunction(2, 4 /* 4 != 3 */), nullptr, x),
      "different block sizes");
}

TEST(Problem, AddResidualWithDuplicateParametersDies) {
  double x[3], z[5];

  Problem problem;
  EXPECT_DEATH_IF_SUPPORTED(
      problem.AddResidualBlock(new BinaryCostFunction(2, 3, 3), nullptr, x, x),
      "Duplicate parameter blocks");
  EXPECT_DEATH_IF_SUPPORTED(
      problem.AddResidualBlock(
          new TernaryCostFunction(1, 5, 3, 5), nullptr, z, x, z),
      "Duplicate parameter blocks");
}

TEST(Problem, AddResidualWithIncorrectSizesOfParameterBlockDies) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);
  problem.AddParameterBlock(z, 5);

  // The cost function expects the size of the second parameter, z, to be 4
  // instead of 5 as declared above. This is fatal.
  EXPECT_DEATH_IF_SUPPORTED(
      problem.AddResidualBlock(new BinaryCostFunction(2, 3, 4), nullptr, x, z),
      "different block sizes");
}

TEST(Problem, AddResidualAddsDuplicatedParametersOnlyOnce) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), nullptr, x);
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), nullptr, x);
  problem.AddResidualBlock(new UnaryCostFunction(2, 4), nullptr, y);
  problem.AddResidualBlock(new UnaryCostFunction(2, 5), nullptr, z);

  EXPECT_EQ(3, problem.NumParameterBlocks());
  EXPECT_EQ(12, problem.NumParameters());
}

TEST(Problem, AddParameterWithDifferentSizesOnTheSameVariableDies) {
  double x[3], y[4];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);

  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(x, 4),
                            "different block sizes");
}

static double* IntToPtr(int i) {
  return reinterpret_cast<double*>(sizeof(double) * i);  // NOLINT
}

TEST(Problem, AddParameterWithAliasedParametersDies) {
  // Layout is
  //
  //   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
  //                 [x] x  x  x  x          [y] y  y
  //         o==o==o                 o==o==o           o==o
  //               o--o--o     o--o--o     o--o  o--o--o
  //
  // Parameter block additions are tested as listed above; expected successful
  // ones marked with o==o and aliasing ones marked with o--o.

  Problem problem;
  problem.AddParameterBlock(IntToPtr(5), 5);   // x
  problem.AddParameterBlock(IntToPtr(13), 3);  // y

  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(4), 2),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(4), 3),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(4), 9),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(8), 3),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(12), 2),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(14), 3),
                            "Aliasing detected");

  // These ones should work.
  problem.AddParameterBlock(IntToPtr(2), 3);
  problem.AddParameterBlock(IntToPtr(10), 3);
  problem.AddParameterBlock(IntToPtr(16), 2);

  ASSERT_EQ(5, problem.NumParameterBlocks());
}

TEST(Problem, AddParameterIgnoresDuplicateCalls) {
  double x[3], y[4];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);

  // Creating parameter blocks multiple times is ignored.
  problem.AddParameterBlock(x, 3);
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), nullptr, x);

  // ... even repeatedly.
  problem.AddParameterBlock(x, 3);
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), nullptr, x);

  // More parameters are fine.
  problem.AddParameterBlock(y, 4);
  problem.AddResidualBlock(new UnaryCostFunction(2, 4), nullptr, y);

  EXPECT_EQ(2, problem.NumParameterBlocks());
  EXPECT_EQ(7, problem.NumParameters());
}

class DestructorCountingCostFunction : public SizedCostFunction<3, 4, 5> {
 public:
  explicit DestructorCountingCostFunction(int* num_destructions)
      : num_destructions_(num_destructions) {}

  ~DestructorCountingCostFunction() override { *num_destructions_ += 1; }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    return true;
  }

 private:
  int* num_destructions_;
};

TEST(Problem, ReusedCostFunctionsAreOnlyDeletedOnce) {
  double y[4], z[5];
  int num_destructions = 0;

  // Add a cost function multiple times and check to make sure that
  // the destructor on the cost function is only called once.
  {
    Problem problem;
    problem.AddParameterBlock(y, 4);
    problem.AddParameterBlock(z, 5);

    CostFunction* cost = new DestructorCountingCostFunction(&num_destructions);
    problem.AddResidualBlock(cost, nullptr, y, z);
    problem.AddResidualBlock(cost, nullptr, y, z);
    problem.AddResidualBlock(cost, nullptr, y, z);
    EXPECT_EQ(3, problem.NumResidualBlocks());
  }

  // Check that the destructor was called only once.
  CHECK_EQ(num_destructions, 1);
}

TEST(Problem, GetCostFunctionForResidualBlock) {
  double x[3];
  Problem problem;
  CostFunction* cost_function = new UnaryCostFunction(2, 3);
  const ResidualBlockId residual_block =
      problem.AddResidualBlock(cost_function, nullptr, x);
  EXPECT_EQ(problem.GetCostFunctionForResidualBlock(residual_block),
            cost_function);
  EXPECT_TRUE(problem.GetLossFunctionForResidualBlock(residual_block) ==
              nullptr);
}

TEST(Problem, GetLossFunctionForResidualBlock) {
  double x[3];
  Problem problem;
  CostFunction* cost_function = new UnaryCostFunction(2, 3);
  LossFunction* loss_function = new TrivialLoss();
  const ResidualBlockId residual_block =
      problem.AddResidualBlock(cost_function, loss_function, x);
  EXPECT_EQ(problem.GetCostFunctionForResidualBlock(residual_block),
            cost_function);
  EXPECT_EQ(problem.GetLossFunctionForResidualBlock(residual_block),
            loss_function);
}

TEST(Problem, CostFunctionsAreDeletedEvenWithRemovals) {
  double y[4], z[5], w[4];
  int num_destructions = 0;
  {
    Problem problem;
    problem.AddParameterBlock(y, 4);
    problem.AddParameterBlock(z, 5);

    CostFunction* cost_yz =
        new DestructorCountingCostFunction(&num_destructions);
    CostFunction* cost_wz =
        new DestructorCountingCostFunction(&num_destructions);
    ResidualBlock* r_yz = problem.AddResidualBlock(cost_yz, nullptr, y, z);
    ResidualBlock* r_wz = problem.AddResidualBlock(cost_wz, nullptr, w, z);
    EXPECT_EQ(2, problem.NumResidualBlocks());

    problem.RemoveResidualBlock(r_yz);
    CHECK_EQ(num_destructions, 1);
    problem.RemoveResidualBlock(r_wz);
    CHECK_EQ(num_destructions, 2);

    EXPECT_EQ(0, problem.NumResidualBlocks());
  }
  CHECK_EQ(num_destructions, 2);
}

// Make the dynamic problem tests (e.g. for removing residual blocks)
// parameterized on whether the low-latency mode is enabled or not.
//
// This tests against ProblemImpl instead of Problem in order to inspect the
// state of the resulting Program; this is difficult with only the thin Problem
// interface.
struct DynamicProblem : public ::testing::TestWithParam<bool> {
  DynamicProblem() {
    Problem::Options options;
    options.enable_fast_removal = GetParam();
    problem = std::make_unique<ProblemImpl>(options);
  }

  ParameterBlock* GetParameterBlock(int block) {
    return problem->program().parameter_blocks()[block];
  }
  ResidualBlock* GetResidualBlock(int block) {
    return problem->program().residual_blocks()[block];
  }

  bool HasResidualBlock(ResidualBlock* residual_block) {
    bool have_residual_block = true;
    if (GetParam()) {
      have_residual_block &=
          (problem->residual_block_set().find(residual_block) !=
           problem->residual_block_set().end());
    }
    have_residual_block &=
        find(problem->program().residual_blocks().begin(),
             problem->program().residual_blocks().end(),
             residual_block) != problem->program().residual_blocks().end();
    return have_residual_block;
  }

  int NumResidualBlocks() {
    // Verify that the hash set of residuals is maintained consistently.
    if (GetParam()) {
      EXPECT_EQ(problem->residual_block_set().size(),
                problem->NumResidualBlocks());
    }
    return problem->NumResidualBlocks();
  }

  // The next block of functions until the end are only for testing the
  // residual block removals.
  void ExpectParameterBlockContainsResidualBlock(
      double* values, ResidualBlock* residual_block) {
    ParameterBlock* parameter_block =
        FindOrDie(problem->parameter_map(), values);
    EXPECT_TRUE(ContainsKey(*(parameter_block->mutable_residual_blocks()),
                            residual_block));
  }

  void ExpectSize(double* values, int size) {
    ParameterBlock* parameter_block =
        FindOrDie(problem->parameter_map(), values);
    EXPECT_EQ(size, parameter_block->mutable_residual_blocks()->size());
  }

  // Degenerate case.
  void ExpectParameterBlockContains(double* values) { ExpectSize(values, 0); }

  void ExpectParameterBlockContains(double* values, ResidualBlock* r1) {
    ExpectSize(values, 1);
    ExpectParameterBlockContainsResidualBlock(values, r1);
  }

  void ExpectParameterBlockContains(double* values,
                                    ResidualBlock* r1,
                                    ResidualBlock* r2) {
    ExpectSize(values, 2);
    ExpectParameterBlockContainsResidualBlock(values, r1);
    ExpectParameterBlockContainsResidualBlock(values, r2);
  }

  void ExpectParameterBlockContains(double* values,
                                    ResidualBlock* r1,
                                    ResidualBlock* r2,
                                    ResidualBlock* r3) {
    ExpectSize(values, 3);
    ExpectParameterBlockContainsResidualBlock(values, r1);
    ExpectParameterBlockContainsResidualBlock(values, r2);
    ExpectParameterBlockContainsResidualBlock(values, r3);
  }

  void ExpectParameterBlockContains(double* values,
                                    ResidualBlock* r1,
                                    ResidualBlock* r2,
                                    ResidualBlock* r3,
                                    ResidualBlock* r4) {
    ExpectSize(values, 4);
    ExpectParameterBlockContainsResidualBlock(values, r1);
    ExpectParameterBlockContainsResidualBlock(values, r2);
    ExpectParameterBlockContainsResidualBlock(values, r3);
    ExpectParameterBlockContainsResidualBlock(values, r4);
  }

  std::unique_ptr<ProblemImpl> problem;
  double y[4], z[5], w[3];
};

TEST(Problem, SetParameterBlockConstantWithUnknownPtrDies) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);

  EXPECT_DEATH_IF_SUPPORTED(problem.SetParameterBlockConstant(y),
                            "Parameter block not found:");
}

TEST(Problem, SetParameterBlockVariableWithUnknownPtrDies) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);

  EXPECT_DEATH_IF_SUPPORTED(problem.SetParameterBlockVariable(y),
                            "Parameter block not found:");
}

TEST(Problem, IsParameterBlockConstant) {
  double x1[3];
  double x2[3];

  Problem problem;
  problem.AddParameterBlock(x1, 3);
  problem.AddParameterBlock(x2, 3);

  EXPECT_FALSE(problem.IsParameterBlockConstant(x1));
  EXPECT_FALSE(problem.IsParameterBlockConstant(x2));

  problem.SetParameterBlockConstant(x1);
  EXPECT_TRUE(problem.IsParameterBlockConstant(x1));
  EXPECT_FALSE(problem.IsParameterBlockConstant(x2));

  problem.SetParameterBlockConstant(x2);
  EXPECT_TRUE(problem.IsParameterBlockConstant(x1));
  EXPECT_TRUE(problem.IsParameterBlockConstant(x2));

  problem.SetParameterBlockVariable(x1);
  EXPECT_FALSE(problem.IsParameterBlockConstant(x1));
  EXPECT_TRUE(problem.IsParameterBlockConstant(x2));
}

TEST(Problem, IsParameterBlockConstantWithUnknownPtrDies) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);

  EXPECT_DEATH_IF_SUPPORTED(problem.IsParameterBlockConstant(y),
                            "Parameter block not found:");
}

TEST(Problem, SetManifoldWithUnknownPtrDies) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);

  EXPECT_DEATH_IF_SUPPORTED(problem.SetManifold(y, new EuclideanManifold<3>),
                            "Parameter block not found:");
}

TEST(Problem, RemoveParameterBlockWithUnknownPtrDies) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);

  EXPECT_DEATH_IF_SUPPORTED(problem.RemoveParameterBlock(y),
                            "Parameter block not found:");
}

TEST(Problem, GetManifold) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 2);

  Manifold* manifold = new EuclideanManifold<3>;
  problem.SetManifold(x, manifold);
  EXPECT_EQ(problem.GetManifold(x), manifold);
  EXPECT_TRUE(problem.GetManifold(y) == nullptr);
}

TEST(Problem, HasManifold) {
  double x[3];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 2);

  Manifold* manifold = new EuclideanManifold<3>;
  problem.SetManifold(x, manifold);
  EXPECT_TRUE(problem.HasManifold(x));
  EXPECT_FALSE(problem.HasManifold(y));
}

TEST(Problem, RepeatedAddParameterBlockResetsManifold) {
  double x[4];
  double y[2];

  Problem problem;
  problem.AddParameterBlock(x, 4, new SubsetManifold(4, {0, 1}));
  problem.AddParameterBlock(y, 2);

  EXPECT_FALSE(problem.HasManifold(y));

  EXPECT_TRUE(problem.HasManifold(x));
  EXPECT_EQ(problem.ParameterBlockSize(x), 4);
  EXPECT_EQ(problem.ParameterBlockTangentSize(x), 2);
  EXPECT_EQ(problem.GetManifold(x)->AmbientSize(), 4);
  EXPECT_EQ(problem.GetManifold(x)->TangentSize(), 2);

  problem.AddParameterBlock(x, 4, static_cast<Manifold*>(nullptr));
  EXPECT_FALSE(problem.HasManifold(x));
  EXPECT_EQ(problem.ParameterBlockSize(x), 4);
  EXPECT_EQ(problem.ParameterBlockTangentSize(x), 4);
  EXPECT_EQ(problem.GetManifold(x), nullptr);

  problem.AddParameterBlock(x, 4, new SubsetManifold(4, {0, 1, 2}));
  problem.AddParameterBlock(y, 2);
  EXPECT_TRUE(problem.HasManifold(x));
  EXPECT_EQ(problem.ParameterBlockSize(x), 4);
  EXPECT_EQ(problem.ParameterBlockTangentSize(x), 1);
  EXPECT_EQ(problem.GetManifold(x)->AmbientSize(), 4);
  EXPECT_EQ(problem.GetManifold(x)->TangentSize(), 1);
}

TEST(Problem, ParameterBlockQueryTestUsingManifold) {
  double x[3];
  double y[4];
  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);

  std::vector<int> constant_parameters;
  constant_parameters.push_back(0);
  problem.SetManifold(x, new SubsetManifold(3, constant_parameters));
  EXPECT_EQ(problem.ParameterBlockSize(x), 3);
  EXPECT_EQ(problem.ParameterBlockTangentSize(x), 2);
  EXPECT_EQ(problem.ParameterBlockTangentSize(y), 4);

  std::vector<double*> parameter_blocks;
  problem.GetParameterBlocks(&parameter_blocks);
  EXPECT_EQ(parameter_blocks.size(), 2);
  EXPECT_NE(parameter_blocks[0], parameter_blocks[1]);
  EXPECT_TRUE(parameter_blocks[0] == x || parameter_blocks[0] == y);
  EXPECT_TRUE(parameter_blocks[1] == x || parameter_blocks[1] == y);

  EXPECT_TRUE(problem.HasParameterBlock(x));
  problem.RemoveParameterBlock(x);
  EXPECT_FALSE(problem.HasParameterBlock(x));
  problem.GetParameterBlocks(&parameter_blocks);
  EXPECT_EQ(parameter_blocks.size(), 1);
  EXPECT_TRUE(parameter_blocks[0] == y);
}

TEST(Problem, ParameterBlockQueryTest) {
  double x[3];
  double y[4];
  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);

  std::vector<int> constant_parameters;
  constant_parameters.push_back(0);
  problem.SetManifold(x, new SubsetManifold(3, constant_parameters));
  EXPECT_EQ(problem.ParameterBlockSize(x), 3);
  EXPECT_EQ(problem.ParameterBlockTangentSize(x), 2);
  EXPECT_EQ(problem.ParameterBlockTangentSize(y), 4);

  std::vector<double*> parameter_blocks;
  problem.GetParameterBlocks(&parameter_blocks);
  EXPECT_EQ(parameter_blocks.size(), 2);
  EXPECT_NE(parameter_blocks[0], parameter_blocks[1]);
  EXPECT_TRUE(parameter_blocks[0] == x || parameter_blocks[0] == y);
  EXPECT_TRUE(parameter_blocks[1] == x || parameter_blocks[1] == y);

  EXPECT_TRUE(problem.HasParameterBlock(x));
  problem.RemoveParameterBlock(x);
  EXPECT_FALSE(problem.HasParameterBlock(x));
  problem.GetParameterBlocks(&parameter_blocks);
  EXPECT_EQ(parameter_blocks.size(), 1);
  EXPECT_TRUE(parameter_blocks[0] == y);
}

TEST_P(DynamicProblem, RemoveParameterBlockWithNoResiduals) {
  problem->AddParameterBlock(y, 4);
  problem->AddParameterBlock(z, 5);
  problem->AddParameterBlock(w, 3);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(y, GetParameterBlock(0)->user_state());
  EXPECT_EQ(z, GetParameterBlock(1)->user_state());
  EXPECT_EQ(w, GetParameterBlock(2)->user_state());

  // w is at the end, which might break the swapping logic so try adding and
  // removing it.
  problem->RemoveParameterBlock(w);
  ASSERT_EQ(2, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(y, GetParameterBlock(0)->user_state());
  EXPECT_EQ(z, GetParameterBlock(1)->user_state());
  problem->AddParameterBlock(w, 3);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(y, GetParameterBlock(0)->user_state());
  EXPECT_EQ(z, GetParameterBlock(1)->user_state());
  EXPECT_EQ(w, GetParameterBlock(2)->user_state());

  // Now remove z, which is in the middle, and add it back.
  problem->RemoveParameterBlock(z);
  ASSERT_EQ(2, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(y, GetParameterBlock(0)->user_state());
  EXPECT_EQ(w, GetParameterBlock(1)->user_state());
  problem->AddParameterBlock(z, 5);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(y, GetParameterBlock(0)->user_state());
  EXPECT_EQ(w, GetParameterBlock(1)->user_state());
  EXPECT_EQ(z, GetParameterBlock(2)->user_state());

  // Now remove everything.
  // y
  problem->RemoveParameterBlock(y);
  ASSERT_EQ(2, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(z, GetParameterBlock(0)->user_state());
  EXPECT_EQ(w, GetParameterBlock(1)->user_state());

  // z
  problem->RemoveParameterBlock(z);
  ASSERT_EQ(1, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(w, GetParameterBlock(0)->user_state());

  // w
  problem->RemoveParameterBlock(w);
  EXPECT_EQ(0, problem->NumParameterBlocks());
  EXPECT_EQ(0, NumResidualBlocks());
}

TEST_P(DynamicProblem, RemoveParameterBlockWithResiduals) {
  problem->AddParameterBlock(y, 4);
  problem->AddParameterBlock(z, 5);
  problem->AddParameterBlock(w, 3);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  EXPECT_EQ(y, GetParameterBlock(0)->user_state());
  EXPECT_EQ(z, GetParameterBlock(1)->user_state());
  EXPECT_EQ(w, GetParameterBlock(2)->user_state());

  // clang-format off

  // Add all combinations of cost functions.
  CostFunction* cost_yzw = new TernaryCostFunction(1, 4, 5, 3);
  CostFunction* cost_yz  = new BinaryCostFunction (1, 4, 5);
  CostFunction* cost_yw  = new BinaryCostFunction (1, 4, 3);
  CostFunction* cost_zw  = new BinaryCostFunction (1, 5, 3);
  CostFunction* cost_y   = new UnaryCostFunction  (1, 4);
  CostFunction* cost_z   = new UnaryCostFunction  (1, 5);
  CostFunction* cost_w   = new UnaryCostFunction  (1, 3);

  ResidualBlock* r_yzw = problem->AddResidualBlock(cost_yzw, nullptr, y, z, w);
  ResidualBlock* r_yz  = problem->AddResidualBlock(cost_yz,  nullptr, y, z);
  ResidualBlock* r_yw  = problem->AddResidualBlock(cost_yw,  nullptr, y, w);
  ResidualBlock* r_zw  = problem->AddResidualBlock(cost_zw,  nullptr, z, w);
  ResidualBlock* r_y   = problem->AddResidualBlock(cost_y,   nullptr, y);
  ResidualBlock* r_z   = problem->AddResidualBlock(cost_z,   nullptr, z);
  ResidualBlock* r_w   = problem->AddResidualBlock(cost_w,   nullptr, w);

  EXPECT_EQ(3, problem->NumParameterBlocks());
  EXPECT_EQ(7, NumResidualBlocks());

  // Remove w, which should remove r_yzw, r_yw, r_zw, r_w.
  problem->RemoveParameterBlock(w);
  ASSERT_EQ(2, problem->NumParameterBlocks());
  ASSERT_EQ(3, NumResidualBlocks());

  ASSERT_FALSE(HasResidualBlock(r_yzw));
  ASSERT_TRUE (HasResidualBlock(r_yz ));
  ASSERT_FALSE(HasResidualBlock(r_yw ));
  ASSERT_FALSE(HasResidualBlock(r_zw ));
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_TRUE (HasResidualBlock(r_z  ));
  ASSERT_FALSE(HasResidualBlock(r_w  ));

  // Remove z, which will remove almost everything else.
  problem->RemoveParameterBlock(z);
  ASSERT_EQ(1, problem->NumParameterBlocks());
  ASSERT_EQ(1, NumResidualBlocks());

  ASSERT_FALSE(HasResidualBlock(r_yzw));
  ASSERT_FALSE(HasResidualBlock(r_yz ));
  ASSERT_FALSE(HasResidualBlock(r_yw ));
  ASSERT_FALSE(HasResidualBlock(r_zw ));
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_FALSE(HasResidualBlock(r_z  ));
  ASSERT_FALSE(HasResidualBlock(r_w  ));

  // Remove y; all gone.
  problem->RemoveParameterBlock(y);
  EXPECT_EQ(0, problem->NumParameterBlocks());
  EXPECT_EQ(0, NumResidualBlocks());

  // clang-format on
}

TEST_P(DynamicProblem, RemoveResidualBlock) {
  problem->AddParameterBlock(y, 4);
  problem->AddParameterBlock(z, 5);
  problem->AddParameterBlock(w, 3);

  // clang-format off

  // Add all combinations of cost functions.
  CostFunction* cost_yzw = new TernaryCostFunction(1, 4, 5, 3);
  CostFunction* cost_yz  = new BinaryCostFunction (1, 4, 5);
  CostFunction* cost_yw  = new BinaryCostFunction (1, 4, 3);
  CostFunction* cost_zw  = new BinaryCostFunction (1, 5, 3);
  CostFunction* cost_y   = new UnaryCostFunction  (1, 4);
  CostFunction* cost_z   = new UnaryCostFunction  (1, 5);
  CostFunction* cost_w   = new UnaryCostFunction  (1, 3);

  ResidualBlock* r_yzw = problem->AddResidualBlock(cost_yzw, nullptr, y, z, w);
  ResidualBlock* r_yz  = problem->AddResidualBlock(cost_yz,  nullptr, y, z);
  ResidualBlock* r_yw  = problem->AddResidualBlock(cost_yw,  nullptr, y, w);
  ResidualBlock* r_zw  = problem->AddResidualBlock(cost_zw,  nullptr, z, w);
  ResidualBlock* r_y   = problem->AddResidualBlock(cost_y,   nullptr, y);
  ResidualBlock* r_z   = problem->AddResidualBlock(cost_z,   nullptr, z);
  ResidualBlock* r_w   = problem->AddResidualBlock(cost_w,   nullptr, w);

  if (GetParam()) {
    // In this test parameterization, there should be back-pointers from the
    // parameter blocks to the residual blocks.
    ExpectParameterBlockContains(y, r_yzw, r_yz, r_yw, r_y);
    ExpectParameterBlockContains(z, r_yzw, r_yz, r_zw, r_z);
    ExpectParameterBlockContains(w, r_yzw, r_yw, r_zw, r_w);
  } else {
    // Otherwise, nothing.
    EXPECT_TRUE(GetParameterBlock(0)->mutable_residual_blocks() == nullptr);
    EXPECT_TRUE(GetParameterBlock(1)->mutable_residual_blocks() == nullptr);
    EXPECT_TRUE(GetParameterBlock(2)->mutable_residual_blocks() == nullptr);
  }
  EXPECT_EQ(3, problem->NumParameterBlocks());
  EXPECT_EQ(7, NumResidualBlocks());

  // Remove each residual and check the state after each removal.

  // Remove r_yzw.
  problem->RemoveResidualBlock(r_yzw);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(6, NumResidualBlocks());
  if (GetParam()) {
    ExpectParameterBlockContains(y, r_yz, r_yw, r_y);
    ExpectParameterBlockContains(z, r_yz, r_zw, r_z);
    ExpectParameterBlockContains(w, r_yw, r_zw, r_w);
  }
  ASSERT_TRUE (HasResidualBlock(r_yz ));
  ASSERT_TRUE (HasResidualBlock(r_yw ));
  ASSERT_TRUE (HasResidualBlock(r_zw ));
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_TRUE (HasResidualBlock(r_z  ));
  ASSERT_TRUE (HasResidualBlock(r_w  ));

  // Remove r_yw.
  problem->RemoveResidualBlock(r_yw);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(5, NumResidualBlocks());
  if (GetParam()) {
    ExpectParameterBlockContains(y, r_yz, r_y);
    ExpectParameterBlockContains(z, r_yz, r_zw, r_z);
    ExpectParameterBlockContains(w, r_zw, r_w);
  }
  ASSERT_TRUE (HasResidualBlock(r_yz ));
  ASSERT_TRUE (HasResidualBlock(r_zw ));
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_TRUE (HasResidualBlock(r_z  ));
  ASSERT_TRUE (HasResidualBlock(r_w  ));

  // Remove r_zw.
  problem->RemoveResidualBlock(r_zw);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(4, NumResidualBlocks());
  if (GetParam()) {
    ExpectParameterBlockContains(y, r_yz, r_y);
    ExpectParameterBlockContains(z, r_yz, r_z);
    ExpectParameterBlockContains(w, r_w);
  }
  ASSERT_TRUE (HasResidualBlock(r_yz ));
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_TRUE (HasResidualBlock(r_z  ));
  ASSERT_TRUE (HasResidualBlock(r_w  ));

  // Remove r_w.
  problem->RemoveResidualBlock(r_w);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(3, NumResidualBlocks());
  if (GetParam()) {
    ExpectParameterBlockContains(y, r_yz, r_y);
    ExpectParameterBlockContains(z, r_yz, r_z);
    ExpectParameterBlockContains(w);
  }
  ASSERT_TRUE (HasResidualBlock(r_yz ));
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_TRUE (HasResidualBlock(r_z  ));

  // Remove r_yz.
  problem->RemoveResidualBlock(r_yz);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(2, NumResidualBlocks());
  if (GetParam()) {
    ExpectParameterBlockContains(y, r_y);
    ExpectParameterBlockContains(z, r_z);
    ExpectParameterBlockContains(w);
  }
  ASSERT_TRUE (HasResidualBlock(r_y  ));
  ASSERT_TRUE (HasResidualBlock(r_z  ));

  // Remove the last two.
  problem->RemoveResidualBlock(r_z);
  problem->RemoveResidualBlock(r_y);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(0, NumResidualBlocks());
  if (GetParam()) {
    ExpectParameterBlockContains(y);
    ExpectParameterBlockContains(z);
    ExpectParameterBlockContains(w);
  }

  // clang-format on
}

TEST_P(DynamicProblem, RemoveInvalidResidualBlockDies) {
  problem->AddParameterBlock(y, 4);
  problem->AddParameterBlock(z, 5);
  problem->AddParameterBlock(w, 3);

  // clang-format off

  // Add all combinations of cost functions.
  CostFunction* cost_yzw = new TernaryCostFunction(1, 4, 5, 3);
  CostFunction* cost_yz  = new BinaryCostFunction (1, 4, 5);
  CostFunction* cost_yw  = new BinaryCostFunction (1, 4, 3);
  CostFunction* cost_zw  = new BinaryCostFunction (1, 5, 3);
  CostFunction* cost_y   = new UnaryCostFunction  (1, 4);
  CostFunction* cost_z   = new UnaryCostFunction  (1, 5);
  CostFunction* cost_w   = new UnaryCostFunction  (1, 3);

  ResidualBlock* r_yzw = problem->AddResidualBlock(cost_yzw, nullptr, y, z, w);
  ResidualBlock* r_yz  = problem->AddResidualBlock(cost_yz,  nullptr, y, z);
  ResidualBlock* r_yw  = problem->AddResidualBlock(cost_yw,  nullptr, y, w);
  ResidualBlock* r_zw  = problem->AddResidualBlock(cost_zw,  nullptr, z, w);
  ResidualBlock* r_y   = problem->AddResidualBlock(cost_y,   nullptr, y);
  ResidualBlock* r_z   = problem->AddResidualBlock(cost_z,   nullptr, z);
  ResidualBlock* r_w   = problem->AddResidualBlock(cost_w,   nullptr, w);

  // clang-format on

  // Remove r_yzw.
  problem->RemoveResidualBlock(r_yzw);
  ASSERT_EQ(3, problem->NumParameterBlocks());
  ASSERT_EQ(6, NumResidualBlocks());
  // Attempt to remove r_yzw again.
  EXPECT_DEATH_IF_SUPPORTED(problem->RemoveResidualBlock(r_yzw), "not found");

  // Attempt to remove a cast pointer never added as a residual.
  int trash_memory = 1234;
  auto* invalid_residual = reinterpret_cast<ResidualBlock*>(&trash_memory);
  EXPECT_DEATH_IF_SUPPORTED(problem->RemoveResidualBlock(invalid_residual),
                            "not found");

  // Remove a parameter block, which in turn removes the dependent residuals
  // then attempt to remove them directly.
  problem->RemoveParameterBlock(z);
  ASSERT_EQ(2, problem->NumParameterBlocks());
  ASSERT_EQ(3, NumResidualBlocks());
  EXPECT_DEATH_IF_SUPPORTED(problem->RemoveResidualBlock(r_yz), "not found");
  EXPECT_DEATH_IF_SUPPORTED(problem->RemoveResidualBlock(r_zw), "not found");
  EXPECT_DEATH_IF_SUPPORTED(problem->RemoveResidualBlock(r_z), "not found");

  problem->RemoveResidualBlock(r_yw);
  problem->RemoveResidualBlock(r_w);
  problem->RemoveResidualBlock(r_y);
}

// Check that a null-terminated array, a, has the same elements as b.
template <typename T>
void ExpectVectorContainsUnordered(const T* a, const std::vector<T>& b) {
  // Compute the size of a.
  int size = 0;
  while (a[size]) {
    ++size;
  }
  ASSERT_EQ(size, b.size());

  // Sort a.
  std::vector<T> a_sorted(size);
  copy(a, a + size, a_sorted.begin());
  sort(a_sorted.begin(), a_sorted.end());

  // Sort b.
  std::vector<T> b_sorted(b);
  sort(b_sorted.begin(), b_sorted.end());

  // Compare.
  for (int i = 0; i < size; ++i) {
    EXPECT_EQ(a_sorted[i], b_sorted[i]);
  }
}

static void ExpectProblemHasResidualBlocks(
    const ProblemImpl& problem,
    const ResidualBlockId* expected_residual_blocks) {
  std::vector<ResidualBlockId> residual_blocks;
  problem.GetResidualBlocks(&residual_blocks);
  ExpectVectorContainsUnordered(expected_residual_blocks, residual_blocks);
}

TEST_P(DynamicProblem, GetXXXBlocksForYYYBlock) {
  problem->AddParameterBlock(y, 4);
  problem->AddParameterBlock(z, 5);
  problem->AddParameterBlock(w, 3);

  // clang-format off

  // Add all combinations of cost functions.
  CostFunction* cost_yzw = new TernaryCostFunction(1, 4, 5, 3);
  CostFunction* cost_yz  = new BinaryCostFunction (1, 4, 5);
  CostFunction* cost_yw  = new BinaryCostFunction (1, 4, 3);
  CostFunction* cost_zw  = new BinaryCostFunction (1, 5, 3);
  CostFunction* cost_y   = new UnaryCostFunction  (1, 4);
  CostFunction* cost_z   = new UnaryCostFunction  (1, 5);
  CostFunction* cost_w   = new UnaryCostFunction  (1, 3);

  ResidualBlock* r_yzw = problem->AddResidualBlock(cost_yzw, nullptr, y, z, w);
  {
    ResidualBlockId expected_residuals[] = {r_yzw, nullptr};
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }
  ResidualBlock* r_yz  = problem->AddResidualBlock(cost_yz,  nullptr, y, z);
  {
    ResidualBlockId expected_residuals[] = {r_yzw, r_yz, nullptr};
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }
  ResidualBlock* r_yw  = problem->AddResidualBlock(cost_yw,  nullptr, y, w);
  {
    ResidualBlock *expected_residuals[] = {r_yzw, r_yz, r_yw, nullptr};
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }
  ResidualBlock* r_zw  = problem->AddResidualBlock(cost_zw,  nullptr, z, w);
  {
    ResidualBlock *expected_residuals[] = {r_yzw, r_yz, r_yw, r_zw, nullptr};
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }
  ResidualBlock* r_y   = problem->AddResidualBlock(cost_y,   nullptr, y);
  {
    ResidualBlock *expected_residuals[] = {r_yzw, r_yz, r_yw, r_zw, r_y, nullptr};
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }
  ResidualBlock* r_z   = problem->AddResidualBlock(cost_z,   nullptr, z);
  {
    ResidualBlock *expected_residuals[] = {
      r_yzw, r_yz, r_yw, r_zw, r_y, r_z, nullptr
    };
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }
  ResidualBlock* r_w   = problem->AddResidualBlock(cost_w,   nullptr, w);
  {
    ResidualBlock *expected_residuals[] = {
      r_yzw, r_yz, r_yw, r_zw, r_y, r_z, r_w, nullptr
    };
    ExpectProblemHasResidualBlocks(*problem, expected_residuals);
  }

  std::vector<double*> parameter_blocks;
  std::vector<ResidualBlockId> residual_blocks;

  // Check GetResidualBlocksForParameterBlock() for all parameter blocks.
  struct GetResidualBlocksForParameterBlockTestCase {
    double* parameter_block;
    ResidualBlockId expected_residual_blocks[10];
  };
  GetResidualBlocksForParameterBlockTestCase get_residual_blocks_cases[] = {
    { y, { r_yzw, r_yz, r_yw, r_y, nullptr} },
    { z, { r_yzw, r_yz, r_zw, r_z, nullptr} },
    { w, { r_yzw, r_yw, r_zw, r_w, nullptr} },
    { nullptr, { nullptr } }
  };
  for (int i = 0; get_residual_blocks_cases[i].parameter_block; ++i) {
    problem->GetResidualBlocksForParameterBlock(
        get_residual_blocks_cases[i].parameter_block,
        &residual_blocks);
    ExpectVectorContainsUnordered(
        get_residual_blocks_cases[i].expected_residual_blocks,
        residual_blocks);
  }

  // Check GetParameterBlocksForResidualBlock() for all residual blocks.
  struct GetParameterBlocksForResidualBlockTestCase {
    ResidualBlockId residual_block;
    double* expected_parameter_blocks[10];
  };
  GetParameterBlocksForResidualBlockTestCase get_parameter_blocks_cases[] = {
    { r_yzw, { y, z, w, nullptr } },
    { r_yz , { y, z, nullptr } },
    { r_yw , { y, w, nullptr } },
    { r_zw , { z, w, nullptr } },
    { r_y  , { y, nullptr } },
    { r_z  , { z, nullptr } },
    { r_w  , { w, nullptr } },
    { nullptr, { nullptr } }
  };
  for (int i = 0; get_parameter_blocks_cases[i].residual_block; ++i) {
    problem->GetParameterBlocksForResidualBlock(
        get_parameter_blocks_cases[i].residual_block,
        &parameter_blocks);
    ExpectVectorContainsUnordered(
        get_parameter_blocks_cases[i].expected_parameter_blocks,
        parameter_blocks);
  }

  // clang-format on
}

INSTANTIATE_TEST_SUITE_P(OptionsInstantiation,
                         DynamicProblem,
                         ::testing::Values(true, false));

// Test for Problem::Evaluate

// r_i = i - (j + 1) * x_ij^2
template <int kNumResiduals, int kNumParameterBlocks>
class QuadraticCostFunction : public CostFunction {
 public:
  QuadraticCostFunction() {
    CHECK_GT(kNumResiduals, 0);
    CHECK_GT(kNumParameterBlocks, 0);
    set_num_residuals(kNumResiduals);
    for (int i = 0; i < kNumParameterBlocks; ++i) {
      mutable_parameter_block_sizes()->push_back(kNumResiduals);
    }
  }

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const final {
    for (int i = 0; i < kNumResiduals; ++i) {
      residuals[i] = i;
      for (int j = 0; j < kNumParameterBlocks; ++j) {
        residuals[i] -= (j + 1.0) * parameters[j][i] * parameters[j][i];
      }
    }

    if (jacobians == nullptr) {
      return true;
    }

    for (int j = 0; j < kNumParameterBlocks; ++j) {
      if (jacobians[j] != nullptr) {
        MatrixRef(jacobians[j], kNumResiduals, kNumResiduals) =
            (-2.0 * (j + 1.0) * ConstVectorRef(parameters[j], kNumResiduals))
                .asDiagonal();
      }
    }

    return true;
  }
};

// Convert a CRSMatrix to a dense Eigen matrix.
static void CRSToDenseMatrix(const CRSMatrix& input, Matrix* output) {
  CHECK(output != nullptr);
  Matrix& m = *output;
  m.resize(input.num_rows, input.num_cols);
  m.setZero();
  for (int row = 0; row < input.num_rows; ++row) {
    for (int j = input.rows[row]; j < input.rows[row + 1]; ++j) {
      const int col = input.cols[j];
      m(row, col) = input.values[j];
    }
  }
}

class ProblemEvaluateTest : public ::testing::Test {
 protected:
  void SetUp() override {
    for (int i = 0; i < 6; ++i) {
      parameters_[i] = static_cast<double>(i + 1);
    }

    parameter_blocks_.push_back(parameters_);
    parameter_blocks_.push_back(parameters_ + 2);
    parameter_blocks_.push_back(parameters_ + 4);

    CostFunction* cost_function = new QuadraticCostFunction<2, 2>;

    // f(x, y)
    residual_blocks_.push_back(problem_.AddResidualBlock(
        cost_function, nullptr, parameters_, parameters_ + 2));
    // g(y, z)
    residual_blocks_.push_back(problem_.AddResidualBlock(
        cost_function, nullptr, parameters_ + 2, parameters_ + 4));
    // h(z, x)
    residual_blocks_.push_back(problem_.AddResidualBlock(
        cost_function, nullptr, parameters_ + 4, parameters_));
  }

  void TearDown() override { EXPECT_TRUE(problem_.program().IsValid()); }

  void EvaluateAndCompare(const Problem::EvaluateOptions& options,
                          const int expected_num_rows,
                          const int expected_num_cols,
                          const double expected_cost,
                          const double* expected_residuals,
                          const double* expected_gradient,
                          const double* expected_jacobian) {
    double cost;
    std::vector<double> residuals;
    std::vector<double> gradient;
    CRSMatrix jacobian;

    EXPECT_TRUE(
        problem_.Evaluate(options,
                          &cost,
                          expected_residuals != nullptr ? &residuals : nullptr,
                          expected_gradient != nullptr ? &gradient : nullptr,
                          expected_jacobian != nullptr ? &jacobian : nullptr));

    if (expected_residuals != nullptr) {
      EXPECT_EQ(residuals.size(), expected_num_rows);
    }

    if (expected_gradient != nullptr) {
      EXPECT_EQ(gradient.size(), expected_num_cols);
    }

    if (expected_jacobian != nullptr) {
      EXPECT_EQ(jacobian.num_rows, expected_num_rows);
      EXPECT_EQ(jacobian.num_cols, expected_num_cols);
    }

    Matrix dense_jacobian;
    if (expected_jacobian != nullptr) {
      CRSToDenseMatrix(jacobian, &dense_jacobian);
    }

    CompareEvaluations(expected_num_rows,
                       expected_num_cols,
                       expected_cost,
                       expected_residuals,
                       expected_gradient,
                       expected_jacobian,
                       cost,
                       !residuals.empty() ? &residuals[0] : nullptr,
                       !gradient.empty() ? &gradient[0] : nullptr,
                       dense_jacobian.data());
  }

  void CheckAllEvaluationCombinations(const Problem::EvaluateOptions& options,
                                      const ExpectedEvaluation& expected) {
    for (int i = 0; i < 8; ++i) {
      EvaluateAndCompare(options,
                         expected.num_rows,
                         expected.num_cols,
                         expected.cost,
                         (i & 1) ? expected.residuals : nullptr,
                         (i & 2) ? expected.gradient : nullptr,
                         (i & 4) ? expected.jacobian : nullptr);
    }
  }

  ProblemImpl problem_;
  double parameters_[6];
  std::vector<double*> parameter_blocks_;
  std::vector<ResidualBlockId> residual_blocks_;
};

TEST_F(ProblemEvaluateTest, MultipleParameterAndResidualBlocks) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 6,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -59.0, -87.0,  // g
      -27.0, -43.0   // h
    },
    // Gradient
    {  146.0,  484.0,   // x
       582.0, 1256.0,   // y
      1450.0, 2604.0,   // z
    },
    // Jacobian
    //                       x             y             z
    { /* f(x, y) */ -2.0,  0.0, -12.0,   0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0, -16.0,   0.0,   0.0,
      /* g(y, z) */  0.0,  0.0,  -6.0,   0.0, -20.0,   0.0,
                     0.0,  0.0,   0.0,  -8.0,   0.0, -24.0,
      /* h(z, x) */ -4.0,  0.0,   0.0,   0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0,   0.0,   0.0, -12.0
    }
  };
  // clang-format on

  CheckAllEvaluationCombinations(Problem::EvaluateOptions(), expected);
}

TEST_F(ProblemEvaluateTest, ParameterAndResidualBlocksPassedInOptions) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 6,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -59.0, -87.0,  // g
      -27.0, -43.0   // h
    },
    // Gradient
    {  146.0,  484.0,   // x
       582.0, 1256.0,   // y
      1450.0, 2604.0,   // z
    },
    // Jacobian
    //                       x             y             z
    { /* f(x, y) */ -2.0,  0.0, -12.0,   0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0, -16.0,   0.0,   0.0,
      /* g(y, z) */  0.0,  0.0,  -6.0,   0.0, -20.0,   0.0,
                     0.0,  0.0,   0.0,  -8.0,   0.0, -24.0,
      /* h(z, x) */ -4.0,  0.0,   0.0,   0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0,   0.0,   0.0, -12.0
    }
  };
  // clang-format on

  Problem::EvaluateOptions evaluate_options;
  evaluate_options.parameter_blocks = parameter_blocks_;
  evaluate_options.residual_blocks = residual_blocks_;
  CheckAllEvaluationCombinations(evaluate_options, expected);
}

TEST_F(ProblemEvaluateTest, ReorderedResidualBlocks) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 6,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -27.0, -43.0,  // h
      -59.0, -87.0   // g
    },
    // Gradient
    {  146.0,  484.0,   // x
       582.0, 1256.0,   // y
      1450.0, 2604.0,   // z
    },
    // Jacobian
    //                       x             y             z
    { /* f(x, y) */ -2.0,  0.0, -12.0,   0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0, -16.0,   0.0,   0.0,
      /* h(z, x) */ -4.0,  0.0,   0.0,   0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0,   0.0,   0.0, -12.0,
      /* g(y, z) */  0.0,  0.0,  -6.0,   0.0, -20.0,   0.0,
                     0.0,  0.0,   0.0,  -8.0,   0.0, -24.0
    }
  };
  // clang-format on

  Problem::EvaluateOptions evaluate_options;
  evaluate_options.parameter_blocks = parameter_blocks_;

  // f, h, g
  evaluate_options.residual_blocks.push_back(residual_blocks_[0]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[2]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[1]);

  CheckAllEvaluationCombinations(evaluate_options, expected);
}

TEST_F(ProblemEvaluateTest,
       ReorderedResidualBlocksAndReorderedParameterBlocks) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 6,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -27.0, -43.0,  // h
      -59.0, -87.0   // g
    },
    // Gradient
    {  1450.0, 2604.0,   // z
        582.0, 1256.0,   // y
        146.0,  484.0,   // x
    },
     // Jacobian
    //                       z             y             x
    { /* f(x, y) */   0.0,   0.0, -12.0,   0.0,  -2.0,   0.0,
                      0.0,   0.0,   0.0, -16.0,   0.0,  -4.0,
      /* h(z, x) */ -10.0,   0.0,   0.0,   0.0,  -4.0,   0.0,
                      0.0, -12.0,   0.0,   0.0,   0.0,  -8.0,
      /* g(y, z) */ -20.0,   0.0,  -6.0,   0.0,   0.0,   0.0,
                      0.0, -24.0,   0.0,  -8.0,   0.0,   0.0
    }
  };
  // clang-format on

  Problem::EvaluateOptions evaluate_options;
  // z, y, x
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[2]);
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[1]);
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[0]);

  // f, h, g
  evaluate_options.residual_blocks.push_back(residual_blocks_[0]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[2]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[1]);

  CheckAllEvaluationCombinations(evaluate_options, expected);
}

TEST_F(ProblemEvaluateTest, ConstantParameterBlock) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 6,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -59.0, -87.0,  // g
      -27.0, -43.0   // h
    },

    // Gradient
    {  146.0,  484.0,  // x
         0.0,    0.0,  // y
      1450.0, 2604.0,  // z
    },

    // Jacobian
    //                       x             y             z
    { /* f(x, y) */ -2.0,  0.0,   0.0,   0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0,   0.0,   0.0,   0.0,
      /* g(y, z) */  0.0,  0.0,   0.0,   0.0, -20.0,   0.0,
                     0.0,  0.0,   0.0,   0.0,   0.0, -24.0,
      /* h(z, x) */ -4.0,  0.0,   0.0,   0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0,   0.0,   0.0, -12.0
    }
  };
  // clang-format on

  problem_.SetParameterBlockConstant(parameters_ + 2);
  CheckAllEvaluationCombinations(Problem::EvaluateOptions(), expected);
}

TEST_F(ProblemEvaluateTest, ExcludedAResidualBlock) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    4, 6,
    // Cost
    2082.0,
    // Residuals
    { -19.0, -35.0,  // f
      -27.0, -43.0   // h
    },
    // Gradient
    {  146.0,  484.0,   // x
       228.0,  560.0,   // y
       270.0,  516.0,   // z
    },
    // Jacobian
    //                       x             y             z
    { /* f(x, y) */ -2.0,  0.0, -12.0,   0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0, -16.0,   0.0,   0.0,
      /* h(z, x) */ -4.0,  0.0,   0.0,   0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0,   0.0,   0.0, -12.0
    }
  };
  // clang-format on

  Problem::EvaluateOptions evaluate_options;
  evaluate_options.residual_blocks.push_back(residual_blocks_[0]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[2]);

  CheckAllEvaluationCombinations(evaluate_options, expected);
}

TEST_F(ProblemEvaluateTest, ExcludedParameterBlock) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 4,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -59.0, -87.0,  // g
      -27.0, -43.0   // h
    },

    // Gradient
    {  146.0,  484.0,  // x
      1450.0, 2604.0,  // z
    },

    // Jacobian
    //                       x             z
    { /* f(x, y) */ -2.0,  0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0,   0.0,
      /* g(y, z) */  0.0,  0.0, -20.0,   0.0,
                     0.0,  0.0,   0.0, -24.0,
      /* h(z, x) */ -4.0,  0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0, -12.0
    }
  };
  // clang-format on

  Problem::EvaluateOptions evaluate_options;
  // x, z
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[0]);
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[2]);
  evaluate_options.residual_blocks = residual_blocks_;
  CheckAllEvaluationCombinations(evaluate_options, expected);
}

TEST_F(ProblemEvaluateTest, ExcludedParameterBlockAndExcludedResidualBlock) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    4, 4,
    // Cost
    6318.0,
    // Residuals
    { -19.0, -35.0,  // f
      -59.0, -87.0,  // g
    },

    // Gradient
    {   38.0,  140.0,  // x
      1180.0, 2088.0,  // z
    },

    // Jacobian
    //                       x             z
    { /* f(x, y) */ -2.0,  0.0,   0.0,   0.0,
                     0.0, -4.0,   0.0,   0.0,
      /* g(y, z) */  0.0,  0.0, -20.0,   0.0,
                     0.0,  0.0,   0.0, -24.0,
    }
  };
  // clang-format on

  Problem::EvaluateOptions evaluate_options;
  // x, z
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[0]);
  evaluate_options.parameter_blocks.push_back(parameter_blocks_[2]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[0]);
  evaluate_options.residual_blocks.push_back(residual_blocks_[1]);

  CheckAllEvaluationCombinations(evaluate_options, expected);
}

TEST_F(ProblemEvaluateTest, Manifold) {
  // clang-format off
  ExpectedEvaluation expected = {
    // Rows/columns
    6, 5,
    // Cost
    7607.0,
    // Residuals
    { -19.0, -35.0,  // f
      -59.0, -87.0,  // g
      -27.0, -43.0   // h
    },
    // Gradient
    {  146.0,  484.0,  // x
      1256.0,          // y with SubsetManifold
      1450.0, 2604.0,  // z
    },
    // Jacobian
    //                       x      y             z
    { /* f(x, y) */ -2.0,  0.0,   0.0,   0.0,   0.0,
                     0.0, -4.0, -16.0,   0.0,   0.0,
      /* g(y, z) */  0.0,  0.0,   0.0, -20.0,   0.0,
                     0.0,  0.0,  -8.0,   0.0, -24.0,
      /* h(z, x) */ -4.0,  0.0,   0.0, -10.0,   0.0,
                     0.0, -8.0,   0.0,   0.0, -12.0
    }
  };
  // clang-format on

  std::vector<int> constant_parameters;
  constant_parameters.push_back(0);
  problem_.SetManifold(parameters_ + 2,
                       new SubsetManifold(2, constant_parameters));

  CheckAllEvaluationCombinations(Problem::EvaluateOptions(), expected);
}

struct IdentityFunctor {
  template <typename T>
  bool operator()(const T* x, const T* y, T* residuals) const {
    residuals[0] = x[0];
    residuals[1] = x[1];
    residuals[2] = y[0];
    residuals[3] = y[1];
    residuals[4] = y[2];
    return true;
  }

  static CostFunction* Create() {
    return new AutoDiffCostFunction<IdentityFunctor, 5, 2, 3>(
        new IdentityFunctor);
  }
};

class ProblemEvaluateResidualBlockTest : public ::testing::Test {
 public:
  static constexpr bool kApplyLossFunction = true;
  static constexpr bool kDoNotApplyLossFunction = false;
  static constexpr bool kNewPoint = true;
  static constexpr bool kNotNewPoint = false;
  static double loss_function_scale_;

 protected:
  ProblemImpl problem_;
  double x_[2] = {1, 2};
  double y_[3] = {1, 2, 3};
};

double ProblemEvaluateResidualBlockTest::loss_function_scale_ = 2.0;

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockNoLossFunctionFullEval) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Matrix expected_dfdx = Matrix::Zero(5, 2);
  expected_dfdx.block(0, 0, 2, 2) = Matrix::Identity(2, 2);
  Matrix expected_dfdy = Matrix::Zero(5, 3);
  expected_dfdy.block(2, 0, 3, 3) = Matrix::Identity(3, 3);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdx - actual_dfdx).norm() / actual_dfdx.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdx;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockNoLossFunctionNullEval) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             nullptr,
                                             nullptr,
                                             nullptr));
}

TEST_F(ProblemEvaluateResidualBlockTest, OneResidualBlockNoLossFunctionCost) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             nullptr,
                                             nullptr));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockNoLossFunctionCostAndResidual) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             nullptr));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockNoLossFunctionCostResidualAndOneJacobian) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Matrix expected_dfdx = Matrix::Zero(5, 2);
  expected_dfdx.block(0, 0, 2, 2) = Matrix::Identity(2, 2);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  double* jacobians[2] = {actual_dfdx.data(), nullptr};
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdx - actual_dfdx).norm() / actual_dfdx.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdx;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockNoLossFunctionResidual) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Vector actual_f(5);
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             nullptr,
                                             actual_f.data(),
                                             nullptr));

  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
}

TEST_F(ProblemEvaluateResidualBlockTest, OneResidualBlockWithLossFunction) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(),
                                new ScaledLoss(nullptr, 2.0, TAKE_OWNERSHIP),
                                x_,
                                y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  expected_f *= std::sqrt(loss_function_scale_);
  Matrix expected_dfdx = Matrix::Zero(5, 2);
  expected_dfdx.block(0, 0, 2, 2) = Matrix::Identity(2, 2);
  expected_dfdx *= std::sqrt(loss_function_scale_);
  Matrix expected_dfdy = Matrix::Zero(5, 3);
  expected_dfdy.block(2, 0, 3, 3) = Matrix::Identity(3, 3);
  expected_dfdy *= std::sqrt(loss_function_scale_);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdx - actual_dfdx).norm() / actual_dfdx.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdx;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockWithLossFunctionDisabled) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(),
                                new ScaledLoss(nullptr, 2.0, TAKE_OWNERSHIP),
                                x_,
                                y_);
  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Matrix expected_dfdx = Matrix::Zero(5, 2);
  expected_dfdx.block(0, 0, 2, 2) = Matrix::Identity(2, 2);
  Matrix expected_dfdy = Matrix::Zero(5, 3);
  expected_dfdy.block(2, 0, 3, 3) = Matrix::Identity(3, 3);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kDoNotApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdx - actual_dfdx).norm() / actual_dfdx.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdx;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST_F(ProblemEvaluateResidualBlockTest, OneResidualBlockWithOneManifold) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  problem_.SetManifold(x_, new SubsetManifold(2, {1}));

  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Matrix expected_dfdx = Matrix::Zero(5, 1);
  expected_dfdx.block(0, 0, 1, 1) = Matrix::Identity(1, 1);
  Matrix expected_dfdy = Matrix::Zero(5, 3);
  expected_dfdy.block(2, 0, 3, 3) = Matrix::Identity(3, 3);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 1);
  Matrix actual_dfdy(5, 3);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdx - actual_dfdx).norm() / actual_dfdx.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdx;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST_F(ProblemEvaluateResidualBlockTest, OneResidualBlockWithTwoManifolds) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  problem_.SetManifold(x_, new SubsetManifold(2, {1}));
  problem_.SetManifold(y_, new SubsetManifold(3, {2}));

  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Matrix expected_dfdx = Matrix::Zero(5, 1);
  expected_dfdx.block(0, 0, 1, 1) = Matrix::Identity(1, 1);
  Matrix expected_dfdy = Matrix::Zero(5, 2);
  expected_dfdy.block(2, 0, 2, 2) = Matrix::Identity(2, 2);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 1);
  Matrix actual_dfdy(5, 2);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdx - actual_dfdx).norm() / actual_dfdx.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdx;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockWithOneConstantParameterBlock) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  problem_.SetParameterBlockConstant(x_);

  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  Matrix expected_dfdy = Matrix::Zero(5, 3);
  expected_dfdy.block(2, 0, 3, 3) = Matrix::Identity(3, 3);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);

  // Try evaluating both Jacobians, this should fail.
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_FALSE(problem_.EvaluateResidualBlock(residual_block_id,
                                              kApplyLossFunction,
                                              kNewPoint,
                                              &actual_cost,
                                              actual_f.data(),
                                              jacobians));

  jacobians[0] = nullptr;
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockWithAllConstantParameterBlocks) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  problem_.SetParameterBlockConstant(x_);
  problem_.SetParameterBlockConstant(y_);

  Vector expected_f(5);
  expected_f << 1, 2, 1, 2, 3;
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);

  // Try evaluating with one or more Jacobians, this should fail.
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_FALSE(problem_.EvaluateResidualBlock(residual_block_id,
                                              kApplyLossFunction,
                                              kNewPoint,
                                              &actual_cost,
                                              actual_f.data(),
                                              jacobians));

  jacobians[0] = nullptr;
  EXPECT_FALSE(problem_.EvaluateResidualBlock(residual_block_id,
                                              kApplyLossFunction,
                                              kNewPoint,
                                              &actual_cost,
                                              actual_f.data(),
                                              jacobians));
  jacobians[1] = nullptr;
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));

  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
}

TEST_F(ProblemEvaluateResidualBlockTest,
       OneResidualBlockWithOneParameterBlockConstantAndParameterBlockChanged) {
  ResidualBlockId residual_block_id =
      problem_.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);
  problem_.SetParameterBlockConstant(x_);

  x_[0] = 2;
  y_[2] = 1;
  Vector expected_f(5);
  expected_f << 2, 2, 1, 2, 1;
  Matrix expected_dfdy = Matrix::Zero(5, 3);
  expected_dfdy.block(2, 0, 3, 3) = Matrix::Identity(3, 3);
  double expected_cost = expected_f.squaredNorm() / 2.0;

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);

  // Try evaluating with one or more Jacobians, this should fail.
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_FALSE(problem_.EvaluateResidualBlock(residual_block_id,
                                              kApplyLossFunction,
                                              kNewPoint,
                                              &actual_cost,
                                              actual_f.data(),
                                              jacobians));

  jacobians[0] = nullptr;
  EXPECT_TRUE(problem_.EvaluateResidualBlock(residual_block_id,
                                             kApplyLossFunction,
                                             kNewPoint,
                                             &actual_cost,
                                             actual_f.data(),
                                             jacobians));
  EXPECT_NEAR(std::abs(expected_cost - actual_cost) / actual_cost,
              0,
              std::numeric_limits<double>::epsilon())
      << actual_cost;
  EXPECT_NEAR((expected_f - actual_f).norm() / actual_f.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_f;
  EXPECT_NEAR((expected_dfdy - actual_dfdy).norm() / actual_dfdy.norm(),
              0,
              std::numeric_limits<double>::epsilon())
      << actual_dfdy;
}

TEST(Problem, SetAndGetParameterLowerBound) {
  Problem problem;
  double x[] = {1.0, 2.0};
  problem.AddParameterBlock(x, 2);

  EXPECT_EQ(problem.GetParameterLowerBound(x, 0),
            -std::numeric_limits<double>::max());
  EXPECT_EQ(problem.GetParameterLowerBound(x, 1),
            -std::numeric_limits<double>::max());

  problem.SetParameterLowerBound(x, 0, -1.0);
  EXPECT_EQ(problem.GetParameterLowerBound(x, 0), -1.0);
  EXPECT_EQ(problem.GetParameterLowerBound(x, 1),
            -std::numeric_limits<double>::max());

  problem.SetParameterLowerBound(x, 0, -2.0);
  EXPECT_EQ(problem.GetParameterLowerBound(x, 0), -2.0);
  EXPECT_EQ(problem.GetParameterLowerBound(x, 1),
            -std::numeric_limits<double>::max());

  problem.SetParameterLowerBound(x, 0, -std::numeric_limits<double>::max());
  EXPECT_EQ(problem.GetParameterLowerBound(x, 0),
            -std::numeric_limits<double>::max());
  EXPECT_EQ(problem.GetParameterLowerBound(x, 1),
            -std::numeric_limits<double>::max());
}

TEST(Problem, SetAndGetParameterUpperBound) {
  Problem problem;
  double x[] = {1.0, 2.0};
  problem.AddParameterBlock(x, 2);

  EXPECT_EQ(problem.GetParameterUpperBound(x, 0),
            std::numeric_limits<double>::max());
  EXPECT_EQ(problem.GetParameterUpperBound(x, 1),
            std::numeric_limits<double>::max());

  problem.SetParameterUpperBound(x, 0, -1.0);
  EXPECT_EQ(problem.GetParameterUpperBound(x, 0), -1.0);
  EXPECT_EQ(problem.GetParameterUpperBound(x, 1),
            std::numeric_limits<double>::max());

  problem.SetParameterUpperBound(x, 0, -2.0);
  EXPECT_EQ(problem.GetParameterUpperBound(x, 0), -2.0);
  EXPECT_EQ(problem.GetParameterUpperBound(x, 1),
            std::numeric_limits<double>::max());

  problem.SetParameterUpperBound(x, 0, std::numeric_limits<double>::max());
  EXPECT_EQ(problem.GetParameterUpperBound(x, 0),
            std::numeric_limits<double>::max());
  EXPECT_EQ(problem.GetParameterUpperBound(x, 1),
            std::numeric_limits<double>::max());
}

TEST(Problem, SetManifoldTwice) {
  Problem problem;
  double x[] = {1.0, 2.0, 3.0};
  problem.AddParameterBlock(x, 3);
  problem.SetManifold(x, new SubsetManifold(3, {1}));
  EXPECT_EQ(problem.GetManifold(x)->AmbientSize(), 3);
  EXPECT_EQ(problem.GetManifold(x)->TangentSize(), 2);

  problem.SetManifold(x, new SubsetManifold(3, {0, 1}));
  EXPECT_EQ(problem.GetManifold(x)->AmbientSize(), 3);
  EXPECT_EQ(problem.GetManifold(x)->TangentSize(), 1);
}

TEST(Problem, SetManifoldAndThenClearItWithNull) {
  Problem problem;
  double x[] = {1.0, 2.0, 3.0};
  problem.AddParameterBlock(x, 3);
  problem.SetManifold(x, new SubsetManifold(3, {1}));
  EXPECT_EQ(problem.GetManifold(x)->AmbientSize(), 3);
  EXPECT_EQ(problem.GetManifold(x)->TangentSize(), 2);

  problem.SetManifold(x, nullptr);
  EXPECT_EQ(problem.GetManifold(x), nullptr);
  EXPECT_EQ(problem.ParameterBlockTangentSize(x), 3);
  EXPECT_EQ(problem.ParameterBlockSize(x), 3);
}

TEST(Solver, ZeroTangentSizedManifoldMeansParameterBlockIsConstant) {
  double x = 0.0;
  double y = 1.0;
  Problem problem;
  problem.AddResidualBlock(new BinaryCostFunction(1, 1, 1), nullptr, &x, &y);
  problem.SetManifold(&y, new SubsetManifold(1, {0}));
  EXPECT_TRUE(problem.IsParameterBlockConstant(&y));
}

class MockEvaluationCallback : public EvaluationCallback {
 public:
  MOCK_METHOD2(PrepareForEvaluation, void(bool, bool));
};

TEST(ProblemEvaluate, CallsEvaluationCallbackWithoutJacobian) {
  constexpr bool kDoNotComputeJacobians = false;
  constexpr bool kNewPoint = true;

  MockEvaluationCallback evaluation_callback;
  EXPECT_CALL(evaluation_callback,
              PrepareForEvaluation(kDoNotComputeJacobians, kNewPoint))
      .Times(1);

  Problem::Options options;
  options.evaluation_callback = &evaluation_callback;
  ProblemImpl problem(options);
  double x_[2] = {1, 2};
  double y_[3] = {1, 2, 3};
  problem.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);

  double actual_cost;
  EXPECT_TRUE(problem.Evaluate(
      Problem::EvaluateOptions(), &actual_cost, nullptr, nullptr, nullptr));
}

TEST(ProblemEvaluate, CallsEvaluationCallbackWithJacobian) {
  constexpr bool kComputeJacobians = true;
  constexpr bool kNewPoint = true;

  MockEvaluationCallback evaluation_callback;
  EXPECT_CALL(evaluation_callback,
              PrepareForEvaluation(kComputeJacobians, kNewPoint))
      .Times(1);

  Problem::Options options;
  options.evaluation_callback = &evaluation_callback;
  ProblemImpl problem(options);
  double x_[2] = {1, 2};
  double y_[3] = {1, 2, 3};
  problem.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);

  double actual_cost;
  ceres::CRSMatrix jacobian;
  EXPECT_TRUE(problem.Evaluate(
      Problem::EvaluateOptions(), &actual_cost, nullptr, nullptr, &jacobian));
}

TEST(ProblemEvaluateResidualBlock, NewPointCallsEvaluationCallback) {
  constexpr bool kComputeJacobians = true;
  constexpr bool kNewPoint = true;

  MockEvaluationCallback evaluation_callback;
  EXPECT_CALL(evaluation_callback,
              PrepareForEvaluation(kComputeJacobians, kNewPoint))
      .Times(1);

  Problem::Options options;
  options.evaluation_callback = &evaluation_callback;
  ProblemImpl problem(options);
  double x_[2] = {1, 2};
  double y_[3] = {1, 2, 3};
  ResidualBlockId residual_block_id =
      problem.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem.EvaluateResidualBlock(
      residual_block_id, true, true, &actual_cost, actual_f.data(), jacobians));
}

TEST(ProblemEvaluateResidualBlock, OldPointCallsEvaluationCallback) {
  constexpr bool kComputeJacobians = true;
  constexpr bool kOldPoint = false;

  MockEvaluationCallback evaluation_callback;
  EXPECT_CALL(evaluation_callback,
              PrepareForEvaluation(kComputeJacobians, kOldPoint))
      .Times(1);

  Problem::Options options;
  options.evaluation_callback = &evaluation_callback;
  ProblemImpl problem(options);
  double x_[2] = {1, 2};
  double y_[3] = {1, 2, 3};
  ResidualBlockId residual_block_id =
      problem.AddResidualBlock(IdentityFunctor::Create(), nullptr, x_, y_);

  double actual_cost;
  Vector actual_f(5);
  Matrix actual_dfdx(5, 2);
  Matrix actual_dfdy(5, 3);
  double* jacobians[2] = {actual_dfdx.data(), actual_dfdy.data()};
  EXPECT_TRUE(problem.EvaluateResidualBlock(residual_block_id,
                                            true,
                                            false,
                                            &actual_cost,
                                            actual_f.data(),
                                            jacobians));
}

}  // namespace ceres::internal
