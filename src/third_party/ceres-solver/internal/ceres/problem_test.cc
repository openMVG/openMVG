// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2010, 2011, 2012 Google Inc. All rights reserved.
// http://code.google.com/p/ceres-solver/
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

#include "gtest/gtest.h"
#include "ceres/cost_function.h"
#include "ceres/local_parameterization.h"
#include "ceres/sized_cost_function.h"
#include "ceres/internal/scoped_ptr.h"

namespace ceres {
namespace internal {

// The following three classes are for the purposes of defining
// function signatures. They have dummy Evaluate functions.

// Trivial cost function that accepts a single argument.
class UnaryCostFunction : public CostFunction {
 public:
  UnaryCostFunction(int num_residuals, int16 parameter_block_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block_size);
  }
  virtual ~UnaryCostFunction() {}

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 1;
    }
    return true;
  }
};

// Trivial cost function that accepts two arguments.
class BinaryCostFunction: public CostFunction {
 public:
  BinaryCostFunction(int num_residuals,
                     int16 parameter_block1_size,
                     int16 parameter_block2_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block1_size);
    mutable_parameter_block_sizes()->push_back(parameter_block2_size);
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 2;
    }
    return true;
  }
};

// Trivial cost function that accepts three arguments.
class TernaryCostFunction: public CostFunction {
 public:
  TernaryCostFunction(int num_residuals,
                      int16 parameter_block1_size,
                      int16 parameter_block2_size,
                      int16 parameter_block3_size) {
    set_num_residuals(num_residuals);
    mutable_parameter_block_sizes()->push_back(parameter_block1_size);
    mutable_parameter_block_sizes()->push_back(parameter_block2_size);
    mutable_parameter_block_sizes()->push_back(parameter_block3_size);
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    for (int i = 0; i < num_residuals(); ++i) {
      residuals[i] = 3;
    }
    return true;
  }
};

TEST(Problem, AddResidualWithNullCostFunctionDies) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);
  problem.AddParameterBlock(z, 5);

  EXPECT_DEATH_IF_SUPPORTED(problem.AddResidualBlock(NULL, NULL, x),
                            "'cost_function' Must be non NULL");
}

TEST(Problem, AddResidualWithIncorrectNumberOfParameterBlocksDies) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  problem.AddParameterBlock(y, 4);
  problem.AddParameterBlock(z, 5);

  // UnaryCostFunction takes only one parameter, but two are passed.
  EXPECT_DEATH_IF_SUPPORTED(
      problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x, y),
      "parameter_blocks.size()");
}

TEST(Problem, AddResidualWithDifferentSizesOnTheSameVariableDies) {
  double x[3];

  Problem problem;
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x);
  EXPECT_DEATH_IF_SUPPORTED(problem.AddResidualBlock(
                                new UnaryCostFunction(
                                    2, 4 /* 4 != 3 */), NULL, x),
                            "different block sizes");
}

TEST(Problem, AddResidualWithDuplicateParametersDies) {
  double x[3], z[5];

  Problem problem;
  EXPECT_DEATH_IF_SUPPORTED(problem.AddResidualBlock(
                                new BinaryCostFunction(2, 3, 3), NULL, x, x),
                            "Duplicate parameter blocks");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddResidualBlock(
                                new TernaryCostFunction(1, 5, 3, 5),
                                NULL, z, x, z),
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
  EXPECT_DEATH_IF_SUPPORTED(problem.AddResidualBlock(
      new BinaryCostFunction(2, 3, 4), NULL, x, z),
               "different block sizes");
}

TEST(Problem, AddResidualAddsDuplicatedParametersOnlyOnce) {
  double x[3], y[4], z[5];

  Problem problem;
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x);
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x);
  problem.AddResidualBlock(new UnaryCostFunction(2, 4), NULL, y);
  problem.AddResidualBlock(new UnaryCostFunction(2, 5), NULL, z);

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

static double *IntToPtr(int i) {
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
  problem.AddParameterBlock(IntToPtr(5),  5);  // x
  problem.AddParameterBlock(IntToPtr(13), 3);  // y

  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr( 4), 2),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr( 4), 3),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr( 4), 9),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr( 8), 3),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(12), 2),
                            "Aliasing detected");
  EXPECT_DEATH_IF_SUPPORTED(problem.AddParameterBlock(IntToPtr(14), 3),
                            "Aliasing detected");

  // These ones should work.
  problem.AddParameterBlock(IntToPtr( 2), 3);
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
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x);

  // ... even repeatedly.
  problem.AddParameterBlock(x, 3);
  problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x);

  // More parameters are fine.
  problem.AddParameterBlock(y, 4);
  problem.AddResidualBlock(new UnaryCostFunction(2, 4), NULL, y);

  EXPECT_EQ(2, problem.NumParameterBlocks());
  EXPECT_EQ(7, problem.NumParameters());
}

TEST(Problem, AddingParametersAndResidualsResultsInExpectedProblem) {
  double x[3], y[4], z[5], w[4];

  Problem problem;
  problem.AddParameterBlock(x, 3);
  EXPECT_EQ(1, problem.NumParameterBlocks());
  EXPECT_EQ(3, problem.NumParameters());

  problem.AddParameterBlock(y, 4);
  EXPECT_EQ(2, problem.NumParameterBlocks());
  EXPECT_EQ(7, problem.NumParameters());

  problem.AddParameterBlock(z, 5);
  EXPECT_EQ(3,  problem.NumParameterBlocks());
  EXPECT_EQ(12, problem.NumParameters());

  // Add a parameter that has a local parameterization.
  w[0] = 1.0; w[1] = 0.0; w[2] = 0.0; w[3] = 0.0;
  problem.AddParameterBlock(w, 4, new QuaternionParameterization);
  EXPECT_EQ(4,  problem.NumParameterBlocks());
  EXPECT_EQ(16, problem.NumParameters());

  problem.AddResidualBlock(new UnaryCostFunction(2, 3), NULL, x);
  problem.AddResidualBlock(new BinaryCostFunction(6, 5, 4) , NULL, z, y);
  problem.AddResidualBlock(new BinaryCostFunction(3, 3, 5), NULL, x, z);
  problem.AddResidualBlock(new BinaryCostFunction(7, 5, 3), NULL, z, x);
  problem.AddResidualBlock(new TernaryCostFunction(1, 5, 3, 4), NULL, z, x, y);

  const int total_residuals = 2 + 6 + 3 + 7 + 1;
  EXPECT_EQ(problem.NumResidualBlocks(), 5);
  EXPECT_EQ(problem.NumResiduals(), total_residuals);
}

class DestructorCountingCostFunction : public SizedCostFunction<3, 4, 5> {
 public:
  explicit DestructorCountingCostFunction(int *counter)
      : counter_(counter) {}

  virtual ~DestructorCountingCostFunction() {
    *counter_ += 1;
  }

  virtual bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const {
    return true;
  }

 private:
  int* counter_;
};

TEST(Problem, ReusedCostFunctionsAreOnlyDeletedOnce) {
  double y[4], z[5];
  int counter = 0;

  // Add a cost function multiple times and check to make sure that
  // the destructor on the cost function is only called once.
  {
    Problem problem;
    problem.AddParameterBlock(y, 4);
    problem.AddParameterBlock(z, 5);

    CostFunction* cost = new DestructorCountingCostFunction(&counter);
    problem.AddResidualBlock(cost, NULL, y, z);
    problem.AddResidualBlock(cost, NULL, y, z);
    problem.AddResidualBlock(cost, NULL, y, z);
  }

  // Check that the destructor was called only once.
  CHECK_EQ(counter, 1);
}

}  // namespace internal
}  // namespace ceres
