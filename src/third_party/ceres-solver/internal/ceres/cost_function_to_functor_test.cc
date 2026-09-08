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

#include "ceres/cost_function_to_functor.h"

#include <cstdint>
#include <memory>
#include <vector>

#include "ceres/autodiff_cost_function.h"
#include "ceres/dynamic_autodiff_cost_function.h"
#include "ceres/dynamic_cost_function_to_functor.h"
#include "gtest/gtest.h"

namespace ceres::internal {

const double kTolerance = 1e-18;

static void ExpectCostFunctionsAreEqual(
    const CostFunction& cost_function,
    const CostFunction& actual_cost_function) {
  EXPECT_EQ(cost_function.num_residuals(),
            actual_cost_function.num_residuals());
  const int num_residuals = cost_function.num_residuals();
  const std::vector<int32_t>& parameter_block_sizes =
      cost_function.parameter_block_sizes();
  const std::vector<int32_t>& actual_parameter_block_sizes =
      actual_cost_function.parameter_block_sizes();
  EXPECT_EQ(parameter_block_sizes.size(), actual_parameter_block_sizes.size());

  int num_parameters = 0;
  for (int i = 0; i < parameter_block_sizes.size(); ++i) {
    EXPECT_EQ(parameter_block_sizes[i], actual_parameter_block_sizes[i]);
    num_parameters += parameter_block_sizes[i];
  }

  std::unique_ptr<double[]> parameters(new double[num_parameters]);
  for (int i = 0; i < num_parameters; ++i) {
    parameters[i] = static_cast<double>(i) + 1.0;
  }

  std::unique_ptr<double[]> residuals(new double[num_residuals]);
  std::unique_ptr<double[]> jacobians(
      new double[num_parameters * num_residuals]);

  std::unique_ptr<double[]> actual_residuals(new double[num_residuals]);
  std::unique_ptr<double[]> actual_jacobians(
      new double[num_parameters * num_residuals]);

  std::unique_ptr<double*[]> parameter_blocks(
      new double*[parameter_block_sizes.size()]);
  std::unique_ptr<double*[]> jacobian_blocks(
      new double*[parameter_block_sizes.size()]);
  std::unique_ptr<double*[]> actual_jacobian_blocks(
      new double*[parameter_block_sizes.size()]);

  num_parameters = 0;
  for (int i = 0; i < parameter_block_sizes.size(); ++i) {
    parameter_blocks[i] = parameters.get() + num_parameters;
    jacobian_blocks[i] = jacobians.get() + num_parameters * num_residuals;
    actual_jacobian_blocks[i] =
        actual_jacobians.get() + num_parameters * num_residuals;
    num_parameters += parameter_block_sizes[i];
  }

  EXPECT_TRUE(
      cost_function.Evaluate(parameter_blocks.get(), residuals.get(), nullptr));
  EXPECT_TRUE(actual_cost_function.Evaluate(
      parameter_blocks.get(), actual_residuals.get(), nullptr));
  for (int i = 0; i < num_residuals; ++i) {
    EXPECT_NEAR(residuals[i], actual_residuals[i], kTolerance)
        << "residual id: " << i;
  }

  EXPECT_TRUE(cost_function.Evaluate(
      parameter_blocks.get(), residuals.get(), jacobian_blocks.get()));
  EXPECT_TRUE(actual_cost_function.Evaluate(parameter_blocks.get(),
                                            actual_residuals.get(),
                                            actual_jacobian_blocks.get()));
  for (int i = 0; i < num_residuals; ++i) {
    EXPECT_NEAR(residuals[i], actual_residuals[i], kTolerance)
        << "residual : " << i;
  }

  for (int i = 0; i < num_residuals * num_parameters; ++i) {
    EXPECT_NEAR(jacobians[i], actual_jacobians[i], kTolerance)
        << "jacobian : " << i << " " << jacobians[i] << " "
        << actual_jacobians[i];
  }
}

struct OneParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1, T* residuals) const {
    residuals[0] = x1[0] * x1[0];
    residuals[1] = x1[1] * x1[1];
    return true;
  }
};

struct TwoParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1, const T* x2, T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1];
    return true;
  }
};

struct ThreeParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1, const T* x2, const T* x3, T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1];
    return true;
  }
};

struct FourParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(
      const T* x1, const T* x2, const T* x3, const T* x4, T* residuals) const {
    residuals[0] =
        x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] + x4[0] * x4[0];
    residuals[1] =
        x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] + x4[1] * x4[1];
    return true;
  }
};

struct FiveParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1,
                  const T* x2,
                  const T* x3,
                  const T* x4,
                  const T* x5,
                  T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] +
                   x4[0] * x4[0] + x5[0] * x5[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] +
                   x4[1] * x4[1] + x5[1] * x5[1];
    return true;
  }
};

struct SixParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1,
                  const T* x2,
                  const T* x3,
                  const T* x4,
                  const T* x5,
                  const T* x6,
                  T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] +
                   x4[0] * x4[0] + x5[0] * x5[0] + x6[0] * x6[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] +
                   x4[1] * x4[1] + x5[1] * x5[1] + x6[1] * x6[1];
    return true;
  }
};

struct SevenParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1,
                  const T* x2,
                  const T* x3,
                  const T* x4,
                  const T* x5,
                  const T* x6,
                  const T* x7,
                  T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] +
                   x4[0] * x4[0] + x5[0] * x5[0] + x6[0] * x6[0] +
                   x7[0] * x7[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] +
                   x4[1] * x4[1] + x5[1] * x5[1] + x6[1] * x6[1] +
                   x7[1] * x7[1];
    return true;
  }
};

struct EightParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1,
                  const T* x2,
                  const T* x3,
                  const T* x4,
                  const T* x5,
                  const T* x6,
                  const T* x7,
                  const T* x8,
                  T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] +
                   x4[0] * x4[0] + x5[0] * x5[0] + x6[0] * x6[0] +
                   x7[0] * x7[0] + x8[0] * x8[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] +
                   x4[1] * x4[1] + x5[1] * x5[1] + x6[1] * x6[1] +
                   x7[1] * x7[1] + x8[1] * x8[1];
    return true;
  }
};

struct NineParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1,
                  const T* x2,
                  const T* x3,
                  const T* x4,
                  const T* x5,
                  const T* x6,
                  const T* x7,
                  const T* x8,
                  const T* x9,
                  T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] +
                   x4[0] * x4[0] + x5[0] * x5[0] + x6[0] * x6[0] +
                   x7[0] * x7[0] + x8[0] * x8[0] + x9[0] * x9[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] +
                   x4[1] * x4[1] + x5[1] * x5[1] + x6[1] * x6[1] +
                   x7[1] * x7[1] + x8[1] * x8[1] + x9[1] * x9[1];
    return true;
  }
};

struct TenParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(const T* x1,
                  const T* x2,
                  const T* x3,
                  const T* x4,
                  const T* x5,
                  const T* x6,
                  const T* x7,
                  const T* x8,
                  const T* x9,
                  const T* x10,
                  T* residuals) const {
    residuals[0] = x1[0] * x1[0] + x2[0] * x2[0] + x3[0] * x3[0] +
                   x4[0] * x4[0] + x5[0] * x5[0] + x6[0] * x6[0] +
                   x7[0] * x7[0] + x8[0] * x8[0] + x9[0] * x9[0] +
                   x10[0] * x10[0];
    residuals[1] = x1[1] * x1[1] + x2[1] * x2[1] + x3[1] * x3[1] +
                   x4[1] * x4[1] + x5[1] * x5[1] + x6[1] * x6[1] +
                   x7[1] * x7[1] + x8[1] * x8[1] + x9[1] * x9[1] +
                   x10[1] * x10[1];
    return true;
  }
};

class DynamicTwoParameterBlockFunctor {
 public:
  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const {
    for (int i = 0; i < 2; ++i) {
      residuals[0] = parameters[i][0] * parameters[i][0];
      residuals[1] = parameters[i][1] * parameters[i][1];
    }
    return true;
  }
};

// Check that AutoDiff(Functor1) == AutoDiff(CostToFunctor(AutoDiff(Functor1)))
#define TEST_BODY(Functor1)                                                    \
  TEST(CostFunctionToFunctor, Functor1) {                                      \
    using CostFunction1 =                                                      \
        AutoDiffCostFunction<Functor1, 2, PARAMETER_BLOCK_SIZES>;              \
    using FunctionToFunctor = CostFunctionToFunctor<2, PARAMETER_BLOCK_SIZES>; \
    using CostFunction2 =                                                      \
        AutoDiffCostFunction<FunctionToFunctor, 2, PARAMETER_BLOCK_SIZES>;     \
                                                                               \
    std::unique_ptr<CostFunction> cost_function(new CostFunction2(             \
        new FunctionToFunctor(new CostFunction1(new Functor1))));              \
                                                                               \
    std::unique_ptr<CostFunction> actual_cost_function(                        \
        new CostFunction1(new Functor1));                                      \
    ExpectCostFunctionsAreEqual(*cost_function, *actual_cost_function);        \
  }

#define PARAMETER_BLOCK_SIZES 2
TEST_BODY(OneParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2
TEST_BODY(TwoParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2
TEST_BODY(ThreeParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2
TEST_BODY(FourParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2, 2
TEST_BODY(FiveParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2, 2, 2
TEST_BODY(SixParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2, 2, 2, 2
TEST_BODY(SevenParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2, 2, 2, 2, 2
TEST_BODY(EightParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2, 2, 2, 2, 2, 2
TEST_BODY(NineParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#define PARAMETER_BLOCK_SIZES 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
TEST_BODY(TenParameterBlockFunctor)
#undef PARAMETER_BLOCK_SIZES

#undef TEST_BODY

TEST(CostFunctionToFunctor, DynamicNumberOfResiduals) {
  std::unique_ptr<CostFunction> cost_function(
      new AutoDiffCostFunction<CostFunctionToFunctor<ceres::DYNAMIC, 2, 2>,
                               ceres::DYNAMIC,
                               2,
                               2>(
          new CostFunctionToFunctor<ceres::DYNAMIC, 2, 2>(
              new AutoDiffCostFunction<TwoParameterBlockFunctor, 2, 2, 2>(
                  new TwoParameterBlockFunctor)),
          2));

  std::unique_ptr<CostFunction> actual_cost_function(
      new AutoDiffCostFunction<TwoParameterBlockFunctor, 2, 2, 2>(
          new TwoParameterBlockFunctor));
  ExpectCostFunctionsAreEqual(*cost_function, *actual_cost_function);
}

TEST(CostFunctionToFunctor, DynamicCostFunctionToFunctor) {
  auto* actual_cost_function(
      new DynamicAutoDiffCostFunction<DynamicTwoParameterBlockFunctor>(
          new DynamicTwoParameterBlockFunctor));
  actual_cost_function->AddParameterBlock(2);
  actual_cost_function->AddParameterBlock(2);
  actual_cost_function->SetNumResiduals(2);

  DynamicAutoDiffCostFunction<DynamicCostFunctionToFunctor> cost_function(
      new DynamicCostFunctionToFunctor(actual_cost_function));
  cost_function.AddParameterBlock(2);
  cost_function.AddParameterBlock(2);
  cost_function.SetNumResiduals(2);

  ExpectCostFunctionsAreEqual(cost_function, *actual_cost_function);
}

}  // namespace ceres::internal
