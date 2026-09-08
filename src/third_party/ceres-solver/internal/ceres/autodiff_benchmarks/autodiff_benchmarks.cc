// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2020 Google Inc. All rights reserved.
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
// Author: darius.rueckert@fau.de (Darius Rueckert)

#include <memory>
#include <random>
#include <utility>

#include "benchmark/benchmark.h"
#include "ceres/autodiff_benchmarks/brdf_cost_function.h"
#include "ceres/autodiff_benchmarks/constant_cost_function.h"
#include "ceres/autodiff_benchmarks/linear_cost_functions.h"
#include "ceres/autodiff_benchmarks/photometric_error.h"
#include "ceres/autodiff_benchmarks/relative_pose_error.h"
#include "ceres/autodiff_benchmarks/snavely_reprojection_error.h"
#include "ceres/ceres.h"

namespace ceres {

enum Dynamic { kNotDynamic, kDynamic };

// Transforms a static functor into a dynamic one.
template <typename CostFunctionType, int kNumParameterBlocks>
class ToDynamic {
 public:
  template <typename... _Args>
  explicit ToDynamic(_Args&&... __args)
      : cost_function_(std::forward<_Args>(__args)...) {}

  template <typename T>
  bool operator()(const T* const* parameters, T* residuals) const {
    return Apply(
        parameters, residuals, std::make_index_sequence<kNumParameterBlocks>());
  }

 private:
  template <typename T, size_t... Indices>
  bool Apply(const T* const* parameters,
             T* residuals,
             std::index_sequence<Indices...>) const {
    return cost_function_(parameters[Indices]..., residuals);
  }

  CostFunctionType cost_function_;
};

template <int kParameterBlockSize>
static void BM_ConstantAnalytic(benchmark::State& state) {
  constexpr int num_residuals = 1;
  std::array<double, kParameterBlockSize> parameters_values;
  std::iota(parameters_values.begin(), parameters_values.end(), 0);
  double* parameters[] = {parameters_values.data()};

  std::array<double, num_residuals> residuals;

  std::array<double, num_residuals * kParameterBlockSize> jacobian_values;
  double* jacobians[] = {jacobian_values.data()};

  std::unique_ptr<ceres::CostFunction> cost_function(
      new AnalyticConstantCostFunction<kParameterBlockSize>());

  for (auto _ : state) {
    cost_function->Evaluate(parameters, residuals.data(), jacobians);
  }
}

// Helpers for CostFunctionFactory.
template <typename DynamicCostFunctionType>
void AddParameterBlocks(DynamicCostFunctionType*) {}

template <int HeadN, int... TailNs, typename DynamicCostFunctionType>
void AddParameterBlocks(DynamicCostFunctionType* dynamic_function) {
  dynamic_function->AddParameterBlock(HeadN);
  AddParameterBlocks<TailNs...>(dynamic_function);
}

// Creates an autodiff cost function wrapping `CostFunctor`, with
// `kNumResiduals` residuals and parameter blocks with sized `Ns..`.
// Depending on `kIsDynamic`, either a static or dynamic cost function is
// created.
// `args` are forwarded to the `CostFunctor` constructor.
template <Dynamic kIsDynamic>
struct CostFunctionFactory {};

template <>
struct CostFunctionFactory<kNotDynamic> {
  template <typename CostFunctor,
            int kNumResiduals,
            int... Ns,
            typename... Args>
  static std::unique_ptr<ceres::CostFunction> Create(Args&&... args) {
    return std::make_unique<
        ceres::AutoDiffCostFunction<CostFunctor, kNumResiduals, Ns...>>(
        new CostFunctor(std::forward<Args>(args)...));
  }
};

template <>
struct CostFunctionFactory<kDynamic> {
  template <typename CostFunctor,
            int kNumResiduals,
            int... Ns,
            typename... Args>
  static std::unique_ptr<ceres::CostFunction> Create(Args&&... args) {
    constexpr const int kNumParameterBlocks = sizeof...(Ns);
    auto dynamic_function = std::make_unique<ceres::DynamicAutoDiffCostFunction<
        ToDynamic<CostFunctor, kNumParameterBlocks>>>(
        new ToDynamic<CostFunctor, kNumParameterBlocks>(
            std::forward<Args>(args)...));
    dynamic_function->SetNumResiduals(kNumResiduals);
    AddParameterBlocks<Ns...>(dynamic_function.get());
    return dynamic_function;
  }
};

template <int kParameterBlockSize, Dynamic kIsDynamic>
static void BM_ConstantAutodiff(benchmark::State& state) {
  constexpr int num_residuals = 1;
  std::array<double, kParameterBlockSize> parameters_values;
  std::iota(parameters_values.begin(), parameters_values.end(), 0);
  double* parameters[] = {parameters_values.data()};

  std::array<double, num_residuals> residuals;

  std::array<double, num_residuals * kParameterBlockSize> jacobian_values;
  double* jacobians[] = {jacobian_values.data()};

  std::unique_ptr<ceres::CostFunction> cost_function =
      CostFunctionFactory<kIsDynamic>::
          template Create<ConstantCostFunction<kParameterBlockSize>, 1, 1>();

  for (auto _ : state) {
    cost_function->Evaluate(parameters, residuals.data(), jacobians);
  }
}

BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 1);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 1, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 1, kDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 10);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 10, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 10, kDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 20);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 20, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 20, kDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 30);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 30, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 30, kDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 40);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 40, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 40, kDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 50);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 50, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 50, kDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAnalytic, 60);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 60, kNotDynamic);
BENCHMARK_TEMPLATE(BM_ConstantAutodiff, 60, kDynamic);

template <Dynamic kIsDynamic>
static void BM_Linear1AutoDiff(benchmark::State& state) {
  double parameter_block1[] = {1.};
  double* parameters[] = {parameter_block1};

  double jacobian1[1];
  double residuals[1];
  double* jacobians[] = {jacobian1};

  std::unique_ptr<ceres::CostFunction> cost_function = CostFunctionFactory<
      kIsDynamic>::template Create<Linear1CostFunction, 1, 1>();

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, residuals, state.range(0) ? jacobians : nullptr);
  }
}
BENCHMARK_TEMPLATE(BM_Linear1AutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_Linear1AutoDiff, kDynamic)->Arg(0)->Arg(1);

template <Dynamic kIsDynamic>
static void BM_Linear10AutoDiff(benchmark::State& state) {
  double parameter_block1[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  double* parameters[] = {parameter_block1};

  double jacobian1[10 * 10];
  double residuals[10];
  double* jacobians[] = {jacobian1};

  std::unique_ptr<ceres::CostFunction> cost_function = CostFunctionFactory<
      kIsDynamic>::template Create<Linear10CostFunction, 10, 10>();

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, residuals, state.range(0) ? jacobians : nullptr);
  }
}
BENCHMARK_TEMPLATE(BM_Linear10AutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_Linear10AutoDiff, kDynamic)->Arg(0)->Arg(1);

// From the NIST problem collection.
struct Rat43CostFunctor {
  Rat43CostFunctor(const double x, const double y) : x_(x), y_(y) {}

  template <typename T>
  inline bool operator()(const T* parameters, T* residuals) const {
    const T& b1 = parameters[0];
    const T& b2 = parameters[1];
    const T& b3 = parameters[2];
    const T& b4 = parameters[3];
    residuals[0] = b1 * pow(1.0 + exp(b2 - b3 * x_), -1.0 / b4) - y_;
    return true;
  }

  static constexpr int kNumParameterBlocks = 1;

 private:
  const double x_;
  const double y_;
};

template <Dynamic kIsDynamic>
static void BM_Rat43AutoDiff(benchmark::State& state) {
  double parameter_block1[] = {1., 2., 3., 4.};
  double* parameters[] = {parameter_block1};

  double jacobian1[] = {0.0, 0.0, 0.0, 0.0};
  double residuals;
  double* jacobians[] = {jacobian1};
  const double x = 0.2;
  const double y = 0.3;
  std::unique_ptr<ceres::CostFunction> cost_function =
      CostFunctionFactory<kIsDynamic>::template Create<Rat43CostFunctor, 1, 4>(
          x, y);

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, &residuals, state.range(0) ? jacobians : nullptr);
  }
}
BENCHMARK_TEMPLATE(BM_Rat43AutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_Rat43AutoDiff, kDynamic)->Arg(0)->Arg(1);

template <Dynamic kIsDynamic>
static void BM_SnavelyReprojectionAutoDiff(benchmark::State& state) {
  double parameter_block1[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
  double parameter_block2[] = {1., 2., 3.};
  double* parameters[] = {parameter_block1, parameter_block2};

  double jacobian1[2 * 9];
  double jacobian2[2 * 3];
  double residuals[2];
  double* jacobians[] = {jacobian1, jacobian2};

  const double x = 0.2;
  const double y = 0.3;
  std::unique_ptr<ceres::CostFunction> cost_function = CostFunctionFactory<
      kIsDynamic>::template Create<SnavelyReprojectionError, 2, 9, 3>(x, y);

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, residuals, state.range(0) ? jacobians : nullptr);
  }
}

BENCHMARK_TEMPLATE(BM_SnavelyReprojectionAutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_SnavelyReprojectionAutoDiff, kDynamic)->Arg(0)->Arg(1);

template <Dynamic kIsDynamic>
static void BM_PhotometricAutoDiff(benchmark::State& state) {
  constexpr int PATCH_SIZE = 8;

  using FunctorType = PhotometricError<PATCH_SIZE>;
  using ImageType = Eigen::Matrix<uint8_t, 128, 128, Eigen::RowMajor>;

  // Prepare parameter / residual / jacobian blocks.
  double parameter_block1[] = {1., 2., 3., 4., 5., 6., 7.};
  double parameter_block2[] = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1};
  double parameter_block3[] = {1.};
  double* parameters[] = {parameter_block1, parameter_block2, parameter_block3};

  Eigen::Map<Eigen::Quaterniond>(parameter_block1).normalize();
  Eigen::Map<Eigen::Quaterniond>(parameter_block2).normalize();

  double jacobian1[FunctorType::PATCH_SIZE * FunctorType::POSE_SIZE];
  double jacobian2[FunctorType::PATCH_SIZE * FunctorType::POSE_SIZE];
  double jacobian3[FunctorType::PATCH_SIZE * FunctorType::POINT_SIZE];
  double residuals[FunctorType::PATCH_SIZE];
  double* jacobians[] = {jacobian1, jacobian2, jacobian3};

  // Prepare data (fixed seed for repeatability).
  std::mt19937::result_type seed = 42;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::uniform_int_distribution<unsigned int> uniform0255(0, 255);

  FunctorType::Patch<double> intensities_host =
      FunctorType::Patch<double>::NullaryExpr(
          [&]() { return uniform0255(gen); });

  // Set bearing vector's z component to 1, i.e. pointing away from the camera,
  // to ensure they are (likely) in the domain of the projection function (given
  // a small rotation between host and target frame).
  FunctorType::PatchVectors<double> bearings_host =
      FunctorType::PatchVectors<double>::NullaryExpr(
          [&]() { return uniform01(gen); });
  bearings_host.row(2).array() = 1;
  bearings_host.colwise().normalize();

  ImageType image = ImageType::NullaryExpr(
      [&]() { return static_cast<uint8_t>(uniform0255(gen)); });
  FunctorType::Grid grid(image.data(), 0, image.rows(), 0, image.cols());
  FunctorType::Interpolator image_target(grid);

  FunctorType::Intrinsics intrinsics;
  intrinsics << 128, 128, 1, -1, 0.5, 0.5;

  std::unique_ptr<ceres::CostFunction> cost_function =
      CostFunctionFactory<kIsDynamic>::template Create<FunctorType,
                                                       FunctorType::PATCH_SIZE,
                                                       FunctorType::POSE_SIZE,
                                                       FunctorType::POSE_SIZE,
                                                       FunctorType::POINT_SIZE>(
          intensities_host, bearings_host, image_target, intrinsics);

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, residuals, state.range(0) ? jacobians : nullptr);
  }
}

BENCHMARK_TEMPLATE(BM_PhotometricAutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_PhotometricAutoDiff, kDynamic)->Arg(0)->Arg(1);

template <Dynamic kIsDynamic>
static void BM_RelativePoseAutoDiff(benchmark::State& state) {
  using FunctorType = RelativePoseError;

  double parameter_block1[] = {1., 2., 3., 4., 5., 6., 7.};
  double parameter_block2[] = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1};
  double* parameters[] = {parameter_block1, parameter_block2};

  Eigen::Map<Eigen::Quaterniond>(parameter_block1).normalize();
  Eigen::Map<Eigen::Quaterniond>(parameter_block2).normalize();

  double jacobian1[6 * 7];
  double jacobian2[6 * 7];
  double residuals[6];
  double* jacobians[] = {jacobian1, jacobian2};

  Eigen::Quaterniond q_i_j = Eigen::Quaterniond(1, 2, 3, 4).normalized();
  Eigen::Vector3d t_i_j(1, 2, 3);

  std::unique_ptr<ceres::CostFunction> cost_function =
      CostFunctionFactory<kIsDynamic>::template Create<FunctorType, 6, 7, 7>(
          q_i_j, t_i_j);

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, residuals, state.range(0) ? jacobians : nullptr);
  }
}

BENCHMARK_TEMPLATE(BM_RelativePoseAutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_RelativePoseAutoDiff, kDynamic)->Arg(0)->Arg(1);

template <Dynamic kIsDynamic>
static void BM_BrdfAutoDiff(benchmark::State& state) {
  using FunctorType = Brdf;

  double material[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  auto c = Eigen::Vector3d(0.1, 0.2, 0.3);
  auto n = Eigen::Vector3d(-0.1, 0.5, 0.2).normalized();
  auto v = Eigen::Vector3d(0.5, -0.2, 0.9).normalized();
  auto l = Eigen::Vector3d(-0.3, 0.4, -0.3).normalized();
  auto x = Eigen::Vector3d(0.5, 0.7, -0.1).normalized();
  auto y = Eigen::Vector3d(0.2, -0.2, -0.2).normalized();

  double* parameters[7] = {
      material, c.data(), n.data(), v.data(), l.data(), x.data(), y.data()};

  double jacobian[(10 + 6 * 3) * 3];
  double residuals[3];
  // clang-format off
  double* jacobians[7] = {
      jacobian + 0,      jacobian + 10 * 3, jacobian + 13 * 3,
      jacobian + 16 * 3, jacobian + 19 * 3, jacobian + 22 * 3,
      jacobian + 25 * 3,
  };
  // clang-format on

  std::unique_ptr<ceres::CostFunction> cost_function = CostFunctionFactory<
      kIsDynamic>::template Create<FunctorType, 3, 10, 3, 3, 3, 3, 3, 3>();

  for (auto _ : state) {
    cost_function->Evaluate(
        parameters, residuals, state.range(0) ? jacobians : nullptr);
  }
}

BENCHMARK_TEMPLATE(BM_BrdfAutoDiff, kNotDynamic)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(BM_BrdfAutoDiff, kDynamic)->Arg(0)->Arg(1);

}  // namespace ceres

BENCHMARK_MAIN();
