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
// Authors: joydeepb@cs.utexas.edu (Joydeep Biswas)

#include <memory>
#include <random>
#include <string>

#include "Eigen/Dense"
#include "benchmark/benchmark.h"
#include "ceres/block_jacobi_preconditioner.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/context_impl.h"
#include "ceres/cuda_sparse_matrix.h"
#include "ceres/cuda_vector.h"
#include "ceres/fake_bundle_adjustment_jacobian.h"
#include "ceres/internal/config.h"
#include "ceres/internal/eigen.h"
#include "ceres/linear_solver.h"

#ifndef CERES_NO_CUDA
#include "cuda_runtime.h"
#endif

namespace ceres::internal {

constexpr int kNumCameras = 1000;
constexpr int kNumPoints = 10000;
constexpr int kCameraSize = 6;
constexpr int kPointSize = 3;
constexpr double kVisibility = 0.1;

constexpr int kNumRowBlocks = 100000;
constexpr int kNumColBlocks = 10000;
constexpr int kMinRowBlockSize = 1;
constexpr int kMaxRowBlockSize = 5;
constexpr int kMinColBlockSize = 1;
constexpr int kMaxColBlockSize = 15;
constexpr double kBlockDensity = 5.0 / kNumColBlocks;

static void BM_BlockSparseRightMultiplyAndAccumulateBA(
    benchmark::State& state) {
  const int num_threads = static_cast<int>(state.range(0));
  std::mt19937 prng;
  auto jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);

  ContextImpl context;
  context.EnsureMinimumThreads(num_threads);

  Vector x(jacobian->num_cols());
  Vector y(jacobian->num_rows());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(
        x.data(), y.data(), &context, num_threads);
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_BlockSparseRightMultiplyAndAccumulateBA)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_BlockSparseRightMultiplyAndAccumulateUnstructured(
    benchmark::State& state) {
  const int num_threads = static_cast<int>(state.range(0));
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_row_blocks = kNumRowBlocks;
  options.num_col_blocks = kNumColBlocks;
  options.min_row_block_size = kMinRowBlockSize;
  options.min_col_block_size = kMinColBlockSize;
  options.max_row_block_size = kMaxRowBlockSize;
  options.max_col_block_size = kMaxColBlockSize;
  options.block_density = kBlockDensity;
  std::mt19937 prng;

  auto jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);

  ContextImpl context;
  context.EnsureMinimumThreads(num_threads);

  Vector x(jacobian->num_cols());
  Vector y(jacobian->num_rows());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(
        x.data(), y.data(), &context, num_threads);
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_BlockSparseRightMultiplyAndAccumulateUnstructured)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_BlockSparseLeftMultiplyAndAccumulateBA(benchmark::State& state) {
  std::mt19937 prng;
  auto jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);
  Vector x(jacobian->num_rows());
  Vector y(jacobian->num_cols());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    jacobian->LeftMultiplyAndAccumulate(x.data(), y.data());
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_BlockSparseLeftMultiplyAndAccumulateBA);

static void BM_BlockSparseLeftMultiplyAndAccumulateUnstructured(
    benchmark::State& state) {
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_row_blocks = 100000;
  options.num_col_blocks = 10000;
  options.min_row_block_size = 1;
  options.min_col_block_size = 1;
  options.max_row_block_size = 10;
  options.max_col_block_size = 15;
  options.block_density = 5.0 / options.num_col_blocks;
  std::mt19937 prng;

  auto jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);
  Vector x(jacobian->num_rows());
  Vector y(jacobian->num_cols());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    jacobian->LeftMultiplyAndAccumulate(x.data(), y.data());
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_BlockSparseLeftMultiplyAndAccumulateUnstructured);

static void BM_CRSRightMultiplyAndAccumulateBA(benchmark::State& state) {
  const int num_threads = static_cast<int>(state.range(0));
  std::mt19937 prng;
  auto bsm_jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);

  auto jacobian = bsm_jacobian->ToCompressedRowSparseMatrix();

  ContextImpl context;
  context.EnsureMinimumThreads(num_threads);

  Vector x(jacobian->num_cols());
  Vector y(jacobian->num_rows());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(
        x.data(), y.data(), &context, num_threads);
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CRSRightMultiplyAndAccumulateBA)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_CRSRightMultiplyAndAccumulateUnstructured(
    benchmark::State& state) {
  const int num_threads = static_cast<int>(state.range(0));
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_row_blocks = kNumRowBlocks;
  options.num_col_blocks = kNumColBlocks;
  options.min_row_block_size = kMinRowBlockSize;
  options.min_col_block_size = kMinColBlockSize;
  options.max_row_block_size = kMaxRowBlockSize;
  options.max_col_block_size = kMaxColBlockSize;
  options.block_density = kBlockDensity;
  std::mt19937 prng;

  auto bsm_jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);
  auto jacobian = bsm_jacobian->ToCompressedRowSparseMatrix();

  ContextImpl context;
  context.EnsureMinimumThreads(num_threads);

  Vector x(jacobian->num_cols());
  Vector y(jacobian->num_rows());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(
        x.data(), y.data(), &context, num_threads);
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CRSRightMultiplyAndAccumulateUnstructured)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_CRSLeftMultiplyAndAccumulateBA(benchmark::State& state) {
  std::mt19937 prng;
  // Perform setup here
  auto bsm_jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);
  auto jacobian = bsm_jacobian->ToCompressedRowSparseMatrix();

  Vector x(jacobian->num_rows());
  Vector y(jacobian->num_cols());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    // This code gets timed
    jacobian->LeftMultiplyAndAccumulate(x.data(), y.data());
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CRSLeftMultiplyAndAccumulateBA);

static void BM_CRSLeftMultiplyAndAccumulateUnstructured(
    benchmark::State& state) {
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_row_blocks = kNumRowBlocks;
  options.num_col_blocks = kNumColBlocks;
  options.min_row_block_size = kMinRowBlockSize;
  options.min_col_block_size = kMinColBlockSize;
  options.max_row_block_size = kMaxRowBlockSize;
  options.max_col_block_size = kMaxColBlockSize;
  options.block_density = kBlockDensity;
  std::mt19937 prng;

  auto bsm_jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);
  auto jacobian = bsm_jacobian->ToCompressedRowSparseMatrix();

  Vector x(jacobian->num_rows());
  Vector y(jacobian->num_cols());
  x.setRandom();
  y.setRandom();
  double sum = 0;
  for (auto _ : state) {
    // This code gets timed
    jacobian->LeftMultiplyAndAccumulate(x.data(), y.data());
    sum += y.norm();
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CRSLeftMultiplyAndAccumulateUnstructured);

#ifndef CERES_NO_CUDA
static void BM_CudaRightMultiplyAndAccumulateBA(benchmark::State& state) {
  std::mt19937 prng;
  auto jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);
  ContextImpl context;
  std::string message;
  context.InitCuda(&message);
  auto jacobian_crs = jacobian->ToCompressedRowSparseMatrix();
  CudaSparseMatrix cuda_jacobian(&context, *jacobian_crs);
  CudaVector cuda_x(&context, 0);
  CudaVector cuda_y(&context, 0);

  Vector x(jacobian->num_cols());
  Vector y(jacobian->num_rows());
  x.setRandom();
  y.setRandom();

  cuda_x.CopyFromCpu(x);
  cuda_y.CopyFromCpu(y);
  double sum = 0;
  for (auto _ : state) {
    cuda_jacobian.RightMultiplyAndAccumulate(cuda_x, &cuda_y);
    sum += cuda_y.Norm();
    CHECK_EQ(cudaDeviceSynchronize(), cudaSuccess);
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CudaRightMultiplyAndAccumulateBA);

static void BM_CudaRightMultiplyAndAccumulateUnstructured(
    benchmark::State& state) {
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_row_blocks = kNumRowBlocks;
  options.num_col_blocks = kNumColBlocks;
  options.min_row_block_size = kMinRowBlockSize;
  options.min_col_block_size = kMinColBlockSize;
  options.max_row_block_size = kMaxRowBlockSize;
  options.max_col_block_size = kMaxColBlockSize;
  options.block_density = kBlockDensity;
  std::mt19937 prng;

  auto jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);
  ContextImpl context;
  std::string message;
  context.InitCuda(&message);
  auto jacobian_crs = jacobian->ToCompressedRowSparseMatrix();
  CudaSparseMatrix cuda_jacobian(&context, *jacobian_crs);
  CudaVector cuda_x(&context, 0);
  CudaVector cuda_y(&context, 0);

  Vector x(jacobian->num_cols());
  Vector y(jacobian->num_rows());
  x.setRandom();
  y.setRandom();

  cuda_x.CopyFromCpu(x);
  cuda_y.CopyFromCpu(y);
  double sum = 0;
  for (auto _ : state) {
    cuda_jacobian.RightMultiplyAndAccumulate(cuda_x, &cuda_y);
    sum += cuda_y.Norm();
    CHECK_EQ(cudaDeviceSynchronize(), cudaSuccess);
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CudaRightMultiplyAndAccumulateUnstructured);

static void BM_CudaLeftMultiplyAndAccumulateBA(benchmark::State& state) {
  std::mt19937 prng;
  auto jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);
  ContextImpl context;
  std::string message;
  context.InitCuda(&message);
  auto jacobian_crs = jacobian->ToCompressedRowSparseMatrix();
  CudaSparseMatrix cuda_jacobian(&context, *jacobian_crs);
  CudaVector cuda_x(&context, 0);
  CudaVector cuda_y(&context, 0);

  Vector x(jacobian->num_rows());
  Vector y(jacobian->num_cols());
  x.setRandom();
  y.setRandom();

  cuda_x.CopyFromCpu(x);
  cuda_y.CopyFromCpu(y);
  double sum = 0;
  for (auto _ : state) {
    cuda_jacobian.LeftMultiplyAndAccumulate(cuda_x, &cuda_y);
    sum += cuda_y.Norm();
    CHECK_EQ(cudaDeviceSynchronize(), cudaSuccess);
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CudaLeftMultiplyAndAccumulateBA);

static void BM_CudaLeftMultiplyAndAccumulateUnstructured(
    benchmark::State& state) {
  BlockSparseMatrix::RandomMatrixOptions options;
  options.num_row_blocks = kNumRowBlocks;
  options.num_col_blocks = kNumColBlocks;
  options.min_row_block_size = kMinRowBlockSize;
  options.min_col_block_size = kMinColBlockSize;
  options.max_row_block_size = kMaxRowBlockSize;
  options.max_col_block_size = kMaxColBlockSize;
  options.block_density = kBlockDensity;
  std::mt19937 prng;

  auto jacobian = BlockSparseMatrix::CreateRandomMatrix(options, prng);
  ContextImpl context;
  std::string message;
  context.InitCuda(&message);
  auto jacobian_crs = jacobian->ToCompressedRowSparseMatrix();
  CudaSparseMatrix cuda_jacobian(&context, *jacobian_crs);
  CudaVector cuda_x(&context, 0);
  CudaVector cuda_y(&context, 0);

  Vector x(jacobian->num_rows());
  Vector y(jacobian->num_cols());
  x.setRandom();
  y.setRandom();

  cuda_x.CopyFromCpu(x);
  cuda_y.CopyFromCpu(y);
  double sum = 0;
  for (auto _ : state) {
    cuda_jacobian.LeftMultiplyAndAccumulate(cuda_x, &cuda_y);
    sum += cuda_y.Norm();
    CHECK_EQ(cudaDeviceSynchronize(), cudaSuccess);
  }
  CHECK_NE(sum, 0.0);
}

BENCHMARK(BM_CudaLeftMultiplyAndAccumulateUnstructured);

#endif

}  // namespace ceres::internal

BENCHMARK_MAIN();
