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
// Authors: sameeragarwal@google.com (Sameer Agarwal)

#include <memory>
#include <random>

#include "Eigen/Dense"
#include "benchmark/benchmark.h"
#include "ceres/block_jacobi_preconditioner.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/fake_bundle_adjustment_jacobian.h"
#include "ceres/internal/config.h"
#include "ceres/internal/eigen.h"

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

static void BM_BlockSparseJacobiPreconditionerBA(benchmark::State& state) {
  std::mt19937 prng;
  auto jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);

  Preconditioner::Options preconditioner_options;
  ContextImpl context;
  preconditioner_options.context = &context;
  preconditioner_options.num_threads = static_cast<int>(state.range(0));
  context.EnsureMinimumThreads(preconditioner_options.num_threads);
  BlockSparseJacobiPreconditioner p(preconditioner_options, *jacobian);

  Vector d = Vector::Ones(jacobian->num_cols());
  for (auto _ : state) {
    p.Update(*jacobian, d.data());
  }
}

BENCHMARK(BM_BlockSparseJacobiPreconditionerBA)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_BlockCRSJacobiPreconditionerBA(benchmark::State& state) {
  std::mt19937 prng;
  auto jacobian = CreateFakeBundleAdjustmentJacobian(
      kNumCameras, kNumPoints, kCameraSize, kPointSize, kVisibility, prng);

  auto jacobian_crs = jacobian->ToCompressedRowSparseMatrix();
  Preconditioner::Options preconditioner_options;
  ContextImpl context;
  preconditioner_options.context = &context;
  preconditioner_options.num_threads = static_cast<int>(state.range(0));
  context.EnsureMinimumThreads(preconditioner_options.num_threads);
  BlockCRSJacobiPreconditioner p(preconditioner_options, *jacobian_crs);

  Vector d = Vector::Ones(jacobian_crs->num_cols());
  for (auto _ : state) {
    p.Update(*jacobian_crs, d.data());
  }
}

BENCHMARK(BM_BlockCRSJacobiPreconditionerBA)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_BlockSparseJacobiPreconditionerUnstructured(
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
  Preconditioner::Options preconditioner_options;
  ContextImpl context;
  preconditioner_options.context = &context;
  preconditioner_options.num_threads = static_cast<int>(state.range(0));
  context.EnsureMinimumThreads(preconditioner_options.num_threads);
  BlockSparseJacobiPreconditioner p(preconditioner_options, *jacobian);

  Vector d = Vector::Ones(jacobian->num_cols());
  for (auto _ : state) {
    p.Update(*jacobian, d.data());
  }
}

BENCHMARK(BM_BlockSparseJacobiPreconditionerUnstructured)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

static void BM_BlockCRSJacobiPreconditionerUnstructured(
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
  auto jacobian_crs = jacobian->ToCompressedRowSparseMatrix();
  Preconditioner::Options preconditioner_options;
  ContextImpl context;
  preconditioner_options.context = &context;
  preconditioner_options.num_threads = static_cast<int>(state.range(0));
  context.EnsureMinimumThreads(preconditioner_options.num_threads);
  BlockCRSJacobiPreconditioner p(preconditioner_options, *jacobian_crs);

  Vector d = Vector::Ones(jacobian_crs->num_cols());
  for (auto _ : state) {
    p.Update(*jacobian_crs, d.data());
  }
}
BENCHMARK(BM_BlockCRSJacobiPreconditionerUnstructured)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16);

}  // namespace ceres::internal

BENCHMARK_MAIN();
