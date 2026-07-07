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

#include "Eigen/Dense"
#include "benchmark/benchmark.h"
#include "ceres/invert_psd_matrix.h"

namespace ceres::internal {

template <int kSize>
void BenchmarkFixedSizedInvertPSDMatrix(benchmark::State& state) {
  using MatrixType = typename EigenTypes<kSize, kSize>::Matrix;
  MatrixType input = MatrixType::Random();
  input += input.transpose() + MatrixType::Identity();

  MatrixType output;
  constexpr bool kAssumeFullRank = true;
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        output = InvertPSDMatrix<kSize>(kAssumeFullRank, input));
  }
}

BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 1);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 2);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 3);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 4);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 5);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 6);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 7);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 8);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 9);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 10);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 11);
BENCHMARK_TEMPLATE(BenchmarkFixedSizedInvertPSDMatrix, 12);

static void BenchmarkDynamicallyInvertPSDMatrix(benchmark::State& state) {
  using MatrixType =
      typename EigenTypes<Eigen::Dynamic, Eigen::Dynamic>::Matrix;
  const int size = static_cast<int>(state.range(0));
  MatrixType input = MatrixType::Random(size, size);
  input += input.transpose() + MatrixType::Identity(size, size);

  MatrixType output;
  constexpr bool kAssumeFullRank = true;
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        output = InvertPSDMatrix<Eigen::Dynamic>(kAssumeFullRank, input));
  }
}

BENCHMARK(BenchmarkDynamicallyInvertPSDMatrix)
    ->Apply([](benchmark::internal::Benchmark* benchmark) {
      for (int i = 1; i < 13; ++i) {
        benchmark->Args({i});
      }
    });

}  // namespace ceres::internal

BENCHMARK_MAIN();
