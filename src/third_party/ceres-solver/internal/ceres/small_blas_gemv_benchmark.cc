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
#include "ceres/small_blas.h"

namespace ceres {

// Benchmarking matrix-vector multiply routines and optimizing memory
// access requires that we make sure that they are not just sitting in
// the cache. So, as the benchmarking routine iterates, we need to
// multiply new/different matrice and vectors. Allocating/creating
// these objects in the benchmarking loop is too heavy duty, so we
// create them before hand and cycle through them in the
// benchmark. This class, given the size of the matrix creates such
// matrix and vector objects for use in the benchmark.
class MatrixVectorMultiplyData {
 public:
  MatrixVectorMultiplyData(int rows, int cols)
      : num_elements_(1000),
        rows_(rows),
        cols_(cols),
        a_(num_elements_ * rows, 1.001),
        b_(num_elements_ * rows * cols, 1.5),
        c_(num_elements_ * cols, 1.00003) {}

  int num_elements() const { return num_elements_; }
  double* GetA(int i) { return &a_[i * rows_]; }
  double* GetB(int i) { return &b_[i * rows_ * cols_]; }
  double* GetC(int i) { return &c_[i * cols_]; }

 private:
  const int num_elements_;
  const int rows_;
  const int cols_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;
};

// Helper function to generate the various matrix sizes for which we
// run the benchmark.
static void MatrixSizeArguments(benchmark::internal::Benchmark* benchmark) {
  std::vector<int> rows = {1, 2, 3, 4, 6, 8};
  std::vector<int> cols = {1, 2, 3, 4, 8, 12, 15};
  for (int r : rows) {
    for (int c : cols) {
      benchmark->Args({r, c});
    }
  }
}

static void BM_MatrixVectorMultiply(benchmark::State& state) {
  const int rows = state.range(0);
  const int cols = state.range(1);
  MatrixVectorMultiplyData data(rows, cols);
  const int num_elements = data.num_elements();
  int iter = 0;
  for (auto _ : state) {
    // A += B * C;
    internal::MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 1>(
        data.GetB(iter), rows, cols, data.GetC(iter), data.GetA(iter));
    iter = (iter + 1) % num_elements;
  }
}

BENCHMARK(BM_MatrixVectorMultiply)->Apply(MatrixSizeArguments);

static void BM_MatrixTransposeVectorMultiply(benchmark::State& state) {
  const int rows = state.range(0);
  const int cols = state.range(1);
  MatrixVectorMultiplyData data(cols, rows);
  const int num_elements = data.num_elements();
  int iter = 0;
  for (auto _ : state) {
    internal::MatrixTransposeVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 1>(
        data.GetB(iter), rows, cols, data.GetC(iter), data.GetA(iter));
    iter = (iter + 1) % num_elements;
  }
}

BENCHMARK(BM_MatrixTransposeVectorMultiply)->Apply(MatrixSizeArguments);

}  // namespace ceres

BENCHMARK_MAIN();
