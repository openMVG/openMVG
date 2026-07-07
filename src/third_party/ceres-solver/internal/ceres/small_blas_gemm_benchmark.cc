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

#include <iostream>

#include "Eigen/Dense"
#include "benchmark/benchmark.h"
#include "ceres/small_blas.h"

namespace ceres::internal {

// Benchmarking matrix-matrix multiply routines and optimizing memory
// access requires that we make sure that they are not just sitting in
// the cache. So, as the benchmarking routine iterates, we need to
// multiply new/different matrice. Allocating/creating these objects
// in the benchmarking loop is too heavy duty, so we create them
// before hand and cycle through them in the benchmark. This class,
// given the size of the matrices creates such objects for use in the
// benchmark.
class MatrixMatrixMultiplyData {
 public:
  MatrixMatrixMultiplyData(
      int a_rows, int a_cols, int b_rows, int b_cols, int c_rows, int c_cols)
      : num_elements_(1000),
        a_size_(num_elements_ * a_rows * a_cols),
        b_size_(b_rows * b_cols),
        c_size_(c_rows * c_cols),
        a_(num_elements_ * a_size_, 1.00001),
        b_(num_elements_ * b_size_, 0.5),
        c_(num_elements_ * c_size_, -1.1) {}

  int num_elements() const { return num_elements_; }
  double* GetA(int i) { return &a_[i * a_size_]; }
  double* GetB(int i) { return &b_[i * b_size_]; }
  double* GetC(int i) { return &c_[i * c_size_]; }

 private:
  int num_elements_;
  int a_size_;
  int b_size_;
  int c_size_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;
};

#define GEMM_KIND_EQ 0
#define GEMM_KIND_ADD 1
#define GEMM_KIND_SUB -1

#define BENCHMARK_MM_FN(FN, M, N, K, NAME, MT, NT, KT)                      \
  void static BM_##FN##_##NAME##_##M##x##N##x##K(benchmark::State& state) { \
    const int b_rows = M;                                                   \
    const int b_cols = N;                                                   \
    const int c_rows = b_cols;                                              \
    const int c_cols = K;                                                   \
    const int a_rows = b_rows;                                              \
    const int a_cols = c_cols;                                              \
    MatrixMatrixMultiplyData data(                                          \
        a_rows, a_cols, b_rows, b_cols, c_rows, c_cols);                    \
    const int num_elements = data.num_elements();                           \
    int iter = 0;                                                           \
    for (auto _ : state) {                                                  \
      FN<MT, KT, KT, NT, GEMM_KIND_ADD>(data.GetB(iter),                    \
                                        b_rows,                             \
                                        b_cols,                             \
                                        data.GetC(iter),                    \
                                        c_rows,                             \
                                        c_cols,                             \
                                        data.GetA(iter),                    \
                                        512,                                \
                                        512,                                \
                                        a_rows,                             \
                                        a_cols);                            \
      iter = (iter + 1) % num_elements;                                     \
    }                                                                       \
  }                                                                         \
  BENCHMARK(BM_##FN##_##NAME##_##M##x##N##x##K);

#define BENCHMARK_STATIC_MM_FN(FN, M, N, K) \
  BENCHMARK_MM_FN(FN, M, N, K, Static, M, N, K)
#define BENCHMARK_DYNAMIC_MM_FN(FN, M, N, K) \
  BENCHMARK_MM_FN(                           \
      FN, M, N, K, Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic)

#define BENCHMARK_MTM_FN(FN, M, N, K, NAME, MT, NT, KT)                     \
  void static BM_##FN##_##NAME##_##M##x##N##x##K(benchmark::State& state) { \
    const int b_rows = M;                                                   \
    const int b_cols = N;                                                   \
    const int c_rows = b_rows;                                              \
    const int c_cols = K;                                                   \
    const int a_rows = b_cols;                                              \
    const int a_cols = c_cols;                                              \
    MatrixMatrixMultiplyData data(                                          \
        a_rows, a_cols, b_rows, b_cols, c_rows, c_cols);                    \
    const int num_elements = data.num_elements();                           \
    int iter = 0;                                                           \
    for (auto _ : state) {                                                  \
      FN<KT, MT, KT, NT, GEMM_KIND_ADD>(data.GetB(iter),                    \
                                        b_rows,                             \
                                        b_cols,                             \
                                        data.GetC(iter),                    \
                                        c_rows,                             \
                                        c_cols,                             \
                                        data.GetA(iter),                    \
                                        0,                                  \
                                        0,                                  \
                                        a_rows,                             \
                                        a_cols);                            \
      iter = (iter + 1) % num_elements;                                     \
    }                                                                       \
  }                                                                         \
  BENCHMARK(BM_##FN##_##NAME##_##M##x##N##x##K);

#define BENCHMARK_STATIC_MMT_FN(FN, M, N, K) \
  BENCHMARK_MTM_FN(FN, M, N, K, Static, M, N, K)
#define BENCHMARK_DYNAMIC_MMT_FN(FN, M, N, K) \
  BENCHMARK_MTM_FN(                           \
      FN, M, N, K, Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic)

BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 2, 3, 4)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 3, 3, 3)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 4, 4, 4)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 8, 8, 8)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 9, 9, 3)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 9, 3, 3)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyEigen, 3, 9, 9)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 2, 3, 4)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 3, 3, 3)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 4, 4, 4)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 8, 8, 8)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 9, 9, 3)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 9, 3, 3)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyEigen, 3, 9, 9)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 2, 3, 4)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 3, 3, 3)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 4, 4, 4)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 8, 8, 8)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 9, 9, 3)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 9, 3, 3)
BENCHMARK_STATIC_MM_FN(MatrixMatrixMultiplyNaive, 3, 9, 9)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 2, 3, 4)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 3, 3, 3)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 4, 4, 4)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 8, 8, 8)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 9, 9, 3)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 9, 3, 3)
BENCHMARK_DYNAMIC_MM_FN(MatrixMatrixMultiplyNaive, 3, 9, 9)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 2, 3, 4)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 3, 3, 3)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 4, 4, 4)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 8, 8, 8)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 9, 9, 3)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 9, 3, 3)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 3, 9, 9)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 2, 3, 4)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 3, 3, 3)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 4, 4, 4)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 8, 8, 8)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 9, 9, 3)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 9, 3, 3)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyEigen, 3, 9, 9)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 2, 3, 4)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 3, 3, 3)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 4, 4, 4)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 8, 8, 8)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 9, 9, 3)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 9, 3, 3)
BENCHMARK_STATIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 3, 9, 9)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 2, 3, 4)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 3, 3, 3)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 4, 4, 4)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 8, 8, 8)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 9, 9, 3)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 9, 3, 3)
BENCHMARK_DYNAMIC_MMT_FN(MatrixTransposeMatrixMultiplyNaive, 3, 9, 9)

#undef GEMM_KIND_EQ
#undef GEMM_KIND_ADD
#undef GEMM_KIND_SUB
#undef BENCHMARK_MM_FN
#undef BENCHMARK_STATIC_MM_FN
#undef BENCHMARK_DYNAMIC_MM_FN
#undef BENCHMARK_MTM_FN
#undef BENCHMARK_DYNAMIC_MMT_FN
#undef BENCHMARK_STATIC_MMT_FN

}  // namespace ceres::internal

BENCHMARK_MAIN();
