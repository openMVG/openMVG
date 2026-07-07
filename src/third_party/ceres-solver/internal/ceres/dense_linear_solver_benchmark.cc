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
#include "ceres/context_impl.h"
#include "ceres/dense_sparse_matrix.h"
#include "ceres/internal/config.h"
#include "ceres/linear_solver.h"

namespace ceres::internal {

template <ceres::DenseLinearAlgebraLibraryType kLibraryType,
          ceres::LinearSolverType kSolverType>
static void BM_DenseSolver(benchmark::State& state) {
  const int num_rows = static_cast<int>(state.range(0));
  const int num_cols = static_cast<int>(state.range(1));
  DenseSparseMatrix jacobian(num_rows, num_cols);
  *jacobian.mutable_matrix() = Eigen::MatrixXd::Random(num_rows, num_cols);
  Eigen::VectorXd rhs = Eigen::VectorXd::Random(num_rows, 1);

  Eigen::VectorXd solution(num_cols);

  LinearSolver::Options options;
  options.type = kSolverType;
  options.dense_linear_algebra_library_type = kLibraryType;
  ContextImpl context;
  options.context = &context;
  auto solver = LinearSolver::Create(options);

  LinearSolver::PerSolveOptions per_solve_options;
  Eigen::VectorXd diagonal = Eigen::VectorXd::Ones(num_cols) * 100;
  per_solve_options.D = diagonal.data();
  for (auto _ : state) {
    solver->Solve(&jacobian, rhs.data(), per_solve_options, solution.data());
  }
}

// Some reasonable matrix sizes. I picked them out of thin air.
static void MatrixSizes(benchmark::internal::Benchmark* b) {
  // {num_rows, num_cols}
  b->Args({1, 1});
  b->Args({2, 1});
  b->Args({3, 1});
  b->Args({6, 2});
  b->Args({10, 3});
  b->Args({12, 4});
  b->Args({20, 5});
  b->Args({40, 5});
  b->Args({100, 10});
  b->Args({150, 15});
  b->Args({200, 16});
  b->Args({225, 18});
  b->Args({300, 20});
  b->Args({400, 20});
  b->Args({600, 22});
  b->Args({800, 25});
}

BENCHMARK_TEMPLATE2(BM_DenseSolver, ceres::EIGEN, ceres::DENSE_QR)
    ->Apply(MatrixSizes);
BENCHMARK_TEMPLATE2(BM_DenseSolver, ceres::EIGEN, ceres::DENSE_NORMAL_CHOLESKY)
    ->Apply(MatrixSizes);

#ifndef CERES_NO_LAPACK
BENCHMARK_TEMPLATE2(BM_DenseSolver, ceres::LAPACK, ceres::DENSE_QR)
    ->Apply(MatrixSizes);
BENCHMARK_TEMPLATE2(BM_DenseSolver, ceres::LAPACK, ceres::DENSE_NORMAL_CHOLESKY)
    ->Apply(MatrixSizes);
#endif  // CERES_NO_LAPACK

#ifndef CERES_NO_CUDA
BENCHMARK_TEMPLATE2(BM_DenseSolver, ceres::CUDA, ceres::DENSE_NORMAL_CHOLESKY)
    ->Apply(MatrixSizes);
BENCHMARK_TEMPLATE2(BM_DenseSolver, ceres::CUDA, ceres::DENSE_QR)
    ->Apply(MatrixSizes);
#endif  // CERES_NO_CUDA

}  // namespace ceres::internal

BENCHMARK_MAIN();
