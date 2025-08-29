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

#include <algorithm>
#include <memory>
#include <random>
#include <vector>

#include "Eigen/Dense"
#include "benchmark/benchmark.h"
#include "ceres/block_random_access_dense_matrix.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/block_structure.h"
#include "ceres/schur_eliminator.h"

namespace ceres::internal {

constexpr int kRowBlockSize = 2;
constexpr int kEBlockSize = 3;
constexpr int kFBlockSize = 6;

class BenchmarkData {
 public:
  explicit BenchmarkData(const int num_e_blocks) {
    auto* bs = new CompressedRowBlockStructure;
    bs->cols.resize(num_e_blocks + 1);
    int col_pos = 0;
    for (int i = 0; i < num_e_blocks; ++i) {
      bs->cols[i].position = col_pos;
      bs->cols[i].size = kEBlockSize;
      col_pos += kEBlockSize;
    }
    bs->cols.back().position = col_pos;
    bs->cols.back().size = kFBlockSize;

    bs->rows.resize(2 * num_e_blocks);
    int row_pos = 0;
    int cell_pos = 0;
    for (int i = 0; i < num_e_blocks; ++i) {
      {
        auto& row = bs->rows[2 * i];
        row.block.position = row_pos;
        row.block.size = kRowBlockSize;
        row_pos += kRowBlockSize;
        auto& cells = row.cells;
        cells.resize(2);
        cells[0].block_id = i;
        cells[0].position = cell_pos;
        cell_pos += kRowBlockSize * kEBlockSize;
        cells[1].block_id = num_e_blocks;
        cells[1].position = cell_pos;
        cell_pos += kRowBlockSize * kFBlockSize;
      }
      {
        auto& row = bs->rows[2 * i + 1];
        row.block.position = row_pos;
        row.block.size = kRowBlockSize;
        row_pos += kRowBlockSize;
        auto& cells = row.cells;
        cells.resize(1);
        cells[0].block_id = i;
        cells[0].position = cell_pos;
        cell_pos += kRowBlockSize * kEBlockSize;
      }
    }

    matrix_ = std::make_unique<BlockSparseMatrix>(bs);
    double* values = matrix_->mutable_values();
    std::generate_n(values, matrix_->num_nonzeros(), [this] {
      return standard_normal_(prng_);
    });

    b_.resize(matrix_->num_rows());
    b_.setRandom();

    std::vector<Block> blocks;
    blocks.emplace_back(kFBlockSize, 0);
    lhs_ = std::make_unique<BlockRandomAccessDenseMatrix>(blocks, &context_, 1);
    diagonal_.resize(matrix_->num_cols());
    diagonal_.setOnes();
    rhs_.resize(kFBlockSize);

    y_.resize(num_e_blocks * kEBlockSize);
    y_.setZero();
    z_.resize(kFBlockSize);
    z_.setOnes();
  }

  const BlockSparseMatrix& matrix() const { return *matrix_; }
  const Vector& b() const { return b_; }
  const Vector& diagonal() const { return diagonal_; }
  BlockRandomAccessDenseMatrix* mutable_lhs() { return lhs_.get(); }
  Vector* mutable_rhs() { return &rhs_; }
  Vector* mutable_y() { return &y_; }
  Vector* mutable_z() { return &z_; }

  ContextImpl* context() { return &context_; }

 private:
  ContextImpl context_;

  std::unique_ptr<BlockSparseMatrix> matrix_;
  Vector b_;
  std::unique_ptr<BlockRandomAccessDenseMatrix> lhs_;
  Vector rhs_;
  Vector diagonal_;
  Vector z_;
  Vector y_;
  std::mt19937 prng_;
  std::normal_distribution<> standard_normal_;
};

static void BM_SchurEliminatorEliminate(benchmark::State& state) {
  const int num_e_blocks = state.range(0);
  BenchmarkData data(num_e_blocks);

  LinearSolver::Options linear_solver_options;
  linear_solver_options.e_block_size = kEBlockSize;
  linear_solver_options.row_block_size = kRowBlockSize;
  linear_solver_options.f_block_size = kFBlockSize;
  linear_solver_options.context = data.context();
  std::unique_ptr<SchurEliminatorBase> eliminator(
      SchurEliminatorBase::Create(linear_solver_options));

  eliminator->Init(num_e_blocks, true, data.matrix().block_structure());
  for (auto _ : state) {
    eliminator->Eliminate(BlockSparseMatrixData(data.matrix()),
                          data.b().data(),
                          data.diagonal().data(),
                          data.mutable_lhs(),
                          data.mutable_rhs()->data());
  }
}

static void BM_SchurEliminatorBackSubstitute(benchmark::State& state) {
  const int num_e_blocks = state.range(0);
  BenchmarkData data(num_e_blocks);

  LinearSolver::Options linear_solver_options;
  linear_solver_options.e_block_size = kEBlockSize;
  linear_solver_options.row_block_size = kRowBlockSize;
  linear_solver_options.f_block_size = kFBlockSize;
  linear_solver_options.context = data.context();
  std::unique_ptr<SchurEliminatorBase> eliminator(
      SchurEliminatorBase::Create(linear_solver_options));

  eliminator->Init(num_e_blocks, true, data.matrix().block_structure());
  eliminator->Eliminate(BlockSparseMatrixData(data.matrix()),
                        data.b().data(),
                        data.diagonal().data(),
                        data.mutable_lhs(),
                        data.mutable_rhs()->data());
  for (auto _ : state) {
    eliminator->BackSubstitute(BlockSparseMatrixData(data.matrix()),
                               data.b().data(),
                               data.diagonal().data(),
                               data.mutable_z()->data(),
                               data.mutable_y()->data());
  }
}

static void BM_SchurEliminatorForOneFBlockEliminate(benchmark::State& state) {
  const int num_e_blocks = state.range(0);
  BenchmarkData data(num_e_blocks);
  SchurEliminatorForOneFBlock<2, 3, 6> eliminator;
  eliminator.Init(num_e_blocks, true, data.matrix().block_structure());
  for (auto _ : state) {
    eliminator.Eliminate(BlockSparseMatrixData(data.matrix()),
                         data.b().data(),
                         data.diagonal().data(),
                         data.mutable_lhs(),
                         data.mutable_rhs()->data());
  }
}

static void BM_SchurEliminatorForOneFBlockBackSubstitute(
    benchmark::State& state) {
  const int num_e_blocks = state.range(0);
  BenchmarkData data(num_e_blocks);
  SchurEliminatorForOneFBlock<2, 3, 6> eliminator;
  eliminator.Init(num_e_blocks, true, data.matrix().block_structure());
  eliminator.Eliminate(BlockSparseMatrixData(data.matrix()),
                       data.b().data(),
                       data.diagonal().data(),
                       data.mutable_lhs(),
                       data.mutable_rhs()->data());
  for (auto _ : state) {
    eliminator.BackSubstitute(BlockSparseMatrixData(data.matrix()),
                              data.b().data(),
                              data.diagonal().data(),
                              data.mutable_z()->data(),
                              data.mutable_y()->data());
  }
}

BENCHMARK(BM_SchurEliminatorEliminate)->Range(10, 10000);
BENCHMARK(BM_SchurEliminatorForOneFBlockEliminate)->Range(10, 10000);
BENCHMARK(BM_SchurEliminatorBackSubstitute)->Range(10, 10000);
BENCHMARK(BM_SchurEliminatorForOneFBlockBackSubstitute)->Range(10, 10000);

}  // namespace ceres::internal

BENCHMARK_MAIN();
