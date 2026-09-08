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
// Authors: dmitriy.korchemkin@gmail.com (Dmitriy Korchemkin)

#include <memory>
#include <random>
#include <string>
#include <vector>

#include "benchmark/benchmark.h"
#include "ceres/block_sparse_matrix.h"
#include "ceres/bundle_adjustment_test_util.h"
#include "ceres/cuda_block_sparse_crs_view.h"
#include "ceres/cuda_partitioned_block_sparse_crs_view.h"
#include "ceres/cuda_sparse_matrix.h"
#include "ceres/cuda_vector.h"
#include "ceres/evaluator.h"
#include "ceres/implicit_schur_complement.h"
#include "ceres/partitioned_matrix_view.h"
#include "ceres/power_series_expansion_preconditioner.h"
#include "ceres/preprocessor.h"
#include "ceres/problem.h"
#include "ceres/problem_impl.h"
#include "ceres/program.h"
#include "ceres/sparse_matrix.h"

namespace ceres::internal {

template <typename Derived, typename Base>
std::unique_ptr<Derived> downcast_unique_ptr(std::unique_ptr<Base>& base) {
  return std::unique_ptr<Derived>(dynamic_cast<Derived*>(base.release()));
}

// Benchmark library might invoke benchmark function multiple times.
// In order to save time required to parse BAL data, we ensure that
// each dataset is being loaded at most once.
// Each type of jacobians is also cached after first creation
struct BALData {
  using PartitionedView = PartitionedMatrixView<2, 3, 9>;
  explicit BALData(const std::string& path) {
    bal_problem = std::make_unique<BundleAdjustmentProblem>(path);
    CHECK(bal_problem != nullptr);

    auto problem_impl = bal_problem->mutable_problem()->mutable_impl();
    auto preprocessor = Preprocessor::Create(MinimizerType::TRUST_REGION);

    preprocessed_problem = std::make_unique<PreprocessedProblem>();
    Solver::Options options = bal_problem->options();
    options.linear_solver_type = ITERATIVE_SCHUR;
    CHECK(preprocessor->Preprocess(
        options, problem_impl, preprocessed_problem.get()));

    auto program = preprocessed_problem->reduced_program.get();

    parameters.resize(program->NumParameters());
    program->ParameterBlocksToStateVector(parameters.data());

    const int num_residuals = program->NumResiduals();
    b.resize(num_residuals);

    std::mt19937 rng;
    std::normal_distribution<double> rnorm;
    for (int i = 0; i < num_residuals; ++i) {
      b[i] = rnorm(rng);
    }

    const int num_parameters = program->NumParameters();
    D.resize(num_parameters);
    for (int i = 0; i < num_parameters; ++i) {
      D[i] = rnorm(rng);
    }
  }

  std::unique_ptr<BlockSparseMatrix> CreateBlockSparseJacobian(
      ContextImpl* context, bool sequential) {
    auto problem = bal_problem->mutable_problem();
    auto problem_impl = problem->mutable_impl();
    CHECK(problem_impl != nullptr);

    Evaluator::Options options;
    options.linear_solver_type = ITERATIVE_SCHUR;
    options.num_threads = 1;
    options.context = context;
    options.num_eliminate_blocks = bal_problem->num_points();

    std::string error;
    auto program = preprocessed_problem->reduced_program.get();
    auto evaluator = Evaluator::Create(options, program, &error);
    CHECK(evaluator != nullptr);

    auto jacobian = evaluator->CreateJacobian();
    auto block_sparse = downcast_unique_ptr<BlockSparseMatrix>(jacobian);
    CHECK(block_sparse != nullptr);

    if (sequential) {
      auto block_structure_sequential =
          std::make_unique<CompressedRowBlockStructure>(
              *block_sparse->block_structure());
      int num_nonzeros = 0;
      for (auto& row_block : block_structure_sequential->rows) {
        const int row_block_size = row_block.block.size;
        for (auto& cell : row_block.cells) {
          const int col_block_size =
              block_structure_sequential->cols[cell.block_id].size;
          cell.position = num_nonzeros;
          num_nonzeros += col_block_size * row_block_size;
        }
      }
      block_sparse = std::make_unique<BlockSparseMatrix>(
          block_structure_sequential.release(),
#ifndef CERES_NO_CUDA
          true
#else
          false
#endif
      );
    }

    std::mt19937 rng;
    std::normal_distribution<double> rnorm;
    const int nnz = block_sparse->num_nonzeros();
    auto values = block_sparse->mutable_values();
    for (int i = 0; i < nnz; ++i) {
      values[i] = rnorm(rng);
    }

    return block_sparse;
  }

  std::unique_ptr<CompressedRowSparseMatrix> CreateCompressedRowSparseJacobian(
      ContextImpl* context) {
    auto block_sparse = BlockSparseJacobian(context);
    return block_sparse->ToCompressedRowSparseMatrix();
  }

  const BlockSparseMatrix* BlockSparseJacobian(ContextImpl* context) {
    if (!block_sparse_jacobian) {
      block_sparse_jacobian = CreateBlockSparseJacobian(context, true);
    }
    return block_sparse_jacobian.get();
  }

  const BlockSparseMatrix* BlockSparseJacobianPartitioned(
      ContextImpl* context) {
    if (!block_sparse_jacobian_partitioned) {
      block_sparse_jacobian_partitioned =
          CreateBlockSparseJacobian(context, false);
    }
    return block_sparse_jacobian_partitioned.get();
  }

  const CompressedRowSparseMatrix* CompressedRowSparseJacobian(
      ContextImpl* context) {
    if (!crs_jacobian) {
      crs_jacobian = CreateCompressedRowSparseJacobian(context);
    }
    return crs_jacobian.get();
  }

  std::unique_ptr<PartitionedView> PartitionedMatrixViewJacobian(
      const LinearSolver::Options& options) {
    auto block_sparse = BlockSparseJacobianPartitioned(options.context);
    return std::make_unique<PartitionedView>(options, *block_sparse);
  }

  BlockSparseMatrix* BlockDiagonalEtE(const LinearSolver::Options& options) {
    if (!block_diagonal_ete) {
      auto partitioned_view = PartitionedMatrixViewJacobian(options);
      block_diagonal_ete = partitioned_view->CreateBlockDiagonalEtE();
    }
    return block_diagonal_ete.get();
  }

  BlockSparseMatrix* BlockDiagonalFtF(const LinearSolver::Options& options) {
    if (!block_diagonal_ftf) {
      auto partitioned_view = PartitionedMatrixViewJacobian(options);
      block_diagonal_ftf = partitioned_view->CreateBlockDiagonalFtF();
    }
    return block_diagonal_ftf.get();
  }

  const ImplicitSchurComplement* ImplicitSchurComplementWithoutDiagonal(
      const LinearSolver::Options& options) {
    auto block_sparse = BlockSparseJacobianPartitioned(options.context);
    implicit_schur_complement =
        std::make_unique<ImplicitSchurComplement>(options);
    implicit_schur_complement->Init(*block_sparse, nullptr, b.data());
    return implicit_schur_complement.get();
  }

  const ImplicitSchurComplement* ImplicitSchurComplementWithDiagonal(
      const LinearSolver::Options& options) {
    auto block_sparse = BlockSparseJacobianPartitioned(options.context);
    implicit_schur_complement_diag =
        std::make_unique<ImplicitSchurComplement>(options);
    implicit_schur_complement_diag->Init(*block_sparse, D.data(), b.data());
    return implicit_schur_complement_diag.get();
  }

  Vector parameters;
  Vector D;
  Vector b;
  std::unique_ptr<BundleAdjustmentProblem> bal_problem;
  std::unique_ptr<PreprocessedProblem> preprocessed_problem;
  std::unique_ptr<BlockSparseMatrix> block_sparse_jacobian_partitioned;
  std::unique_ptr<BlockSparseMatrix> block_sparse_jacobian;
  std::unique_ptr<CompressedRowSparseMatrix> crs_jacobian;
  std::unique_ptr<BlockSparseMatrix> block_diagonal_ete;
  std::unique_ptr<BlockSparseMatrix> block_diagonal_ftf;
  std::unique_ptr<ImplicitSchurComplement> implicit_schur_complement;
  std::unique_ptr<ImplicitSchurComplement> implicit_schur_complement_diag;
};

static void Residuals(benchmark::State& state,
                      BALData* data,
                      ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  Evaluator::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.num_threads = num_threads;
  options.context = context;
  options.num_eliminate_blocks = 0;

  std::string error;
  CHECK(data->preprocessed_problem != nullptr);
  auto program = data->preprocessed_problem->reduced_program.get();
  CHECK(program != nullptr);
  auto evaluator = Evaluator::Create(options, program, &error);
  CHECK(evaluator != nullptr);

  double cost = 0.;
  Vector residuals = Vector::Zero(program->NumResiduals());

  Evaluator::EvaluateOptions eval_options;
  for (auto _ : state) {
    CHECK(evaluator->Evaluate(eval_options,
                              data->parameters.data(),
                              &cost,
                              residuals.data(),
                              nullptr,
                              nullptr));
  }
}

static void ResidualsAndJacobian(benchmark::State& state,
                                 BALData* data,
                                 ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  Evaluator::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.num_threads = num_threads;
  options.context = context;
  options.num_eliminate_blocks = 0;

  std::string error;
  CHECK(data->preprocessed_problem != nullptr);
  auto program = data->preprocessed_problem->reduced_program.get();
  CHECK(program != nullptr);
  auto evaluator = Evaluator::Create(options, program, &error);
  CHECK(evaluator != nullptr);

  double cost = 0.;
  Vector residuals = Vector::Zero(program->NumResiduals());
  auto jacobian = evaluator->CreateJacobian();

  Evaluator::EvaluateOptions eval_options;
  for (auto _ : state) {
    CHECK(evaluator->Evaluate(eval_options,
                              data->parameters.data(),
                              &cost,
                              residuals.data(),
                              nullptr,
                              jacobian.get()));
  }
}

static void Plus(benchmark::State& state, BALData* data, ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  Evaluator::Options options;
  options.linear_solver_type = SPARSE_NORMAL_CHOLESKY;
  options.num_threads = num_threads;
  options.context = context;
  options.num_eliminate_blocks = 0;

  std::string error;
  CHECK(data->preprocessed_problem != nullptr);
  auto program = data->preprocessed_problem->reduced_program.get();
  CHECK(program != nullptr);
  auto evaluator = Evaluator::Create(options, program, &error);
  CHECK(evaluator != nullptr);

  Vector state_plus_delta = Vector::Zero(program->NumParameters());
  Vector delta = Vector::Random(program->NumEffectiveParameters());

  for (auto _ : state) {
    CHECK(evaluator->Plus(
        data->parameters.data(), delta.data(), state_plus_delta.data()));
  }
  CHECK_GT(state_plus_delta.squaredNorm(), 0.);
}

static void PSEPreconditioner(benchmark::State& state,
                              BALData* data,
                              ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;

  auto jacobian = data->ImplicitSchurComplementWithDiagonal(options);
  Preconditioner::Options preconditioner_options(options);

  PowerSeriesExpansionPreconditioner preconditioner(
      jacobian, 10, 0, preconditioner_options);

  Vector y = Vector::Zero(jacobian->num_cols());
  Vector x = Vector::Random(jacobian->num_cols());

  for (auto _ : state) {
    preconditioner.RightMultiplyAndAccumulate(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void PMVRightMultiplyAndAccumulateF(benchmark::State& state,
                                           BALData* data,
                                           ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);

  Vector y = Vector::Zero(jacobian->num_rows());
  Vector x = Vector::Random(jacobian->num_cols_f());

  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulateF(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void PMVLeftMultiplyAndAccumulateF(benchmark::State& state,
                                          BALData* data,
                                          ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);

  Vector y = Vector::Zero(jacobian->num_cols_f());
  Vector x = Vector::Random(jacobian->num_rows());

  for (auto _ : state) {
    jacobian->LeftMultiplyAndAccumulateF(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void PMVRightMultiplyAndAccumulateE(benchmark::State& state,
                                           BALData* data,
                                           ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);

  Vector y = Vector::Zero(jacobian->num_rows());
  Vector x = Vector::Random(jacobian->num_cols_e());

  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulateE(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void PMVLeftMultiplyAndAccumulateE(benchmark::State& state,
                                          BALData* data,
                                          ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);

  Vector y = Vector::Zero(jacobian->num_cols_e());
  Vector x = Vector::Random(jacobian->num_rows());

  for (auto _ : state) {
    jacobian->LeftMultiplyAndAccumulateE(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void PMVUpdateBlockDiagonalEtE(benchmark::State& state,
                                      BALData* data,
                                      ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);
  auto block_diagonal_ete = data->BlockDiagonalEtE(options);

  for (auto _ : state) {
    jacobian->UpdateBlockDiagonalEtE(block_diagonal_ete);
  }
}

static void PMVUpdateBlockDiagonalFtF(benchmark::State& state,
                                      BALData* data,
                                      ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);
  auto block_diagonal_ftf = data->BlockDiagonalFtF(options);

  for (auto _ : state) {
    jacobian->UpdateBlockDiagonalFtF(block_diagonal_ftf);
  }
}

static void ISCRightMultiplyNoDiag(benchmark::State& state,
                                   BALData* data,
                                   ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  auto jacobian = data->ImplicitSchurComplementWithoutDiagonal(options);

  Vector y = Vector::Zero(jacobian->num_rows());
  Vector x = Vector::Random(jacobian->num_cols());
  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void ISCRightMultiplyDiag(benchmark::State& state,
                                 BALData* data,
                                 ContextImpl* context) {
  LinearSolver::Options options;
  options.num_threads = static_cast<int>(state.range(0));
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;

  auto jacobian = data->ImplicitSchurComplementWithDiagonal(options);

  Vector y = Vector::Zero(jacobian->num_rows());
  Vector x = Vector::Random(jacobian->num_cols());
  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(x.data(), y.data());
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void JacobianToCRS(benchmark::State& state,
                          BALData* data,
                          ContextImpl* context) {
  auto jacobian = data->BlockSparseJacobian(context);

  std::unique_ptr<CompressedRowSparseMatrix> matrix;
  for (auto _ : state) {
    matrix = jacobian->ToCompressedRowSparseMatrix();
  }
  CHECK(matrix != nullptr);
}

#ifndef CERES_NO_CUDA
static void PMVRightMultiplyAndAccumulateFCuda(benchmark::State& state,
                                               BALData* data,
                                               ContextImpl* context) {
  LinearSolver::Options options;
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  options.num_threads = 1;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);
  auto underlying_matrix = data->BlockSparseJacobianPartitioned(context);
  CudaPartitionedBlockSparseCRSView view(
      *underlying_matrix, jacobian->num_col_blocks_e(), context);

  Vector x = Vector::Random(jacobian->num_cols_f());
  CudaVector cuda_x(context, x.size());
  CudaVector cuda_y(context, jacobian->num_rows());

  cuda_x.CopyFromCpu(x);
  cuda_y.SetZero();

  auto matrix = view.matrix_f();
  for (auto _ : state) {
    matrix->RightMultiplyAndAccumulate(cuda_x, &cuda_y);
  }
  CHECK_GT(cuda_y.Norm(), 0.);
}

static void PMVLeftMultiplyAndAccumulateFCuda(benchmark::State& state,
                                              BALData* data,
                                              ContextImpl* context) {
  LinearSolver::Options options;
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  options.num_threads = 1;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);
  auto underlying_matrix = data->BlockSparseJacobianPartitioned(context);
  CudaPartitionedBlockSparseCRSView view(
      *underlying_matrix, jacobian->num_col_blocks_e(), context);

  Vector x = Vector::Random(jacobian->num_rows());
  CudaVector cuda_x(context, x.size());
  CudaVector cuda_y(context, jacobian->num_cols_f());

  cuda_x.CopyFromCpu(x);
  cuda_y.SetZero();

  auto matrix = view.matrix_f();
  for (auto _ : state) {
    matrix->LeftMultiplyAndAccumulate(cuda_x, &cuda_y);
  }
  CHECK_GT(cuda_y.Norm(), 0.);
}

static void PMVRightMultiplyAndAccumulateECuda(benchmark::State& state,
                                               BALData* data,
                                               ContextImpl* context) {
  LinearSolver::Options options;
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  options.num_threads = 1;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);
  auto underlying_matrix = data->BlockSparseJacobianPartitioned(context);
  CudaPartitionedBlockSparseCRSView view(
      *underlying_matrix, jacobian->num_col_blocks_e(), context);

  Vector x = Vector::Random(jacobian->num_cols_e());
  CudaVector cuda_x(context, x.size());
  CudaVector cuda_y(context, jacobian->num_rows());

  cuda_x.CopyFromCpu(x);
  cuda_y.SetZero();

  auto matrix = view.matrix_e();
  for (auto _ : state) {
    matrix->RightMultiplyAndAccumulate(cuda_x, &cuda_y);
  }
  CHECK_GT(cuda_y.Norm(), 0.);
}

static void PMVLeftMultiplyAndAccumulateECuda(benchmark::State& state,
                                              BALData* data,
                                              ContextImpl* context) {
  LinearSolver::Options options;
  options.elimination_groups.push_back(data->bal_problem->num_points());
  options.context = context;
  options.num_threads = 1;
  auto jacobian = data->PartitionedMatrixViewJacobian(options);
  auto underlying_matrix = data->BlockSparseJacobianPartitioned(context);
  CudaPartitionedBlockSparseCRSView view(
      *underlying_matrix, jacobian->num_col_blocks_e(), context);

  Vector x = Vector::Random(jacobian->num_rows());
  CudaVector cuda_x(context, x.size());
  CudaVector cuda_y(context, jacobian->num_cols_e());

  cuda_x.CopyFromCpu(x);
  cuda_y.SetZero();

  auto matrix = view.matrix_e();
  for (auto _ : state) {
    matrix->LeftMultiplyAndAccumulate(cuda_x, &cuda_y);
  }
  CHECK_GT(cuda_y.Norm(), 0.);
}

// We want CudaBlockSparseCRSView to be not slower than explicit conversion to
// CRS on CPU
static void JacobianToCRSView(benchmark::State& state,
                              BALData* data,
                              ContextImpl* context) {
  auto jacobian = data->BlockSparseJacobian(context);

  std::unique_ptr<CudaBlockSparseCRSView> matrix;
  for (auto _ : state) {
    matrix = std::make_unique<CudaBlockSparseCRSView>(*jacobian, context);
  }
  CHECK(matrix != nullptr);
}
static void JacobianToCRSMatrix(benchmark::State& state,
                                BALData* data,
                                ContextImpl* context) {
  auto jacobian = data->BlockSparseJacobian(context);

  std::unique_ptr<CudaSparseMatrix> matrix;
  std::unique_ptr<CompressedRowSparseMatrix> matrix_cpu;
  for (auto _ : state) {
    matrix_cpu = jacobian->ToCompressedRowSparseMatrix();
    matrix = std::make_unique<CudaSparseMatrix>(context, *matrix_cpu);
  }
  CHECK(matrix != nullptr);
}
// Updating values in CudaBlockSparseCRSView should be +- as fast as just
// copying values (time spent in value permutation has to be hidden by PCIe
// transfer)
static void JacobianToCRSViewUpdate(benchmark::State& state,
                                    BALData* data,
                                    ContextImpl* context) {
  auto jacobian = data->BlockSparseJacobian(context);

  auto matrix = CudaBlockSparseCRSView(*jacobian, context);
  for (auto _ : state) {
    matrix.UpdateValues(*jacobian);
  }
}
static void JacobianToCRSMatrixUpdate(benchmark::State& state,
                                      BALData* data,
                                      ContextImpl* context) {
  auto jacobian = data->BlockSparseJacobian(context);

  auto matrix_cpu = jacobian->ToCompressedRowSparseMatrix();
  auto matrix = std::make_unique<CudaSparseMatrix>(context, *matrix_cpu);
  for (auto _ : state) {
    CHECK_EQ(cudaSuccess,
             cudaMemcpy(matrix->mutable_values(),
                        matrix_cpu->values(),
                        matrix->num_nonzeros() * sizeof(double),
                        cudaMemcpyHostToDevice));
  }
}
#endif

static void JacobianSquaredColumnNorm(benchmark::State& state,
                                      BALData* data,
                                      ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  auto jacobian = data->BlockSparseJacobian(context);

  Vector x = Vector::Zero(jacobian->num_cols());

  for (auto _ : state) {
    jacobian->SquaredColumnNorm(x.data(), context, num_threads);
  }
  CHECK_GT(x.squaredNorm(), 0.);
}

static void JacobianScaleColumns(benchmark::State& state,
                                 BALData* data,
                                 ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  auto jacobian_const = data->BlockSparseJacobian(context);
  auto jacobian = const_cast<BlockSparseMatrix*>(jacobian_const);

  Vector x = Vector::Ones(jacobian->num_cols());

  for (auto _ : state) {
    jacobian->ScaleColumns(x.data(), context, num_threads);
  }
}

static void JacobianRightMultiplyAndAccumulate(benchmark::State& state,
                                               BALData* data,
                                               ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  auto jacobian = data->BlockSparseJacobian(context);

  Vector y = Vector::Zero(jacobian->num_rows());
  Vector x = Vector::Random(jacobian->num_cols());

  for (auto _ : state) {
    jacobian->RightMultiplyAndAccumulate(
        x.data(), y.data(), context, num_threads);
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

static void JacobianLeftMultiplyAndAccumulate(benchmark::State& state,
                                              BALData* data,
                                              ContextImpl* context) {
  const int num_threads = static_cast<int>(state.range(0));

  auto jacobian = data->BlockSparseJacobian(context);

  Vector y = Vector::Zero(jacobian->num_cols());
  Vector x = Vector::Random(jacobian->num_rows());

  for (auto _ : state) {
    jacobian->LeftMultiplyAndAccumulate(
        x.data(), y.data(), context, num_threads);
  }
  CHECK_GT(y.squaredNorm(), 0.);
}

#ifndef CERES_NO_CUDA
static void JacobianRightMultiplyAndAccumulateCuda(benchmark::State& state,
                                                   BALData* data,
                                                   ContextImpl* context) {
  auto crs_jacobian = data->CompressedRowSparseJacobian(context);
  CudaSparseMatrix cuda_jacobian(context, *crs_jacobian);
  CudaVector cuda_x(context, 0);
  CudaVector cuda_y(context, 0);

  Vector x(crs_jacobian->num_cols());
  Vector y(crs_jacobian->num_rows());
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

static void JacobianLeftMultiplyAndAccumulateCuda(benchmark::State& state,
                                                  BALData* data,
                                                  ContextImpl* context) {
  auto crs_jacobian = data->CompressedRowSparseJacobian(context);
  CudaSparseMatrix cuda_jacobian(context, *crs_jacobian);
  CudaVector cuda_x(context, 0);
  CudaVector cuda_y(context, 0);

  Vector x(crs_jacobian->num_rows());
  Vector y(crs_jacobian->num_cols());
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
#endif

}  // namespace ceres::internal

// Older versions of benchmark library might come without ::benchmark::Shutdown
// function. We provide an empty fallback variant of Shutdown function in
// order to support both older and newer versions
namespace benchmark_shutdown_fallback {
template <typename... Args>
void Shutdown(Args... args) {}
};  // namespace benchmark_shutdown_fallback

int main(int argc, char** argv) {
  ::benchmark::Initialize(&argc, argv);

  std::vector<std::unique_ptr<ceres::internal::BALData>> benchmark_data;
  if (argc == 1) {
    LOG(FATAL) << "No input datasets specified. Usage: " << argv[0]
               << " [benchmark flags] path_to_BAL_data_1.txt ... "
                  "path_to_BAL_data_N.txt";
    return -1;
  }

  ceres::internal::ContextImpl context;
  context.EnsureMinimumThreads(16);
#ifndef CERES_NO_CUDA
  std::string message;
  context.InitCuda(&message);
#endif

  for (int i = 1; i < argc; ++i) {
    const std::string path(argv[i]);
    const std::string name_residuals = "Residuals<" + path + ">";
    benchmark_data.emplace_back(
        std::make_unique<ceres::internal::BALData>(path));
    auto data = benchmark_data.back().get();
    ::benchmark::RegisterBenchmark(
        name_residuals.c_str(), ceres::internal::Residuals, data, &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_jacobians = "ResidualsAndJacobian<" + path + ">";
    ::benchmark::RegisterBenchmark(name_jacobians.c_str(),
                                   ceres::internal::ResidualsAndJacobian,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_plus = "Plus<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_plus.c_str(), ceres::internal::Plus, data, &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_right_product =
        "JacobianRightMultiplyAndAccumulate<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_right_product.c_str(),
        ceres::internal::JacobianRightMultiplyAndAccumulate,
        data,
        &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_right_product_partitioned_f =
        "PMVRightMultiplyAndAccumulateF<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_right_product_partitioned_f.c_str(),
        ceres::internal::PMVRightMultiplyAndAccumulateF,
        data,
        &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

#ifndef CERES_NO_CUDA
    const std::string name_right_product_partitioned_f_cuda =
        "PMVRightMultiplyAndAccumulateFCuda<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_right_product_partitioned_f_cuda.c_str(),
        ceres::internal::PMVRightMultiplyAndAccumulateFCuda,
        data,
        &context);
#endif

    const std::string name_right_product_partitioned_e =
        "PMVRightMultiplyAndAccumulateE<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_right_product_partitioned_e.c_str(),
        ceres::internal::PMVRightMultiplyAndAccumulateE,
        data,
        &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

#ifndef CERES_NO_CUDA
    const std::string name_right_product_partitioned_e_cuda =
        "PMVRightMultiplyAndAccumulateECuda<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_right_product_partitioned_e_cuda.c_str(),
        ceres::internal::PMVRightMultiplyAndAccumulateECuda,
        data,
        &context);
#endif

    const std::string name_update_block_diagonal_ftf =
        "PMVUpdateBlockDiagonalFtF<" + path + ">";
    ::benchmark::RegisterBenchmark(name_update_block_diagonal_ftf.c_str(),
                                   ceres::internal::PMVUpdateBlockDiagonalFtF,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_pse =
        "PSEPreconditionerRightMultiplyAndAccumulate<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_pse.c_str(), ceres::internal::PSEPreconditioner, data, &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_isc_no_diag =
        "ISCRightMultiplyAndAccumulate<" + path + ">";
    ::benchmark::RegisterBenchmark(name_isc_no_diag.c_str(),
                                   ceres::internal::ISCRightMultiplyNoDiag,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_update_block_diagonal_ete =
        "PMVUpdateBlockDiagonalEtE<" + path + ">";
    ::benchmark::RegisterBenchmark(name_update_block_diagonal_ete.c_str(),
                                   ceres::internal::PMVUpdateBlockDiagonalEtE,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);
    const std::string name_isc_diag =
        "ISCRightMultiplyAndAccumulateDiag<" + path + ">";
    ::benchmark::RegisterBenchmark(name_isc_diag.c_str(),
                                   ceres::internal::ISCRightMultiplyDiag,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

#ifndef CERES_NO_CUDA
    const std::string name_right_product_cuda =
        "JacobianRightMultiplyAndAccumulateCuda<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_right_product_cuda.c_str(),
        ceres::internal::JacobianRightMultiplyAndAccumulateCuda,
        data,
        &context)
        ->Arg(1);
#endif

    const std::string name_left_product =
        "JacobianLeftMultiplyAndAccumulate<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_left_product.c_str(),
        ceres::internal::JacobianLeftMultiplyAndAccumulate,
        data,
        &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_left_product_partitioned_f =
        "PMVLeftMultiplyAndAccumulateF<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_left_product_partitioned_f.c_str(),
        ceres::internal::PMVLeftMultiplyAndAccumulateF,
        data,
        &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

#ifndef CERES_NO_CUDA
    const std::string name_left_product_partitioned_f_cuda =
        "PMVLeftMultiplyAndAccumulateFCuda<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_left_product_partitioned_f_cuda.c_str(),
        ceres::internal::PMVLeftMultiplyAndAccumulateFCuda,
        data,
        &context);
#endif

    const std::string name_left_product_partitioned_e =
        "PMVLeftMultiplyAndAccumulateE<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_left_product_partitioned_e.c_str(),
        ceres::internal::PMVLeftMultiplyAndAccumulateE,
        data,
        &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

#ifndef CERES_NO_CUDA
    const std::string name_left_product_partitioned_e_cuda =
        "PMVLeftMultiplyAndAccumulateECuda<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_left_product_partitioned_e_cuda.c_str(),
        ceres::internal::PMVLeftMultiplyAndAccumulateECuda,
        data,
        &context);
#endif

#ifndef CERES_NO_CUDA
    const std::string name_left_product_cuda =
        "JacobianLeftMultiplyAndAccumulateCuda<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_left_product_cuda.c_str(),
        ceres::internal::JacobianLeftMultiplyAndAccumulateCuda,
        data,
        &context)
        ->Arg(1);
#endif

    const std::string name_squared_column_norm =
        "JacobianSquaredColumnNorm<" + path + ">";
    ::benchmark::RegisterBenchmark(name_squared_column_norm.c_str(),
                                   ceres::internal::JacobianSquaredColumnNorm,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_scale_columns = "JacobianScaleColumns<" + path + ">";
    ::benchmark::RegisterBenchmark(name_scale_columns.c_str(),
                                   ceres::internal::JacobianScaleColumns,
                                   data,
                                   &context)
        ->Arg(1)
        ->Arg(2)
        ->Arg(4)
        ->Arg(8)
        ->Arg(16);

    const std::string name_to_crs = "JacobianToCRS<" + path + ">";
    ::benchmark::RegisterBenchmark(
        name_to_crs.c_str(), ceres::internal::JacobianToCRS, data, &context);
#ifndef CERES_NO_CUDA
    const std::string name_to_crs_view = "JacobianToCRSView<" + path + ">";
    ::benchmark::RegisterBenchmark(name_to_crs_view.c_str(),
                                   ceres::internal::JacobianToCRSView,
                                   data,
                                   &context);
    const std::string name_to_crs_matrix = "JacobianToCRSMatrix<" + path + ">";
    ::benchmark::RegisterBenchmark(name_to_crs_matrix.c_str(),
                                   ceres::internal::JacobianToCRSMatrix,
                                   data,
                                   &context);
    const std::string name_to_crs_view_update =
        "JacobianToCRSViewUpdate<" + path + ">";
    ::benchmark::RegisterBenchmark(name_to_crs_view_update.c_str(),
                                   ceres::internal::JacobianToCRSViewUpdate,
                                   data,
                                   &context);
    const std::string name_to_crs_matrix_update =
        "JacobianToCRSMatrixUpdate<" + path + ">";
    ::benchmark::RegisterBenchmark(name_to_crs_matrix_update.c_str(),
                                   ceres::internal::JacobianToCRSMatrixUpdate,
                                   data,
                                   &context);
#endif
  }
  ::benchmark::RunSpecifiedBenchmarks();

  using namespace ::benchmark;
  using namespace benchmark_shutdown_fallback;
  Shutdown();
  return 0;
}
