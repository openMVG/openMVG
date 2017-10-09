// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2017 Google Inc. All rights reserved.
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

#include "ceres/sparse_normal_cholesky_solver.h"

#include <algorithm>
#include <cstring>
#include <ctime>

#include "ceres/block_sparse_matrix.h"
#include "ceres/inner_product_computer.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_solver.h"
#include "ceres/sparse_cholesky.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "ceres/wall_time.h"

namespace ceres {
namespace internal {

SparseNormalCholeskySolver::SparseNormalCholeskySolver(
    const LinearSolver::Options& options)
    : options_(options) {
  sparse_cholesky_.reset(
      SparseCholesky::Create(options_.sparse_linear_algebra_library_type,
                             options_.use_postordering ? AMD : NATURAL));
}

SparseNormalCholeskySolver::~SparseNormalCholeskySolver() {}

LinearSolver::Summary SparseNormalCholeskySolver::SolveImpl(
    BlockSparseMatrix* A,
    const double* b,
    const LinearSolver::PerSolveOptions& per_solve_options,
    double* x) {
  EventLogger event_logger("SparseNormalCholeskySolver::Solve");
  LinearSolver::Summary summary;
  summary.num_iterations = 1;
  summary.termination_type = LINEAR_SOLVER_SUCCESS;
  summary.message = "Success.";

  const int num_cols = A->num_cols();
  VectorRef(x, num_cols).setZero();
  A->LeftMultiply(b, x);
  event_logger.AddEvent("Compute RHS");

  if (per_solve_options.D != NULL) {
    // Temporarily append a diagonal block to the A matrix, but undo
    // it before returning the matrix to the user.
    scoped_ptr<BlockSparseMatrix> regularizer;
    regularizer.reset(BlockSparseMatrix::CreateDiagonalMatrix(
        per_solve_options.D, A->block_structure()->cols));
    event_logger.AddEvent("Diagonal");
    A->AppendRows(*regularizer);
    event_logger.AddEvent("Append");
  }
  event_logger.AddEvent("Append Rows");

  if (inner_product_computer_.get() == NULL) {
    inner_product_computer_.reset(
        InnerProductComputer::Create(*A, sparse_cholesky_->StorageType()));

    event_logger.AddEvent("InnerProductComputer::Create");
  }

  inner_product_computer_->Compute();
  event_logger.AddEvent("InnerProductComputer::Compute");

  // TODO(sameeragarwal):

  if (per_solve_options.D != NULL) {
    A->DeleteRowBlocks(A->block_structure()->cols.size());
  }
  summary.termination_type = sparse_cholesky_->FactorAndSolve(
      inner_product_computer_->mutable_result(), x, x, &summary.message);
  event_logger.AddEvent("Factor & Solve");
  return summary;
}

}  // namespace internal
}  // namespace ceres
