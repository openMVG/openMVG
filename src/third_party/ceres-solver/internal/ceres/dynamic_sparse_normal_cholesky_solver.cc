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

#include "ceres/dynamic_sparse_normal_cholesky_solver.h"

#include <algorithm>
#include <cstring>
#include <ctime>
#include <sstream>

#include "Eigen/SparseCore"
#include "ceres/compressed_row_sparse_matrix.h"
#include "ceres/cxsparse.h"
#include "ceres/internal/eigen.h"
#include "ceres/internal/scoped_ptr.h"
#include "ceres/linear_solver.h"
#include "ceres/suitesparse.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "ceres/wall_time.h"

#ifdef CERES_USE_EIGEN_SPARSE
#include "Eigen/SparseCholesky"
#endif

namespace ceres {
namespace internal {

DynamicSparseNormalCholeskySolver::DynamicSparseNormalCholeskySolver(
    const LinearSolver::Options& options)
    : options_(options) {}

LinearSolver::Summary DynamicSparseNormalCholeskySolver::SolveImpl(
    CompressedRowSparseMatrix* A,
    const double* b,
    const LinearSolver::PerSolveOptions& per_solve_options,
    double* x) {
  const int num_cols = A->num_cols();
  VectorRef(x, num_cols).setZero();
  A->LeftMultiply(b, x);

  if (per_solve_options.D != NULL) {
    // Temporarily append a diagonal block to the A matrix, but undo
    // it before returning the matrix to the user.
    scoped_ptr<CompressedRowSparseMatrix> regularizer;
    if (!A->col_blocks().empty()) {
      regularizer.reset(CompressedRowSparseMatrix::CreateBlockDiagonalMatrix(
          per_solve_options.D, A->col_blocks()));
    } else {
      regularizer.reset(
          new CompressedRowSparseMatrix(per_solve_options.D, num_cols));
    }
    A->AppendRows(*regularizer);
  }

  LinearSolver::Summary summary;
  switch (options_.sparse_linear_algebra_library_type) {
    case SUITE_SPARSE:
      summary = SolveImplUsingSuiteSparse(A, x);
      break;
    case CX_SPARSE:
      summary = SolveImplUsingCXSparse(A, x);
      break;
    case EIGEN_SPARSE:
      summary = SolveImplUsingEigen(A, x);
      break;
    default:
      LOG(FATAL) << "Unknown sparse linear algebra library : "
                 << options_.sparse_linear_algebra_library_type;
  }

  if (per_solve_options.D != NULL) {
    A->DeleteRows(num_cols);
  }

  return summary;
}

LinearSolver::Summary DynamicSparseNormalCholeskySolver::SolveImplUsingEigen(
    CompressedRowSparseMatrix* A, double* rhs_and_solution) {
#ifndef CERES_USE_EIGEN_SPARSE

  LinearSolver::Summary summary;
  summary.num_iterations = 0;
  summary.termination_type = LINEAR_SOLVER_FATAL_ERROR;
  summary.message =
      "SPARSE_NORMAL_CHOLESKY cannot be used with EIGEN_SPARSE "
      "because Ceres was not built with support for "
      "Eigen's SimplicialLDLT decomposition. "
      "This requires enabling building with -DEIGENSPARSE=ON.";
  return summary;

#else

  EventLogger event_logger("DynamicSparseNormalCholeskySolver::Eigen::Solve");

  Eigen::MappedSparseMatrix<double, Eigen::RowMajor> a(A->num_rows(),
                                                       A->num_cols(),
                                                       A->num_nonzeros(),
                                                       A->mutable_rows(),
                                                       A->mutable_cols(),
                                                       A->mutable_values());

  Eigen::SparseMatrix<double> lhs = a.transpose() * a;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

  LinearSolver::Summary summary;
  summary.num_iterations = 1;
  summary.termination_type = LINEAR_SOLVER_SUCCESS;
  summary.message = "Success.";

  solver.analyzePattern(lhs);
  if (VLOG_IS_ON(2)) {
    std::stringstream ss;
    solver.dumpMemory(ss);
    VLOG(2) << "Symbolic Analysis\n" << ss.str();
  }

  event_logger.AddEvent("Analyze");
  if (solver.info() != Eigen::Success) {
    summary.termination_type = LINEAR_SOLVER_FATAL_ERROR;
    summary.message = "Eigen failure. Unable to find symbolic factorization.";
    return summary;
  }

  solver.factorize(lhs);
  event_logger.AddEvent("Factorize");
  if (solver.info() != Eigen::Success) {
    summary.termination_type = LINEAR_SOLVER_FAILURE;
    summary.message = "Eigen failure. Unable to find numeric factorization.";
    return summary;
  }

  const Vector rhs = VectorRef(rhs_and_solution, lhs.cols());
  VectorRef(rhs_and_solution, lhs.cols()) = solver.solve(rhs);
  event_logger.AddEvent("Solve");
  if (solver.info() != Eigen::Success) {
    summary.termination_type = LINEAR_SOLVER_FAILURE;
    summary.message = "Eigen failure. Unable to do triangular solve.";
    return summary;
  }

  return summary;
#endif  // CERES_USE_EIGEN_SPARSE
}

LinearSolver::Summary DynamicSparseNormalCholeskySolver::SolveImplUsingCXSparse(
    CompressedRowSparseMatrix* A, double* rhs_and_solution) {
#ifdef CERES_NO_CXSPARSE

  LinearSolver::Summary summary;
  summary.num_iterations = 0;
  summary.termination_type = LINEAR_SOLVER_FATAL_ERROR;
  summary.message =
      "SPARSE_NORMAL_CHOLESKY cannot be used with CX_SPARSE "
      "because Ceres was not built with support for CXSparse. "
      "This requires enabling building with -DCXSPARSE=ON.";

  return summary;

#else
  EventLogger event_logger(
      "DynamicSparseNormalCholeskySolver::CXSparse::Solve");

  LinearSolver::Summary summary;
  summary.num_iterations = 1;
  summary.termination_type = LINEAR_SOLVER_SUCCESS;
  summary.message = "Success.";

  CXSparse cxsparse;

  // Wrap the augmented Jacobian in a compressed sparse column matrix.
  cs_di a_transpose = cxsparse.CreateSparseMatrixTransposeView(A);

  // Compute the normal equations. J'J delta = J'f and solve them
  // using a sparse Cholesky factorization. Notice that when compared
  // to SuiteSparse we have to explicitly compute the transpose of Jt,
  // and then the normal equations before they can be
  // factorized. CHOLMOD/SuiteSparse on the other hand can just work
  // off of Jt to compute the Cholesky factorization of the normal
  // equations.
  cs_di* a = cxsparse.TransposeMatrix(&a_transpose);
  cs_di* lhs = cxsparse.MatrixMatrixMultiply(&a_transpose, a);
  cxsparse.Free(a);
  event_logger.AddEvent("NormalEquations");

  if (!cxsparse.SolveCholesky(lhs, rhs_and_solution)) {
    summary.termination_type = LINEAR_SOLVER_FAILURE;
    summary.message = "CXSparse::SolveCholesky failed";
  }
  event_logger.AddEvent("Solve");

  cxsparse.Free(lhs);
  event_logger.AddEvent("TearDown");
  return summary;
#endif
}

LinearSolver::Summary
DynamicSparseNormalCholeskySolver::SolveImplUsingSuiteSparse(
    CompressedRowSparseMatrix* A, double* rhs_and_solution) {
#ifdef CERES_NO_SUITESPARSE

  LinearSolver::Summary summary;
  summary.num_iterations = 0;
  summary.termination_type = LINEAR_SOLVER_FATAL_ERROR;
  summary.message =
      "SPARSE_NORMAL_CHOLESKY cannot be used with SUITE_SPARSE "
      "because Ceres was not built with support for SuiteSparse. "
      "This requires enabling building with -DSUITESPARSE=ON.";
  return summary;

#else

  EventLogger event_logger(
      "DynamicSparseNormalCholeskySolver::SuiteSparse::Solve");
  LinearSolver::Summary summary;
  summary.termination_type = LINEAR_SOLVER_SUCCESS;
  summary.num_iterations = 1;
  summary.message = "Success.";

  SuiteSparse ss;
  const int num_cols = A->num_cols();
  cholmod_sparse lhs = ss.CreateSparseMatrixTransposeView(A);
  event_logger.AddEvent("Setup");
  cholmod_factor* factor = ss.AnalyzeCholesky(&lhs, &summary.message);
  event_logger.AddEvent("Analysis");

  if (factor == NULL) {
    summary.termination_type = LINEAR_SOLVER_FATAL_ERROR;
    return summary;
  }

  summary.termination_type = ss.Cholesky(&lhs, factor, &summary.message);
  if (summary.termination_type == LINEAR_SOLVER_SUCCESS) {
    cholmod_dense* rhs =
        ss.CreateDenseVector(rhs_and_solution, num_cols, num_cols);
    cholmod_dense* solution = ss.Solve(factor, rhs, &summary.message);
    event_logger.AddEvent("Solve");
    ss.Free(rhs);
    if (solution != NULL) {
      memcpy(
          rhs_and_solution, solution->x, num_cols * sizeof(*rhs_and_solution));
      ss.Free(solution);
    } else {
      summary.termination_type = LINEAR_SOLVER_FAILURE;
    }
  }

  ss.Free(factor);
  event_logger.AddEvent("Teardown");
  return summary;

#endif
}

}  // namespace internal
}  // namespace ceres
