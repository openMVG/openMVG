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

#include "ceres/sparse_cholesky.h"

#include "ceres/cxsparse.h"
#include "ceres/eigensparse.h"
#include "ceres/suitesparse.h"

namespace ceres {
namespace internal {

SparseCholesky* SparseCholesky::Create(
    SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type,
    OrderingType ordering_type) {
  switch (sparse_linear_algebra_library_type) {
    case SUITE_SPARSE:
#ifndef CERES_NO_SUITESPARSE
      return SuiteSparseCholesky::Create(ordering_type);
#else
      LOG(FATAL) << "Ceres was compiled without support for SuiteSparse.";
      return NULL;
#endif

    case EIGEN_SPARSE:
#ifdef CERES_USE_EIGEN_SPARSE
      return EigenSparseCholesky::Create(ordering_type);
#else
      LOG(FATAL) << "Ceres was compiled without support for "
                 << "Eigen's sparse Cholesky factorization routines.";
      return NULL;
#endif

    case CX_SPARSE:
#ifndef CERES_NO_CXSPARSE
      return CXSparseCholesky::Create(ordering_type);
#else
      LOG(FATAL) << "Ceres was compiled without support for CXSparse.";
      return NULL;
#endif

    default:
      LOG(FATAL) << "Unknown sparse linear algebra library type : "
                 << SparseLinearAlgebraLibraryTypeToString(
                        sparse_linear_algebra_library_type);
  }
  return NULL;
}

SparseCholesky::~SparseCholesky() {}

LinearSolverTerminationType SparseCholesky::FactorAndSolve(
    CompressedRowSparseMatrix* lhs,
    const double* rhs,
    double* solution,
    std::string* message) {
  LinearSolverTerminationType termination_type = Factorize(lhs, message);
  if (termination_type == LINEAR_SOLVER_SUCCESS) {
    termination_type = Solve(rhs, solution, message);
  }
  return termination_type;
}

CompressedRowSparseMatrix::StorageType StorageTypeForSparseLinearAlgebraLibrary(
    SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type) {
  if (sparse_linear_algebra_library_type == SUITE_SPARSE) {
    return CompressedRowSparseMatrix::UPPER_TRIANGULAR;
  }
  return CompressedRowSparseMatrix::LOWER_TRIANGULAR;
}

}  // namespace internal
}  // namespace ceres
