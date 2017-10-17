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
//
// A simple C++ interface to the Eigen's Sparse Cholesky routines.

#ifndef CERES_INTERNAL_EIGENSPARSE_H_
#define CERES_INTERNAL_EIGENSPARSE_H_

// This include must come before any #ifndef check on Ceres compile options.
#include "ceres/internal/port.h"

#ifdef CERES_USE_EIGEN_SPARSE

#include <string>
#include "Eigen/SparseCore"
#include "ceres/linear_solver.h"
#include "ceres/sparse_cholesky.h"

namespace ceres {
namespace internal {

class EigenSparseCholesky : public SparseCholesky {
 public:
  // Factory
  static EigenSparseCholesky* Create(const OrderingType ordering_type);

  // SparseCholesky interface.
  virtual ~EigenSparseCholesky();
  virtual CompressedRowSparseMatrix::StorageType StorageType() const = 0;
  virtual LinearSolverTerminationType Factorize(
      const Eigen::SparseMatrix<double>& lhs, std::string* message) = 0;
  virtual LinearSolverTerminationType Solve(const double* rhs,
                                            double* solution,
                                            std::string* message) = 0;
};

}  // namespace internal
}  // namespace ceres

#endif  // CERES_USE_EIGEN_SPARSE

#endif  // CERES_INTERNAL_EIGENSPARSE_H_
