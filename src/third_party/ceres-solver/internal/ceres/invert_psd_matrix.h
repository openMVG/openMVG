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

#ifndef CERES_INTERNAL_INVERT_PSD_MATRIX_H_
#define CERES_INTERNAL_INVERT_PSD_MATRIX_H_

#include "ceres/internal/eigen.h"
#include "glog/logging.h"
#include "Eigen/Dense"

namespace ceres {
namespace internal {

// Helper routine to compute the inverse or pseudo-inverse of a
// symmetric positive semi-definite matrix.
//
// assume_full_rank controls whether a Cholesky factorization or an
// Singular Value Decomposition is used to compute the inverse and the
// pseudo-inverse respectively.
//
// The template parameter kSize can either be Eigen::Dynamic or a
// positive integer equal to the number of rows of m.
template <int kSize>
typename EigenTypes<kSize, kSize>::Matrix InvertPSDMatrix(
    const bool assume_full_rank,
    const typename EigenTypes<kSize, kSize>::Matrix& m) {
  const int size = m.rows();

  // If the matrix can be assumed to be full rank, then just use the
  // Cholesky factorization to invert it.
  if (assume_full_rank) {
    return m.template selfadjointView<Eigen::Upper>().llt().solve(
        Matrix::Identity(size, size));
  }

  Eigen::JacobiSVD<Matrix> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const double tolerance =
      std::numeric_limits<double>::epsilon() * size * svd.singularValues()(0);

  return svd.matrixV() *
         (svd.singularValues().array() > tolerance)
             .select(svd.singularValues().array().inverse(), 0)
             .matrix()
             .asDiagonal() *
         svd.matrixU().adjoint();
}

}  // namespace internal
}  // namespace ceres

#endif // CERES_INTERNAL_INVERT_PSD_MATRIX_H_
