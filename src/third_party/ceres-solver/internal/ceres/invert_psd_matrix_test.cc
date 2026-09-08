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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "ceres/invert_psd_matrix.h"

#include "ceres/internal/eigen.h"
#include "gtest/gtest.h"

namespace ceres::internal {

static constexpr bool kFullRank = true;
static constexpr bool kRankDeficient = false;

template <int kSize>
typename EigenTypes<kSize, kSize>::Matrix RandomPSDMatrixWithEigenValues(
    const typename EigenTypes<kSize>::Vector& eigenvalues) {
  typename EigenTypes<kSize, kSize>::Matrix m(eigenvalues.rows(),
                                              eigenvalues.rows());
  m.setRandom();
  Eigen::SelfAdjointEigenSolver<typename EigenTypes<kSize, kSize>::Matrix> es(
      m);
  return es.eigenvectors() * eigenvalues.asDiagonal() *
         es.eigenvectors().transpose();
}

TEST(InvertPSDMatrix, Identity3x3) {
  const Matrix m = Matrix::Identity(3, 3);
  const Matrix inverse_m = InvertPSDMatrix<3>(kFullRank, m);
  EXPECT_NEAR((inverse_m - m).norm() / m.norm(),
              0.0,
              std::numeric_limits<double>::epsilon());
}

TEST(InvertPSDMatrix, FullRank5x5) {
  EigenTypes<5>::Vector eigenvalues;
  eigenvalues.setRandom();
  eigenvalues = eigenvalues.array().abs().matrix();
  const Matrix m = RandomPSDMatrixWithEigenValues<5>(eigenvalues);
  const Matrix inverse_m = InvertPSDMatrix<5>(kFullRank, m);
  EXPECT_NEAR((m * inverse_m - Matrix::Identity(5, 5)).norm() / 5.0,
              0.0,
              10 * std::numeric_limits<double>::epsilon());
}

TEST(InvertPSDMatrix, RankDeficient5x5) {
  EigenTypes<5>::Vector eigenvalues;
  eigenvalues.setRandom();
  eigenvalues = eigenvalues.array().abs().matrix();
  eigenvalues(3) = 0.0;
  const Matrix m = RandomPSDMatrixWithEigenValues<5>(eigenvalues);
  const Matrix inverse_m = InvertPSDMatrix<5>(kRankDeficient, m);
  Matrix pseudo_identity = Matrix::Identity(5, 5);
  pseudo_identity(3, 3) = 0.0;
  EXPECT_NEAR((m * inverse_m * m - m).norm() / m.norm(),
              0.0,
              10 * std::numeric_limits<double>::epsilon());
}

TEST(InvertPSDMatrix, DynamicFullRank5x5) {
  EigenTypes<Eigen::Dynamic>::Vector eigenvalues(5);
  eigenvalues.setRandom();
  eigenvalues = eigenvalues.array().abs().matrix();
  const Matrix m = RandomPSDMatrixWithEigenValues<Eigen::Dynamic>(eigenvalues);
  const Matrix inverse_m = InvertPSDMatrix<Eigen::Dynamic>(kFullRank, m);
  EXPECT_NEAR((m * inverse_m - Matrix::Identity(5, 5)).norm() / 5.0,
              0.0,
              10 * std::numeric_limits<double>::epsilon());
}

TEST(InvertPSDMatrix, DynamicRankDeficient5x5) {
  EigenTypes<Eigen::Dynamic>::Vector eigenvalues(5);
  eigenvalues.setRandom();
  eigenvalues = eigenvalues.array().abs().matrix();
  eigenvalues(3) = 0.0;
  const Matrix m = RandomPSDMatrixWithEigenValues<Eigen::Dynamic>(eigenvalues);
  const Matrix inverse_m = InvertPSDMatrix<Eigen::Dynamic>(kRankDeficient, m);
  Matrix pseudo_identity = Matrix::Identity(5, 5);
  pseudo_identity(3, 3) = 0.0;
  EXPECT_NEAR((m * inverse_m * m - m).norm() / m.norm(),
              0.0,
              10 * std::numeric_limits<double>::epsilon());
}

}  // namespace ceres::internal
