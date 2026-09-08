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

#include "ceres/normal_prior.h"

#include <algorithm>
#include <cstddef>
#include <random>

#include "ceres/internal/eigen.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(NormalPriorTest, ResidualAtRandomPosition) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  auto randu = [&distribution, &prng] { return distribution(prng); };
  for (int num_rows = 1; num_rows < 5; ++num_rows) {
    for (int num_cols = 1; num_cols < 5; ++num_cols) {
      Vector b(num_cols);
      b.setRandom();
      Matrix A(num_rows, num_cols);
      A.setRandom();

      auto* x = new double[num_cols];
      std::generate_n(x, num_cols, randu);

      auto* jacobian = new double[num_rows * num_cols];
      Vector residuals(num_rows);

      NormalPrior prior(A, b);
      prior.Evaluate(&x, residuals.data(), &jacobian);

      // Compare the norm of the residual
      double residual_diff_norm =
          (residuals - A * (VectorRef(x, num_cols) - b)).squaredNorm();
      EXPECT_NEAR(residual_diff_norm, 0, 1e-10);

      // Compare the jacobians
      MatrixRef J(jacobian, num_rows, num_cols);
      double jacobian_diff_norm = (J - A).norm();
      EXPECT_NEAR(jacobian_diff_norm, 0.0, 1e-10);

      delete[] x;
      delete[] jacobian;
    }
  }
}

TEST(NormalPriorTest, ResidualAtRandomPositionNullJacobians) {
  std::mt19937 prng;
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  auto randu = [&distribution, &prng] { return distribution(prng); };
  for (int num_rows = 1; num_rows < 5; ++num_rows) {
    for (int num_cols = 1; num_cols < 5; ++num_cols) {
      Vector b(num_cols);
      b.setRandom();
      Matrix A(num_rows, num_cols);
      A.setRandom();

      auto* x = new double[num_cols];
      std::generate_n(x, num_cols, randu);

      double* jacobians[1];
      jacobians[0] = nullptr;

      Vector residuals(num_rows);

      NormalPrior prior(A, b);
      prior.Evaluate(&x, residuals.data(), jacobians);

      // Compare the norm of the residual
      double residual_diff_norm =
          (residuals - A * (VectorRef(x, num_cols) - b)).squaredNorm();
      EXPECT_NEAR(residual_diff_norm, 0, 1e-10);

      prior.Evaluate(&x, residuals.data(), nullptr);
      // Compare the norm of the residual
      residual_diff_norm =
          (residuals - A * (VectorRef(x, num_cols) - b)).squaredNorm();
      EXPECT_NEAR(residual_diff_norm, 0, 1e-10);

      delete[] x;
    }
  }
}

}  // namespace internal
}  // namespace ceres
