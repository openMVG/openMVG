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
// Author: fredp@google.com (Fred Pighin)
//
// TODO(sameeragarwal): More comprehensive testing with larger and
// more badly conditioned problem.

#include "ceres/conjugate_gradients_solver.h"

#include <memory>

#include "ceres/internal/eigen.h"
#include "ceres/linear_solver.h"
#include "ceres/preconditioner.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "gtest/gtest.h"

namespace ceres::internal {

TEST(ConjugateGradientTest, Solves3x3IdentitySystem) {
  double diagonal[] = {1.0, 1.0, 1.0};
  std::unique_ptr<TripletSparseMatrix> A(
      TripletSparseMatrix::CreateSparseDiagonalMatrix(diagonal, 3));
  Vector b(3);
  Vector x(3);

  b(0) = 1.0;
  b(1) = 2.0;
  b(2) = 3.0;

  x(0) = 1;
  x(1) = 1;
  x(2) = 1;

  ConjugateGradientsSolverOptions cg_options;
  cg_options.min_num_iterations = 1;
  cg_options.max_num_iterations = 10;
  cg_options.residual_reset_period = 20;
  cg_options.q_tolerance = 0.0;
  cg_options.r_tolerance = 1e-9;

  Vector scratch[4];
  for (int i = 0; i < 4; ++i) {
    scratch[i] = Vector::Zero(A->num_cols());
  }

  IdentityPreconditioner identity(A->num_cols());
  LinearOperatorAdapter lhs(*A);
  LinearOperatorAdapter preconditioner(identity);
  Vector* scratch_array[4] = {
      &scratch[0], &scratch[1], &scratch[2], &scratch[3]};
  auto summary = ConjugateGradientsSolver(
      cg_options, lhs, b, preconditioner, scratch_array, x);

  EXPECT_EQ(summary.termination_type, LinearSolverTerminationType::SUCCESS);
  ASSERT_EQ(summary.num_iterations, 1);

  ASSERT_DOUBLE_EQ(1, x(0));
  ASSERT_DOUBLE_EQ(2, x(1));
  ASSERT_DOUBLE_EQ(3, x(2));
}

TEST(ConjuateGradientTest, Solves3x3SymmetricSystem) {
  std::unique_ptr<TripletSparseMatrix> A(new TripletSparseMatrix(3, 3, 9));
  Vector b(3);
  Vector x(3);

  //      | 2  -1  0|
  //  A = |-1   2 -1| is symmetric positive definite.
  //      | 0  -1  2|
  int* Ai = A->mutable_rows();
  int* Aj = A->mutable_cols();
  double* Ax = A->mutable_values();
  int counter = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Ai[counter] = i;
      Aj[counter] = j;
      ++counter;
    }
  }
  Ax[0] = 2.;
  Ax[1] = -1.;
  Ax[2] = 0;
  Ax[3] = -1.;
  Ax[4] = 2;
  Ax[5] = -1;
  Ax[6] = 0;
  Ax[7] = -1;
  Ax[8] = 2;
  A->set_num_nonzeros(9);

  b(0) = -1;
  b(1) = 0;
  b(2) = 3;

  x(0) = 1;
  x(1) = 1;
  x(2) = 1;

  ConjugateGradientsSolverOptions cg_options;
  cg_options.min_num_iterations = 1;
  cg_options.max_num_iterations = 10;
  cg_options.residual_reset_period = 20;
  cg_options.q_tolerance = 0.0;
  cg_options.r_tolerance = 1e-9;

  Vector scratch[4];
  for (int i = 0; i < 4; ++i) {
    scratch[i] = Vector::Zero(A->num_cols());
  }
  Vector* scratch_array[4] = {
      &scratch[0], &scratch[1], &scratch[2], &scratch[3]};
  IdentityPreconditioner identity(A->num_cols());
  LinearOperatorAdapter lhs(*A);
  LinearOperatorAdapter preconditioner(identity);

  auto summary = ConjugateGradientsSolver(
      cg_options, lhs, b, preconditioner, scratch_array, x);

  EXPECT_EQ(summary.termination_type, LinearSolverTerminationType::SUCCESS);

  ASSERT_DOUBLE_EQ(0, x(0));
  ASSERT_DOUBLE_EQ(1, x(1));
  ASSERT_DOUBLE_EQ(2, x(2));
}

}  // namespace ceres::internal
