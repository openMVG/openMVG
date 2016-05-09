// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
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
// Author: keir@google.com (Keir Mierle)
//         sameeragarwal@google.com (Sameer Agarwal)
//
// End-to-end tests for Ceres using Powell's function.

#include <cmath>
#include <cstdlib>

#include "ceres/autodiff_cost_function.h"
#include "ceres/problem.h"
#include "ceres/solver.h"
#include "ceres/test_util.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

// This class implements the SystemTestProblem interface and provides
// access to an implementation of Powell's singular function.
//
//   F = 1/2 (f1^2 + f2^2 + f3^2 + f4^2)
//
//   f1 = x1 + 10*x2;
//   f2 = sqrt(5) * (x3 - x4)
//   f3 = (x2 - 2*x3)^2
//   f4 = sqrt(10) * (x1 - x4)^2
//
// The starting values are x1 = 3, x2 = -1, x3 = 0, x4 = 1.
// The minimum is 0 at (x1, x2, x3, x4) = 0.
//
// From: Testing Unconstrained Optimization Software by Jorge J. More, Burton S.
// Garbow and Kenneth E. Hillstrom in ACM Transactions on Mathematical Software,
// Vol 7(1), March 1981.
class PowellsFunction {
 public:
  PowellsFunction() {
    x_[0] =  3.0;
    x_[1] = -1.0;
    x_[2] =  0.0;
    x_[3] =  1.0;

    problem_.AddResidualBlock(
        new AutoDiffCostFunction<F1, 1, 1, 1>(new F1), NULL, &x_[0], &x_[1]);
    problem_.AddResidualBlock(
        new AutoDiffCostFunction<F2, 1, 1, 1>(new F2), NULL, &x_[2], &x_[3]);
    problem_.AddResidualBlock(
        new AutoDiffCostFunction<F3, 1, 1, 1>(new F3), NULL, &x_[1], &x_[2]);
    problem_.AddResidualBlock(
        new AutoDiffCostFunction<F4, 1, 1, 1>(new F4), NULL, &x_[0], &x_[3]);

    // Settings for the reference solution.
    options_.linear_solver_type = ceres::DENSE_QR;
    options_.max_num_iterations = 10;
    options_.num_threads = 1;
  }

  Problem* mutable_problem() { return &problem_; }
  Solver::Options* mutable_solver_options() { return &options_; }

  static double kResidualTolerance;

 private:
  // Templated functions used for automatically differentiated cost
  // functions.
  class F1 {
   public:
    template <typename T> bool operator()(const T* const x1,
                                          const T* const x2,
                                          T* residual) const {
      // f1 = x1 + 10 * x2;
      *residual = *x1 + T(10.0) * *x2;
      return true;
    }
  };

  class F2 {
   public:
    template <typename T> bool operator()(const T* const x3,
                                          const T* const x4,
                                          T* residual) const {
      // f2 = sqrt(5) (x3 - x4)
      *residual = T(sqrt(5.0)) * (*x3 - *x4);
      return true;
    }
  };

  class F3 {
   public:
    template <typename T> bool operator()(const T* const x2,
                                          const T* const x4,
                                          T* residual) const {
      // f3 = (x2 - 2 x3)^2
      residual[0] = (x2[0] - T(2.0) * x4[0]) * (x2[0] - T(2.0) * x4[0]);
      return true;
    }
  };

  class F4 {
   public:
    template <typename T> bool operator()(const T* const x1,
                                          const T* const x4,
                                          T* residual) const {
      // f4 = sqrt(10) (x1 - x4)^2
      residual[0] = T(sqrt(10.0)) * (x1[0] - x4[0]) * (x1[0] - x4[0]);
      return true;
    }
  };

  double x_[4];
  Problem problem_;
  Solver::Options options_;
};

double PowellsFunction::kResidualTolerance = 1e-8;

typedef SystemTest<PowellsFunction> PowellTest;
const bool kAutomaticOrdering = true;

TEST_F(PowellTest, DenseQR) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(DENSE_QR, NO_SPARSE));
}

TEST_F(PowellTest, DenseNormalCholesky) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(DENSE_NORMAL_CHOLESKY));
}

TEST_F(PowellTest, DenseSchur) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(DENSE_SCHUR));
}

TEST_F(PowellTest, IterativeSchurWithJacobi) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(ITERATIVE_SCHUR, NO_SPARSE, kAutomaticOrdering, JACOBI));
}

#ifndef CERES_NO_SUITESPARSE
TEST_F(PowellTest, SparseNormalCholeskyUsingSuiteSparse) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(SPARSE_NORMAL_CHOLESKY, SUITE_SPARSE, kAutomaticOrdering));
}
#endif  // CERES_NO_SUITESPARSE

#ifndef CERES_NO_CXSPARSE
TEST_F(PowellTest, SparseNormalCholeskyUsingCXSparse) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(SPARSE_NORMAL_CHOLESKY, CX_SPARSE, kAutomaticOrdering));
}
#endif  // CERES_NO_CXSPARSE

#ifdef CERES_USE_EIGEN_SPARSE
TEST_F(PowellTest, SparseNormalCholeskyUsingEigenSparse) {
  RunSolverForConfigAndExpectResidualsMatch(
      SolverConfig(SPARSE_NORMAL_CHOLESKY, EIGEN_SPARSE, kAutomaticOrdering));
}
#endif  // CERES_USE_EIGEN_SPARSE

}  // namespace internal
}  // namespace ceres
