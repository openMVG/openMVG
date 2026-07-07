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

#include "ceres/schur_complement_solver.h"

#include <cstddef>
#include <memory>

#include "ceres/block_sparse_matrix.h"
#include "ceres/block_structure.h"
#include "ceres/casts.h"
#include "ceres/context_impl.h"
#include "ceres/detect_structure.h"
#include "ceres/linear_least_squares_problems.h"
#include "ceres/linear_solver.h"
#include "ceres/triplet_sparse_matrix.h"
#include "ceres/types.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

class SchurComplementSolverTest : public ::testing::Test {
 protected:
  void SetUpFromProblemId(int problem_id) {
    std::unique_ptr<LinearLeastSquaresProblem> problem =
        CreateLinearLeastSquaresProblemFromId(problem_id);

    CHECK(problem != nullptr);
    A.reset(down_cast<BlockSparseMatrix*>(problem->A.release()));
    b = std::move(problem->b);
    D = std::move(problem->D);

    num_cols = A->num_cols();
    num_rows = A->num_rows();
    num_eliminate_blocks = problem->num_eliminate_blocks;

    x.resize(num_cols);
    sol.resize(num_cols);
    sol_d.resize(num_cols);

    LinearSolver::Options options;
    options.type = DENSE_QR;
    ContextImpl context;
    options.context = &context;

    std::unique_ptr<LinearSolver> qr(LinearSolver::Create(options));

    TripletSparseMatrix triplet_A(
        A->num_rows(), A->num_cols(), A->num_nonzeros());
    A->ToTripletSparseMatrix(&triplet_A);

    // Gold standard solutions using dense QR factorization.
    DenseSparseMatrix dense_A(triplet_A);
    qr->Solve(&dense_A, b.get(), LinearSolver::PerSolveOptions(), sol.data());

    // Gold standard solution with appended diagonal.
    LinearSolver::PerSolveOptions per_solve_options;
    per_solve_options.D = D.get();
    qr->Solve(&dense_A, b.get(), per_solve_options, sol_d.data());
  }

  void ComputeAndCompareSolutions(
      int problem_id,
      bool regularization,
      ceres::LinearSolverType linear_solver_type,
      ceres::DenseLinearAlgebraLibraryType dense_linear_algebra_library_type,
      ceres::SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type,
      ceres::internal::OrderingType ordering_type) {
    SetUpFromProblemId(problem_id);
    LinearSolver::Options options;
    options.elimination_groups.push_back(num_eliminate_blocks);
    options.elimination_groups.push_back(A->block_structure()->cols.size() -
                                         num_eliminate_blocks);
    options.type = linear_solver_type;
    options.dense_linear_algebra_library_type =
        dense_linear_algebra_library_type;
    options.sparse_linear_algebra_library_type =
        sparse_linear_algebra_library_type;
    options.ordering_type = ordering_type;
    ContextImpl context;
    options.context = &context;
    DetectStructure(*A->block_structure(),
                    num_eliminate_blocks,
                    &options.row_block_size,
                    &options.e_block_size,
                    &options.f_block_size);

    std::unique_ptr<LinearSolver> solver(LinearSolver::Create(options));

    LinearSolver::PerSolveOptions per_solve_options;
    LinearSolver::Summary summary;
    if (regularization) {
      per_solve_options.D = D.get();
    }

    summary = solver->Solve(A.get(), b.get(), per_solve_options, x.data());
    EXPECT_EQ(summary.termination_type, LinearSolverTerminationType::SUCCESS);

    if (regularization) {
      ASSERT_NEAR((sol_d - x).norm() / num_cols, 0, 1e-10)
          << "Regularized Expected solution: " << sol_d.transpose()
          << " Actual solution: " << x.transpose();
    } else {
      ASSERT_NEAR((sol - x).norm() / num_cols, 0, 1e-10)
          << "Unregularized Expected solution: " << sol.transpose()
          << " Actual solution: " << x.transpose();
    }
  }

  int num_rows;
  int num_cols;
  int num_eliminate_blocks;

  std::unique_ptr<BlockSparseMatrix> A;
  std::unique_ptr<double[]> b;
  std::unique_ptr<double[]> D;
  Vector x;
  Vector sol;
  Vector sol_d;
};

// TODO(sameeragarwal): Refactor these using value parameterized tests.
// TODO(sameeragarwal): More extensive tests using random matrices.
TEST_F(SchurComplementSolverTest, DenseSchurWithEigenSmallProblem) {
  ComputeAndCompareSolutions(
      2, false, DENSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      2, true, DENSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
}

TEST_F(SchurComplementSolverTest, DenseSchurWithEigenLargeProblem) {
  ComputeAndCompareSolutions(
      3, false, DENSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      3, true, DENSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
}

TEST_F(SchurComplementSolverTest, DenseSchurWithEigenVaryingFBlockSize) {
  ComputeAndCompareSolutions(
      4, true, DENSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
}

#ifndef CERES_NO_LAPACK
TEST_F(SchurComplementSolverTest, DenseSchurWithLAPACKSmallProblem) {
  ComputeAndCompareSolutions(
      2, false, DENSE_SCHUR, LAPACK, SUITE_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      2, true, DENSE_SCHUR, LAPACK, SUITE_SPARSE, OrderingType::NATURAL);
}

TEST_F(SchurComplementSolverTest, DenseSchurWithLAPACKLargeProblem) {
  ComputeAndCompareSolutions(
      3, false, DENSE_SCHUR, LAPACK, SUITE_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      3, true, DENSE_SCHUR, LAPACK, SUITE_SPARSE, OrderingType::NATURAL);
}
#endif

#ifndef CERES_NO_SUITESPARSE
TEST_F(SchurComplementSolverTest,
       SparseSchurWithSuiteSparseSmallProblemNATURAL) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
}

TEST_F(SchurComplementSolverTest,
       SparseSchurWithSuiteSparseLargeProblemNATURAL) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NATURAL);
}

TEST_F(SchurComplementSolverTest, SparseSchurWithSuiteSparseSmallProblemAMD) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::AMD);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::AMD);
}

TEST_F(SchurComplementSolverTest, SparseSchurWithSuiteSparseLargeProblemAMD) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::AMD);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::AMD);
}

#ifndef CERES_NO_EIGEN_METIS
TEST_F(SchurComplementSolverTest,
       SparseSchurWithSuiteSparseSmallProblemNESDIS) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NESDIS);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NESDIS);
}
TEST_F(SchurComplementSolverTest,
       SparseSchurWithSuiteSparseLargeProblemNESDIS) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NESDIS);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, SUITE_SPARSE, OrderingType::NESDIS);
}
#endif  // CERES_NO_EIGEN_METIS
#endif  // CERES_NO_SUITESPARSE

#ifndef CERES_NO_ACCELERATE_SPARSE
TEST_F(SchurComplementSolverTest,
       SparseSchurWithAccelerateSparseSmallProblemAMD) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::AMD);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::AMD);
}

TEST_F(SchurComplementSolverTest,
       SparseSchurWithAccelerateSparseSmallProblemNESDIS) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::NESDIS);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::NESDIS);
}

TEST_F(SchurComplementSolverTest,
       SparseSchurWithAccelerateSparseLargeProblemAMD) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::AMD);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::AMD);
}

TEST_F(SchurComplementSolverTest,
       SparseSchurWithAccelerateSparseLargeProblemNESDIS) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::NESDIS);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, ACCELERATE_SPARSE, OrderingType::NESDIS);
}
#endif  // CERES_NO_ACCELERATE_SPARSE

#ifdef CERES_USE_EIGEN_SPARSE
TEST_F(SchurComplementSolverTest, SparseSchurWithEigenSparseSmallProblemAMD) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::AMD);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::AMD);
}

#ifndef CERES_NO_EIGEN_METIS
TEST_F(SchurComplementSolverTest,
       SparseSchurWithEigenSparseSmallProblemNESDIS) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NESDIS);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NESDIS);
}
#endif

TEST_F(SchurComplementSolverTest,
       SparseSchurWithEigenSparseSmallProblemNATURAL) {
  ComputeAndCompareSolutions(
      2, false, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      2, true, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NATURAL);
}

TEST_F(SchurComplementSolverTest, SparseSchurWithEigenSparseLargeProblemAMD) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::AMD);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::AMD);
}

#ifndef CERES_NO_EIGEN_METIS
TEST_F(SchurComplementSolverTest,
       SparseSchurWithEigenSparseLargeProblemNESDIS) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NESDIS);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NESDIS);
}
#endif

TEST_F(SchurComplementSolverTest,
       SparseSchurWithEigenSparseLargeProblemNATURAL) {
  ComputeAndCompareSolutions(
      3, false, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NATURAL);
  ComputeAndCompareSolutions(
      3, true, SPARSE_SCHUR, EIGEN, EIGEN_SPARSE, OrderingType::NATURAL);
}
#endif  // CERES_USE_EIGEN_SPARSE

}  // namespace ceres::internal
