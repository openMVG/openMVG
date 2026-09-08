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

#include "ceres/iterative_refiner.h"

#include <utility>

#include "Eigen/Dense"
#include "ceres/dense_cholesky.h"
#include "ceres/internal/eigen.h"
#include "ceres/sparse_cholesky.h"
#include "ceres/sparse_matrix.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// Macros to help us define virtual methods which we do not expect to
// use/call in this test.
#define DO_NOT_CALL \
  { LOG(FATAL) << "DO NOT CALL"; }
#define DO_NOT_CALL_WITH_RETURN(x) \
  {                                \
    LOG(FATAL) << "DO NOT CALL";   \
    return x;                      \
  }

// A fake SparseMatrix, which uses an Eigen matrix to do the real work.
class FakeSparseMatrix : public SparseMatrix {
 public:
  explicit FakeSparseMatrix(Matrix m) : m_(std::move(m)) {}

  // y += Ax
  void RightMultiplyAndAccumulate(const double* x, double* y) const final {
    VectorRef(y, m_.cols()) += m_ * ConstVectorRef(x, m_.cols());
  }
  // y += A'x
  void LeftMultiplyAndAccumulate(const double* x, double* y) const final {
    // We will assume that this is a symmetric matrix.
    RightMultiplyAndAccumulate(x, y);
  }

  double* mutable_values() final { return m_.data(); }
  const double* values() const final { return m_.data(); }
  int num_rows() const final { return m_.cols(); }
  int num_cols() const final { return m_.cols(); }
  int num_nonzeros() const final { return m_.cols() * m_.cols(); }

  // The following methods are not needed for tests in this file.
  void SquaredColumnNorm(double* x) const final DO_NOT_CALL;
  void ScaleColumns(const double* scale) final DO_NOT_CALL;
  void SetZero() final DO_NOT_CALL;
  void ToDenseMatrix(Matrix* dense_matrix) const final DO_NOT_CALL;
  void ToTextFile(FILE* file) const final DO_NOT_CALL;

 private:
  Matrix m_;
};

// A fake SparseCholesky which uses Eigen's Cholesky factorization to
// do the real work. The template parameter allows us to work in
// doubles or floats, even though the source matrix is double.
template <typename Scalar>
class FakeSparseCholesky : public SparseCholesky {
 public:
  explicit FakeSparseCholesky(const Matrix& lhs) { lhs_ = lhs.cast<Scalar>(); }

  LinearSolverTerminationType Solve(const double* rhs_ptr,
                                    double* solution_ptr,
                                    std::string* message) final {
    const int num_cols = lhs_.cols();
    VectorRef solution(solution_ptr, num_cols);
    ConstVectorRef rhs(rhs_ptr, num_cols);
    auto llt = lhs_.llt();
    CHECK_EQ(llt.info(), Eigen::Success);
    solution = llt.solve(rhs.cast<Scalar>()).template cast<double>();
    return LinearSolverTerminationType::SUCCESS;
  }

  // The following methods are not needed for tests in this file.
  CompressedRowSparseMatrix::StorageType StorageType() const final
      DO_NOT_CALL_WITH_RETURN(
          CompressedRowSparseMatrix::StorageType::UPPER_TRIANGULAR);
  LinearSolverTerminationType Factorize(CompressedRowSparseMatrix* lhs,
                                        std::string* message) final
      DO_NOT_CALL_WITH_RETURN(LinearSolverTerminationType::FAILURE);

 private:
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> lhs_;
};

// A fake DenseCholesky which uses Eigen's Cholesky factorization to
// do the real work. The template parameter allows us to work in
// doubles or floats, even though the source matrix is double.
template <typename Scalar>
class FakeDenseCholesky : public DenseCholesky {
 public:
  explicit FakeDenseCholesky(const Matrix& lhs) { lhs_ = lhs.cast<Scalar>(); }

  LinearSolverTerminationType Solve(const double* rhs_ptr,
                                    double* solution_ptr,
                                    std::string* message) final {
    const int num_cols = lhs_.cols();
    VectorRef solution(solution_ptr, num_cols);
    ConstVectorRef rhs(rhs_ptr, num_cols);
    solution = lhs_.llt().solve(rhs.cast<Scalar>()).template cast<double>();
    return LinearSolverTerminationType::SUCCESS;
  }

  LinearSolverTerminationType Factorize(int num_cols,
                                        double* lhs,
                                        std::string* message) final
      DO_NOT_CALL_WITH_RETURN(LinearSolverTerminationType::FAILURE);

 private:
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> lhs_;
};

#undef DO_NOT_CALL
#undef DO_NOT_CALL_WITH_RETURN

class SparseIterativeRefinerTest : public ::testing::Test {
 public:
  void SetUp() override {
    num_cols_ = 5;
    max_num_iterations_ = 30;
    Matrix m(num_cols_, num_cols_);
    m.setRandom();
    lhs_ = m * m.transpose();
    solution_.resize(num_cols_);
    solution_.setRandom();
    rhs_ = lhs_ * solution_;
  };

 protected:
  int num_cols_;
  int max_num_iterations_;
  Matrix lhs_;
  Vector rhs_, solution_;
};

TEST_F(SparseIterativeRefinerTest,
       RandomSolutionWithExactFactorizationConverges) {
  FakeSparseMatrix lhs(lhs_);
  FakeSparseCholesky<double> sparse_cholesky(lhs_);
  SparseIterativeRefiner refiner(max_num_iterations_);
  Vector refined_solution(num_cols_);
  refined_solution.setRandom();
  refiner.Refine(lhs, rhs_.data(), &sparse_cholesky, refined_solution.data());
  EXPECT_NEAR((lhs_ * refined_solution - rhs_).norm(),
              0.0,
              std::numeric_limits<double>::epsilon() * 10);
}

TEST_F(SparseIterativeRefinerTest,
       RandomSolutionWithApproximationFactorizationConverges) {
  FakeSparseMatrix lhs(lhs_);
  // Use a single precision Cholesky factorization of the double
  // precision matrix. This will give us an approximate factorization.
  FakeSparseCholesky<float> sparse_cholesky(lhs_);
  SparseIterativeRefiner refiner(max_num_iterations_);
  Vector refined_solution(num_cols_);
  refined_solution.setRandom();
  refiner.Refine(lhs, rhs_.data(), &sparse_cholesky, refined_solution.data());
  EXPECT_NEAR((lhs_ * refined_solution - rhs_).norm(),
              0.0,
              std::numeric_limits<double>::epsilon() * 10);
}

class DenseIterativeRefinerTest : public ::testing::Test {
 public:
  void SetUp() override {
    num_cols_ = 5;
    max_num_iterations_ = 30;
    Matrix m(num_cols_, num_cols_);
    m.setRandom();
    lhs_ = m * m.transpose();
    solution_.resize(num_cols_);
    solution_.setRandom();
    rhs_ = lhs_ * solution_;
  };

 protected:
  int num_cols_;
  int max_num_iterations_;
  Matrix lhs_;
  Vector rhs_, solution_;
};

TEST_F(DenseIterativeRefinerTest,
       RandomSolutionWithExactFactorizationConverges) {
  Matrix lhs = lhs_;
  FakeDenseCholesky<double> dense_cholesky(lhs);
  DenseIterativeRefiner refiner(max_num_iterations_);
  Vector refined_solution(num_cols_);
  refined_solution.setRandom();
  refiner.Refine(lhs.cols(),
                 lhs.data(),
                 rhs_.data(),
                 &dense_cholesky,
                 refined_solution.data());
  EXPECT_NEAR((lhs_ * refined_solution - rhs_).norm(),
              0.0,
              std::numeric_limits<double>::epsilon() * 10);
}

TEST_F(DenseIterativeRefinerTest,
       RandomSolutionWithApproximationFactorizationConverges) {
  Matrix lhs = lhs_;
  // Use a single precision Cholesky factorization of the double
  // precision matrix. This will give us an approximate factorization.
  FakeDenseCholesky<float> dense_cholesky(lhs_);
  DenseIterativeRefiner refiner(max_num_iterations_);
  Vector refined_solution(num_cols_);
  refined_solution.setRandom();
  refiner.Refine(lhs.cols(),
                 lhs.data(),
                 rhs_.data(),
                 &dense_cholesky,
                 refined_solution.data());
  EXPECT_NEAR((lhs_ * refined_solution - rhs_).norm(),
              0.0,
              std::numeric_limits<double>::epsilon() * 10);
}

}  // namespace ceres::internal
