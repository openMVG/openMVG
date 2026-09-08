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
// Author: keir@google.com (Keir Mierle)

#include "ceres/small_blas.h"

#include <limits>
#include <string>

#include "ceres/internal/eigen.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

const double kTolerance = 5.0 * std::numeric_limits<double>::epsilon();

// Static or dynamic problem types.
enum class DimType { Static, Dynamic };

// Constructs matrix functor type.
#define MATRIX_FUN_TY(FN)                                         \
  template <int kRowA,                                            \
            int kColA,                                            \
            int kRowB,                                            \
            int kColB,                                            \
            int kOperation,                                       \
            DimType kDimType>                                     \
  struct FN##Ty {                                                 \
    void operator()(const double* A,                              \
                    const int num_row_a,                          \
                    const int num_col_a,                          \
                    const double* B,                              \
                    const int num_row_b,                          \
                    const int num_col_b,                          \
                    double* C,                                    \
                    const int start_row_c,                        \
                    const int start_col_c,                        \
                    const int row_stride_c,                       \
                    const int col_stride_c) {                     \
      if (kDimType == DimType::Static) {                          \
        FN<kRowA, kColA, kRowB, kColB, kOperation>(A,             \
                                                   num_row_a,     \
                                                   num_col_a,     \
                                                   B,             \
                                                   num_row_b,     \
                                                   num_col_b,     \
                                                   C,             \
                                                   start_row_c,   \
                                                   start_col_c,   \
                                                   row_stride_c,  \
                                                   col_stride_c); \
      } else {                                                    \
        FN<Eigen::Dynamic,                                        \
           Eigen::Dynamic,                                        \
           Eigen::Dynamic,                                        \
           Eigen::Dynamic,                                        \
           kOperation>(A,                                         \
                       num_row_a,                                 \
                       num_col_a,                                 \
                       B,                                         \
                       num_row_b,                                 \
                       num_col_b,                                 \
                       C,                                         \
                       start_row_c,                               \
                       start_col_c,                               \
                       row_stride_c,                              \
                       col_stride_c);                             \
      }                                                           \
    }                                                             \
  };

MATRIX_FUN_TY(MatrixMatrixMultiply)
MATRIX_FUN_TY(MatrixMatrixMultiplyNaive)
MATRIX_FUN_TY(MatrixTransposeMatrixMultiply)
MATRIX_FUN_TY(MatrixTransposeMatrixMultiplyNaive)

#undef MATRIX_FUN_TY

// Initializes matrix entries.
static void initMatrix(Matrix& mat) {
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = i + j + 1;
    }
  }
}

template <int kRowA,
          int kColA,
          int kColB,
          DimType kDimType,
          template <int, int, int, int, int, DimType>
          class FunctorTy>
struct TestMatrixFunctions {
  void operator()() {
    Matrix A(kRowA, kColA);
    initMatrix(A);
    const int kRowB = kColA;
    Matrix B(kRowB, kColB);
    initMatrix(B);

    for (int row_stride_c = kRowA; row_stride_c < 3 * kRowA; ++row_stride_c) {
      for (int col_stride_c = kColB; col_stride_c < 3 * kColB; ++col_stride_c) {
        Matrix C(row_stride_c, col_stride_c);
        C.setOnes();

        Matrix C_plus = C;
        Matrix C_minus = C;
        Matrix C_assign = C;

        Matrix C_plus_ref = C;
        Matrix C_minus_ref = C;
        Matrix C_assign_ref = C;

        for (int start_row_c = 0; start_row_c + kRowA < row_stride_c;
             ++start_row_c) {
          for (int start_col_c = 0; start_col_c + kColB < col_stride_c;
               ++start_col_c) {
            C_plus_ref.block(start_row_c, start_col_c, kRowA, kColB) += A * B;
            FunctorTy<kRowA, kColA, kRowB, kColB, 1, kDimType>()(A.data(),
                                                                 kRowA,
                                                                 kColA,
                                                                 B.data(),
                                                                 kRowB,
                                                                 kColB,
                                                                 C_plus.data(),
                                                                 start_row_c,
                                                                 start_col_c,
                                                                 row_stride_c,
                                                                 col_stride_c);

            EXPECT_NEAR((C_plus_ref - C_plus).norm(), 0.0, kTolerance)
                << "C += A * B \n"
                << "row_stride_c : " << row_stride_c << "\n"
                << "col_stride_c : " << col_stride_c << "\n"
                << "start_row_c  : " << start_row_c << "\n"
                << "start_col_c  : " << start_col_c << "\n"
                << "Cref : \n"
                << C_plus_ref << "\n"
                << "C: \n"
                << C_plus;

            C_minus_ref.block(start_row_c, start_col_c, kRowA, kColB) -= A * B;
            FunctorTy<kRowA, kColA, kRowB, kColB, -1, kDimType>()(
                A.data(),
                kRowA,
                kColA,
                B.data(),
                kRowB,
                kColB,
                C_minus.data(),
                start_row_c,
                start_col_c,
                row_stride_c,
                col_stride_c);

            EXPECT_NEAR((C_minus_ref - C_minus).norm(), 0.0, kTolerance)
                << "C -= A * B \n"
                << "row_stride_c : " << row_stride_c << "\n"
                << "col_stride_c : " << col_stride_c << "\n"
                << "start_row_c  : " << start_row_c << "\n"
                << "start_col_c  : " << start_col_c << "\n"
                << "Cref : \n"
                << C_minus_ref << "\n"
                << "C: \n"
                << C_minus;

            C_assign_ref.block(start_row_c, start_col_c, kRowA, kColB) = A * B;

            FunctorTy<kRowA, kColA, kRowB, kColB, 0, kDimType>()(
                A.data(),
                kRowA,
                kColA,
                B.data(),
                kRowB,
                kColB,
                C_assign.data(),
                start_row_c,
                start_col_c,
                row_stride_c,
                col_stride_c);

            EXPECT_NEAR((C_assign_ref - C_assign).norm(), 0.0, kTolerance)
                << "C = A * B \n"
                << "row_stride_c : " << row_stride_c << "\n"
                << "col_stride_c : " << col_stride_c << "\n"
                << "start_row_c  : " << start_row_c << "\n"
                << "start_col_c  : " << start_col_c << "\n"
                << "Cref : \n"
                << C_assign_ref << "\n"
                << "C: \n"
                << C_assign;
          }
        }
      }
    }
  }
};

template <int kRowA,
          int kColA,
          int kColB,
          DimType kDimType,
          template <int, int, int, int, int, DimType>
          class FunctorTy>
struct TestMatrixTransposeFunctions {
  void operator()() {
    Matrix A(kRowA, kColA);
    initMatrix(A);
    const int kRowB = kRowA;
    Matrix B(kRowB, kColB);
    initMatrix(B);

    for (int row_stride_c = kColA; row_stride_c < 3 * kColA; ++row_stride_c) {
      for (int col_stride_c = kColB; col_stride_c < 3 * kColB; ++col_stride_c) {
        Matrix C(row_stride_c, col_stride_c);
        C.setOnes();

        Matrix C_plus = C;
        Matrix C_minus = C;
        Matrix C_assign = C;

        Matrix C_plus_ref = C;
        Matrix C_minus_ref = C;
        Matrix C_assign_ref = C;
        for (int start_row_c = 0; start_row_c + kColA < row_stride_c;
             ++start_row_c) {
          for (int start_col_c = 0; start_col_c + kColB < col_stride_c;
               ++start_col_c) {
            C_plus_ref.block(start_row_c, start_col_c, kColA, kColB) +=
                A.transpose() * B;

            FunctorTy<kRowA, kColA, kRowB, kColB, 1, kDimType>()(A.data(),
                                                                 kRowA,
                                                                 kColA,
                                                                 B.data(),
                                                                 kRowB,
                                                                 kColB,
                                                                 C_plus.data(),
                                                                 start_row_c,
                                                                 start_col_c,
                                                                 row_stride_c,
                                                                 col_stride_c);

            EXPECT_NEAR((C_plus_ref - C_plus).norm(), 0.0, kTolerance)
                << "C += A' * B \n"
                << "row_stride_c : " << row_stride_c << "\n"
                << "col_stride_c : " << col_stride_c << "\n"
                << "start_row_c  : " << start_row_c << "\n"
                << "start_col_c  : " << start_col_c << "\n"
                << "Cref : \n"
                << C_plus_ref << "\n"
                << "C: \n"
                << C_plus;

            C_minus_ref.block(start_row_c, start_col_c, kColA, kColB) -=
                A.transpose() * B;

            FunctorTy<kRowA, kColA, kRowB, kColB, -1, kDimType>()(
                A.data(),
                kRowA,
                kColA,
                B.data(),
                kRowB,
                kColB,
                C_minus.data(),
                start_row_c,
                start_col_c,
                row_stride_c,
                col_stride_c);

            EXPECT_NEAR((C_minus_ref - C_minus).norm(), 0.0, kTolerance)
                << "C -= A' * B \n"
                << "row_stride_c : " << row_stride_c << "\n"
                << "col_stride_c : " << col_stride_c << "\n"
                << "start_row_c  : " << start_row_c << "\n"
                << "start_col_c  : " << start_col_c << "\n"
                << "Cref : \n"
                << C_minus_ref << "\n"
                << "C: \n"
                << C_minus;

            C_assign_ref.block(start_row_c, start_col_c, kColA, kColB) =
                A.transpose() * B;

            FunctorTy<kRowA, kColA, kRowB, kColB, 0, kDimType>()(
                A.data(),
                kRowA,
                kColA,
                B.data(),
                kRowB,
                kColB,
                C_assign.data(),
                start_row_c,
                start_col_c,
                row_stride_c,
                col_stride_c);

            EXPECT_NEAR((C_assign_ref - C_assign).norm(), 0.0, kTolerance)
                << "C = A' * B \n"
                << "row_stride_c : " << row_stride_c << "\n"
                << "col_stride_c : " << col_stride_c << "\n"
                << "start_row_c  : " << start_row_c << "\n"
                << "start_col_c  : " << start_col_c << "\n"
                << "Cref : \n"
                << C_assign_ref << "\n"
                << "C: \n"
                << C_assign;
          }
        }
      }
    }
  }
};

TEST(BLAS, MatrixMatrixMultiply_5_3_7) {
  TestMatrixFunctions<5, 3, 7, DimType::Static, MatrixMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixMatrixMultiply_5_3_7_Dynamic) {
  TestMatrixFunctions<5, 3, 7, DimType::Dynamic, MatrixMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixMatrixMultiply_1_1_1) {
  TestMatrixFunctions<1, 1, 1, DimType::Static, MatrixMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixMatrixMultiply_1_1_1_Dynamic) {
  TestMatrixFunctions<1, 1, 1, DimType::Dynamic, MatrixMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixMatrixMultiply_9_9_9) {
  TestMatrixFunctions<9, 9, 9, DimType::Static, MatrixMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixMatrixMultiply_9_9_9_Dynamic) {
  TestMatrixFunctions<9, 9, 9, DimType::Dynamic, MatrixMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixMatrixMultiplyNaive_5_3_7) {
  TestMatrixFunctions<5,
                      3,
                      7,
                      DimType::Static,
                      MatrixMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixMatrixMultiplyNaive_5_3_7_Dynamic) {
  TestMatrixFunctions<5,
                      3,
                      7,
                      DimType::Dynamic,
                      MatrixMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixMatrixMultiplyNaive_1_1_1) {
  TestMatrixFunctions<1,
                      1,
                      1,
                      DimType::Static,
                      MatrixMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixMatrixMultiplyNaive_1_1_1_Dynamic) {
  TestMatrixFunctions<1,
                      1,
                      1,
                      DimType::Dynamic,
                      MatrixMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixMatrixMultiplyNaive_9_9_9) {
  TestMatrixFunctions<9,
                      9,
                      9,
                      DimType::Static,
                      MatrixMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixMatrixMultiplyNaive_9_9_9_Dynamic) {
  TestMatrixFunctions<9,
                      9,
                      9,
                      DimType::Dynamic,
                      MatrixMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiply_5_3_7) {
  TestMatrixTransposeFunctions<5,
                               3,
                               7,
                               DimType::Static,
                               MatrixTransposeMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiply_5_3_7_Dynamic) {
  TestMatrixTransposeFunctions<5,
                               3,
                               7,
                               DimType::Dynamic,
                               MatrixTransposeMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiply_1_1_1) {
  TestMatrixTransposeFunctions<1,
                               1,
                               1,
                               DimType::Static,
                               MatrixTransposeMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiply_1_1_1_Dynamic) {
  TestMatrixTransposeFunctions<1,
                               1,
                               1,
                               DimType::Dynamic,
                               MatrixTransposeMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiply_9_9_9) {
  TestMatrixTransposeFunctions<9,
                               9,
                               9,
                               DimType::Static,
                               MatrixTransposeMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiply_9_9_9_Dynamic) {
  TestMatrixTransposeFunctions<9,
                               9,
                               9,
                               DimType::Dynamic,
                               MatrixTransposeMatrixMultiplyTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiplyNaive_5_3_7) {
  TestMatrixTransposeFunctions<5,
                               3,
                               7,
                               DimType::Static,
                               MatrixTransposeMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiplyNaive_5_3_7_Dynamic) {
  TestMatrixTransposeFunctions<5,
                               3,
                               7,
                               DimType::Dynamic,
                               MatrixTransposeMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiplyNaive_1_1_1) {
  TestMatrixTransposeFunctions<1,
                               1,
                               1,
                               DimType::Static,
                               MatrixTransposeMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiplyNaive_1_1_1_Dynamic) {
  TestMatrixTransposeFunctions<1,
                               1,
                               1,
                               DimType::Dynamic,
                               MatrixTransposeMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiplyNaive_9_9_9) {
  TestMatrixTransposeFunctions<9,
                               9,
                               9,
                               DimType::Static,
                               MatrixTransposeMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixTransposeMatrixMultiplyNaive_9_9_9_Dynamic) {
  TestMatrixTransposeFunctions<9,
                               9,
                               9,
                               DimType::Dynamic,
                               MatrixTransposeMatrixMultiplyNaiveTy>()();
}

TEST(BLAS, MatrixVectorMultiply) {
  for (int num_rows_a = 1; num_rows_a < 10; ++num_rows_a) {
    for (int num_cols_a = 1; num_cols_a < 10; ++num_cols_a) {
      Matrix A(num_rows_a, num_cols_a);
      A.setOnes();

      Vector b(num_cols_a);
      b.setOnes();

      Vector c(num_rows_a);
      c.setOnes();

      Vector c_plus = c;
      Vector c_minus = c;
      Vector c_assign = c;

      Vector c_plus_ref = c;
      Vector c_minus_ref = c;
      Vector c_assign_ref = c;

      // clang-format off
      c_plus_ref += A * b;
      MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 1>(
          A.data(), num_rows_a, num_cols_a,
          b.data(),
          c_plus.data());
      EXPECT_NEAR((c_plus_ref - c_plus).norm(), 0.0, kTolerance)
          << "c += A * b \n"
          << "c_ref : \n" << c_plus_ref << "\n"
          << "c: \n" << c_plus;

      c_minus_ref -= A * b;
      MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, -1>(
          A.data(), num_rows_a, num_cols_a,
          b.data(),
          c_minus.data());
      EXPECT_NEAR((c_minus_ref - c_minus).norm(), 0.0, kTolerance)
          << "c -= A * b \n"
          << "c_ref : \n" << c_minus_ref << "\n"
          << "c: \n" << c_minus;

      c_assign_ref = A * b;
      MatrixVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 0>(
          A.data(), num_rows_a, num_cols_a,
          b.data(),
          c_assign.data());
      EXPECT_NEAR((c_assign_ref - c_assign).norm(), 0.0, kTolerance)
          << "c = A * b \n"
          << "c_ref : \n" << c_assign_ref << "\n"
          << "c: \n" << c_assign;
      // clang-format on
    }
  }
}

TEST(BLAS, MatrixTransposeVectorMultiply) {
  for (int num_rows_a = 1; num_rows_a < 10; ++num_rows_a) {
    for (int num_cols_a = 1; num_cols_a < 10; ++num_cols_a) {
      Matrix A(num_rows_a, num_cols_a);
      A.setRandom();

      Vector b(num_rows_a);
      b.setRandom();

      Vector c(num_cols_a);
      c.setOnes();

      Vector c_plus = c;
      Vector c_minus = c;
      Vector c_assign = c;

      Vector c_plus_ref = c;
      Vector c_minus_ref = c;
      Vector c_assign_ref = c;

      // clang-format off
      c_plus_ref += A.transpose() * b;
      MatrixTransposeVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 1>(
          A.data(), num_rows_a, num_cols_a,
          b.data(),
          c_plus.data());
      EXPECT_NEAR((c_plus_ref - c_plus).norm(), 0.0, kTolerance)
          << "c += A' * b \n"
          << "c_ref : \n" << c_plus_ref << "\n"
          << "c: \n" << c_plus;

      c_minus_ref -= A.transpose() * b;
      MatrixTransposeVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, -1>(
          A.data(), num_rows_a, num_cols_a,
          b.data(),
          c_minus.data());
      EXPECT_NEAR((c_minus_ref - c_minus).norm(), 0.0, kTolerance)
          << "c -= A' * b \n"
          << "c_ref : \n" << c_minus_ref << "\n"
          << "c: \n" << c_minus;

      c_assign_ref = A.transpose() * b;
      MatrixTransposeVectorMultiply<Eigen::Dynamic, Eigen::Dynamic, 0>(
          A.data(), num_rows_a, num_cols_a,
          b.data(),
          c_assign.data());
      EXPECT_NEAR((c_assign_ref - c_assign).norm(), 0.0, kTolerance)
          << "c = A' * b \n"
          << "c_ref : \n" << c_assign_ref << "\n"
          << "c: \n" << c_assign;
      // clang-format on
    }
  }
}

}  // namespace internal
}  // namespace ceres
