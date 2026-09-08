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

#include "ceres/cubic_interpolation.h"

#include <memory>

#include "ceres/jet.h"
#include "glog/logging.h"
#include "gtest/gtest.h"

namespace ceres::internal {

static constexpr double kTolerance = 1e-12;

TEST(Grid1D, OneDataDimension) {
  int x[] = {1, 2, 3};
  Grid1D<int, 1> grid(x, 0, 3);
  for (int i = 0; i < 3; ++i) {
    double value;
    grid.GetValue(i, &value);
    EXPECT_EQ(value, static_cast<double>(i + 1));
  }
}

TEST(Grid1D, OneDataDimensionOutOfBounds) {
  int x[] = {1, 2, 3};
  Grid1D<int, 1> grid(x, 0, 3);
  double value;
  grid.GetValue(-1, &value);
  EXPECT_EQ(value, x[0]);
  grid.GetValue(-2, &value);
  EXPECT_EQ(value, x[0]);
  grid.GetValue(3, &value);
  EXPECT_EQ(value, x[2]);
  grid.GetValue(4, &value);
  EXPECT_EQ(value, x[2]);
}

TEST(Grid1D, TwoDataDimensionIntegerDataInterleaved) {
  // clang-format off
  int x[] = {1, 5,
             2, 6,
             3, 7};
  // clang-format on

  Grid1D<int, 2, true> grid(x, 0, 3);
  for (int i = 0; i < 3; ++i) {
    double value[2];
    grid.GetValue(i, value);
    EXPECT_EQ(value[0], static_cast<double>(i + 1));
    EXPECT_EQ(value[1], static_cast<double>(i + 5));
  }
}

TEST(Grid1D, TwoDataDimensionIntegerDataStacked) {
  // clang-format off
  int x[] = {1, 2, 3,
             5, 6, 7};
  // clang-format on

  Grid1D<int, 2, false> grid(x, 0, 3);
  for (int i = 0; i < 3; ++i) {
    double value[2];
    grid.GetValue(i, value);
    EXPECT_EQ(value[0], static_cast<double>(i + 1));
    EXPECT_EQ(value[1], static_cast<double>(i + 5));
  }
}

TEST(Grid2D, OneDataDimensionRowMajor) {
  // clang-format off
  int x[] = {1, 2, 3,
             2, 3, 4};
  // clang-format on
  Grid2D<int, 1, true, true> grid(x, 0, 2, 0, 3);
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      double value;
      grid.GetValue(r, c, &value);
      EXPECT_EQ(value, static_cast<double>(r + c + 1));
    }
  }
}

TEST(Grid2D, OneDataDimensionRowMajorOutOfBounds) {
  // clang-format off
  int x[] = {1, 2, 3,
             2, 3, 4};
  // clang-format on
  Grid2D<int, 1, true, true> grid(x, 0, 2, 0, 3);
  double value;
  grid.GetValue(-1, -1, &value);
  EXPECT_EQ(value, x[0]);
  grid.GetValue(-1, 0, &value);
  EXPECT_EQ(value, x[0]);
  grid.GetValue(-1, 1, &value);
  EXPECT_EQ(value, x[1]);
  grid.GetValue(-1, 2, &value);
  EXPECT_EQ(value, x[2]);
  grid.GetValue(-1, 3, &value);
  EXPECT_EQ(value, x[2]);
  grid.GetValue(0, 3, &value);
  EXPECT_EQ(value, x[2]);
  grid.GetValue(1, 3, &value);
  EXPECT_EQ(value, x[5]);
  grid.GetValue(2, 3, &value);
  EXPECT_EQ(value, x[5]);
  grid.GetValue(2, 2, &value);
  EXPECT_EQ(value, x[5]);
  grid.GetValue(2, 1, &value);
  EXPECT_EQ(value, x[4]);
  grid.GetValue(2, 0, &value);
  EXPECT_EQ(value, x[3]);
  grid.GetValue(2, -1, &value);
  EXPECT_EQ(value, x[3]);
  grid.GetValue(1, -1, &value);
  EXPECT_EQ(value, x[3]);
  grid.GetValue(0, -1, &value);
  EXPECT_EQ(value, x[0]);
}

TEST(Grid2D, TwoDataDimensionRowMajorInterleaved) {
  // clang-format off
  int x[] = {1, 4, 2, 8, 3, 12,
             2, 8, 3, 12, 4, 16};
  // clang-format on
  Grid2D<int, 2, true, true> grid(x, 0, 2, 0, 3);
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      double value[2];
      grid.GetValue(r, c, value);
      EXPECT_EQ(value[0], static_cast<double>(r + c + 1));
      EXPECT_EQ(value[1], static_cast<double>(4 * (r + c + 1)));
    }
  }
}

TEST(Grid2D, TwoDataDimensionRowMajorStacked) {
  // clang-format off
  int x[] = {1,  2,  3,
             2,  3,  4,
             4,  8, 12,
             8, 12, 16};
  // clang-format on
  Grid2D<int, 2, true, false> grid(x, 0, 2, 0, 3);
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      double value[2];
      grid.GetValue(r, c, value);
      EXPECT_EQ(value[0], static_cast<double>(r + c + 1));
      EXPECT_EQ(value[1], static_cast<double>(4 * (r + c + 1)));
    }
  }
}

TEST(Grid2D, TwoDataDimensionColMajorInterleaved) {
  // clang-format off
  int x[] = { 1,  4, 2,  8,
              2,  8, 3, 12,
              3, 12, 4, 16};
  // clang-format on
  Grid2D<int, 2, false, true> grid(x, 0, 2, 0, 3);
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      double value[2];
      grid.GetValue(r, c, value);
      EXPECT_EQ(value[0], static_cast<double>(r + c + 1));
      EXPECT_EQ(value[1], static_cast<double>(4 * (r + c + 1)));
    }
  }
}

TEST(Grid2D, TwoDataDimensionColMajorStacked) {
  // clang-format off
  int x[] = {1,   2,
             2,   3,
             3,   4,
             4,   8,
             8,  12,
             12, 16};
  // clang-format on
  Grid2D<int, 2, false, false> grid(x, 0, 2, 0, 3);
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      double value[2];
      grid.GetValue(r, c, value);
      EXPECT_EQ(value[0], static_cast<double>(r + c + 1));
      EXPECT_EQ(value[1], static_cast<double>(4 * (r + c + 1)));
    }
  }
}

class CubicInterpolatorTest : public ::testing::Test {
 public:
  template <int kDataDimension>
  void RunPolynomialInterpolationTest(const double a,
                                      const double b,
                                      const double c,
                                      const double d) {
    values_ = std::make_unique<double[]>(kDataDimension * kNumSamples);

    for (int x = 0; x < kNumSamples; ++x) {
      for (int dim = 0; dim < kDataDimension; ++dim) {
        values_[x * kDataDimension + dim] =
            (dim * dim + 1) * (a * x * x * x + b * x * x + c * x + d);
      }
    }

    Grid1D<double, kDataDimension> grid(values_.get(), 0, kNumSamples);
    CubicInterpolator<Grid1D<double, kDataDimension>> interpolator(grid);

    // Check values in the all the cells but the first and the last
    // ones. In these cells, the interpolated function values should
    // match exactly the values of the function being interpolated.
    //
    // On the boundary, we extrapolate the values of the function on
    // the basis of its first derivative, so we do not expect the
    // function values and its derivatives not to match.
    for (int j = 0; j < kNumTestSamples; ++j) {
      const double x = 1.0 + 7.0 / (kNumTestSamples - 1) * j;
      double expected_f[kDataDimension], expected_dfdx[kDataDimension];
      double f[kDataDimension], dfdx[kDataDimension];

      for (int dim = 0; dim < kDataDimension; ++dim) {
        expected_f[dim] =
            (dim * dim + 1) * (a * x * x * x + b * x * x + c * x + d);
        expected_dfdx[dim] =
            (dim * dim + 1) * (3.0 * a * x * x + 2.0 * b * x + c);
      }

      interpolator.Evaluate(x, f, dfdx);
      for (int dim = 0; dim < kDataDimension; ++dim) {
        EXPECT_NEAR(f[dim], expected_f[dim], kTolerance)
            << "x: " << x << " dim: " << dim
            << " actual f(x): " << expected_f[dim]
            << " estimated f(x): " << f[dim];
        EXPECT_NEAR(dfdx[dim], expected_dfdx[dim], kTolerance)
            << "x: " << x << " dim: " << dim
            << " actual df(x)/dx: " << expected_dfdx[dim]
            << " estimated df(x)/dx: " << dfdx[dim];
      }
    }
  }

 private:
  static constexpr int kNumSamples = 10;
  static constexpr int kNumTestSamples = 100;
  std::unique_ptr<double[]> values_;
};

TEST_F(CubicInterpolatorTest, ConstantFunction) {
  RunPolynomialInterpolationTest<1>(0.0, 0.0, 0.0, 0.5);
  RunPolynomialInterpolationTest<2>(0.0, 0.0, 0.0, 0.5);
  RunPolynomialInterpolationTest<3>(0.0, 0.0, 0.0, 0.5);
}

TEST_F(CubicInterpolatorTest, LinearFunction) {
  RunPolynomialInterpolationTest<1>(0.0, 0.0, 1.0, 0.5);
  RunPolynomialInterpolationTest<2>(0.0, 0.0, 1.0, 0.5);
  RunPolynomialInterpolationTest<3>(0.0, 0.0, 1.0, 0.5);
}

TEST_F(CubicInterpolatorTest, QuadraticFunction) {
  RunPolynomialInterpolationTest<1>(0.0, 0.4, 1.0, 0.5);
  RunPolynomialInterpolationTest<2>(0.0, 0.4, 1.0, 0.5);
  RunPolynomialInterpolationTest<3>(0.0, 0.4, 1.0, 0.5);
}

TEST(CubicInterpolator, JetEvaluation) {
  const double values[] = {1.0, 2.0, 2.0, 5.0, 3.0, 9.0, 2.0, 7.0};

  Grid1D<double, 2, true> grid(values, 0, 4);
  CubicInterpolator<Grid1D<double, 2, true>> interpolator(grid);

  double f[2], dfdx[2];
  const double x = 2.5;
  interpolator.Evaluate(x, f, dfdx);

  // Create a Jet with the same scalar part as x, so that the output
  // Jet will be evaluated at x.
  Jet<double, 4> x_jet;
  x_jet.a = x;
  x_jet.v(0) = 1.0;
  x_jet.v(1) = 1.1;
  x_jet.v(2) = 1.2;
  x_jet.v(3) = 1.3;

  Jet<double, 4> f_jets[2];
  interpolator.Evaluate(x_jet, f_jets);

  // Check that the scalar part of the Jet is f(x).
  EXPECT_EQ(f_jets[0].a, f[0]);
  EXPECT_EQ(f_jets[1].a, f[1]);

  // Check that the derivative part of the Jet is dfdx * x_jet.v
  // by the chain rule.
  EXPECT_NEAR((f_jets[0].v - dfdx[0] * x_jet.v).norm(), 0.0, kTolerance);
  EXPECT_NEAR((f_jets[1].v - dfdx[1] * x_jet.v).norm(), 0.0, kTolerance);
}

class BiCubicInterpolatorTest : public ::testing::Test {
 public:
  // This class needs to have an Eigen aligned operator new as it contains
  // fixed-size Eigen types.
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  template <int kDataDimension>
  void RunPolynomialInterpolationTest(const Eigen::Matrix3d& coeff) {
    values_ = std::make_unique<double[]>(kNumRows * kNumCols * kDataDimension);
    coeff_ = coeff;
    double* v = values_.get();
    for (int r = 0; r < kNumRows; ++r) {
      for (int c = 0; c < kNumCols; ++c) {
        for (int dim = 0; dim < kDataDimension; ++dim) {
          *v++ = (dim * dim + 1) * EvaluateF(r, c);
        }
      }
    }

    Grid2D<double, kDataDimension> grid(
        values_.get(), 0, kNumRows, 0, kNumCols);
    BiCubicInterpolator<Grid2D<double, kDataDimension>> interpolator(grid);

    for (int j = 0; j < kNumRowSamples; ++j) {
      const double r = 1.0 + 7.0 / (kNumRowSamples - 1) * j;
      for (int k = 0; k < kNumColSamples; ++k) {
        const double c = 1.0 + 7.0 / (kNumColSamples - 1) * k;
        double f[kDataDimension], dfdr[kDataDimension], dfdc[kDataDimension];
        interpolator.Evaluate(r, c, f, dfdr, dfdc);
        for (int dim = 0; dim < kDataDimension; ++dim) {
          EXPECT_NEAR(f[dim], (dim * dim + 1) * EvaluateF(r, c), kTolerance);
          EXPECT_NEAR(
              dfdr[dim], (dim * dim + 1) * EvaluatedFdr(r, c), kTolerance);
          EXPECT_NEAR(
              dfdc[dim], (dim * dim + 1) * EvaluatedFdc(r, c), kTolerance);
        }
      }
    }
  }

 private:
  double EvaluateF(double r, double c) {
    Eigen::Vector3d x;
    x(0) = r;
    x(1) = c;
    x(2) = 1;
    return x.transpose() * coeff_ * x;
  }

  double EvaluatedFdr(double r, double c) {
    Eigen::Vector3d x;
    x(0) = r;
    x(1) = c;
    x(2) = 1;
    return (coeff_.row(0) + coeff_.col(0).transpose()) * x;
  }

  double EvaluatedFdc(double r, double c) {
    Eigen::Vector3d x;
    x(0) = r;
    x(1) = c;
    x(2) = 1;
    return (coeff_.row(1) + coeff_.col(1).transpose()) * x;
  }

  Eigen::Matrix3d coeff_;
  static constexpr int kNumRows = 10;
  static constexpr int kNumCols = 10;
  static constexpr int kNumRowSamples = 100;
  static constexpr int kNumColSamples = 100;
  std::unique_ptr<double[]> values_;
};

TEST_F(BiCubicInterpolatorTest, ZeroFunction) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree00Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree01Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  coeff(0, 2) = 0.1;
  coeff(2, 0) = 0.1;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree10Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  coeff(0, 1) = 0.1;
  coeff(1, 0) = 0.1;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree11Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  coeff(0, 1) = 0.1;
  coeff(1, 0) = 0.1;
  coeff(0, 2) = 0.2;
  coeff(2, 0) = 0.2;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree12Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  coeff(0, 1) = 0.1;
  coeff(1, 0) = 0.1;
  coeff(0, 2) = 0.2;
  coeff(2, 0) = 0.2;
  coeff(1, 1) = 0.3;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree21Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  coeff(0, 1) = 0.1;
  coeff(1, 0) = 0.1;
  coeff(0, 2) = 0.2;
  coeff(2, 0) = 0.2;
  coeff(0, 0) = 0.3;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST_F(BiCubicInterpolatorTest, Degree22Function) {
  Eigen::Matrix3d coeff = Eigen::Matrix3d::Zero();
  coeff(2, 2) = 1.0;
  coeff(0, 1) = 0.1;
  coeff(1, 0) = 0.1;
  coeff(0, 2) = 0.2;
  coeff(2, 0) = 0.2;
  coeff(0, 0) = 0.3;
  coeff(0, 1) = -0.4;
  coeff(1, 0) = -0.4;
  RunPolynomialInterpolationTest<1>(coeff);
  RunPolynomialInterpolationTest<2>(coeff);
  RunPolynomialInterpolationTest<3>(coeff);
}

TEST(BiCubicInterpolator, JetEvaluation) {
  // clang-format off
  const double values[] = {1.0, 5.0, 2.0, 10.0, 2.0, 6.0, 3.0, 5.0,
                           1.0, 2.0, 2.0,  2.0, 2.0, 2.0, 3.0, 1.0};
  // clang-format on

  Grid2D<double, 2> grid(values, 0, 2, 0, 4);
  BiCubicInterpolator<Grid2D<double, 2>> interpolator(grid);

  double f[2], dfdr[2], dfdc[2];
  const double r = 0.5;
  const double c = 2.5;
  interpolator.Evaluate(r, c, f, dfdr, dfdc);

  // Create a Jet with the same scalar part as x, so that the output
  // Jet will be evaluated at x.
  Jet<double, 4> r_jet;
  r_jet.a = r;
  r_jet.v(0) = 1.0;
  r_jet.v(1) = 1.1;
  r_jet.v(2) = 1.2;
  r_jet.v(3) = 1.3;

  Jet<double, 4> c_jet;
  c_jet.a = c;
  c_jet.v(0) = 2.0;
  c_jet.v(1) = 3.1;
  c_jet.v(2) = 4.2;
  c_jet.v(3) = 5.3;

  Jet<double, 4> f_jets[2];
  interpolator.Evaluate(r_jet, c_jet, f_jets);
  EXPECT_EQ(f_jets[0].a, f[0]);
  EXPECT_EQ(f_jets[1].a, f[1]);
  EXPECT_NEAR((f_jets[0].v - dfdr[0] * r_jet.v - dfdc[0] * c_jet.v).norm(),
              0.0,
              kTolerance);
  EXPECT_NEAR((f_jets[1].v - dfdr[1] * r_jet.v - dfdc[1] * c_jet.v).norm(),
              0.0,
              kTolerance);
}

}  // namespace ceres::internal
