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
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include <cmath>
#include <limits>
#include "Eigen/Geometry"
#include "ceres/autodiff_local_parameterization.h"
#include "ceres/fpclassify.h"
#include "ceres/householder_vector.h"
#include "ceres/internal/autodiff.h"
#include "ceres/internal/eigen.h"
#include "ceres/local_parameterization.h"
#include "ceres/random.h"
#include "ceres/rotation.h"
#include "gtest/gtest.h"

namespace ceres {
namespace internal {

TEST(IdentityParameterization, EverythingTest) {
  IdentityParameterization parameterization(3);
  EXPECT_EQ(parameterization.GlobalSize(), 3);
  EXPECT_EQ(parameterization.LocalSize(), 3);

  double x[3] = {1.0, 2.0, 3.0};
  double delta[3] = {0.0, 1.0, 2.0};
  double x_plus_delta[3] = {0.0, 0.0, 0.0};
  parameterization.Plus(x, delta, x_plus_delta);
  EXPECT_EQ(x_plus_delta[0], 1.0);
  EXPECT_EQ(x_plus_delta[1], 3.0);
  EXPECT_EQ(x_plus_delta[2], 5.0);

  double jacobian[9];
  parameterization.ComputeJacobian(x, jacobian);
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j, ++k) {
      EXPECT_EQ(jacobian[k], (i == j) ? 1.0 : 0.0);
    }
  }

  Matrix global_matrix = Matrix::Ones(10, 3);
  Matrix local_matrix = Matrix::Zero(10, 3);
  parameterization.MultiplyByJacobian(x,
                                      10,
                                      global_matrix.data(),
                                      local_matrix.data());
  EXPECT_EQ((local_matrix - global_matrix).norm(), 0.0);
}


TEST(SubsetParameterization, NegativeParameterIndexDeathTest) {
  std::vector<int> constant_parameters;
  constant_parameters.push_back(-1);
  EXPECT_DEATH_IF_SUPPORTED(
      SubsetParameterization parameterization(2, constant_parameters),
      "greater than zero");
}

TEST(SubsetParameterization, GreaterThanSizeParameterIndexDeathTest) {
  std::vector<int> constant_parameters;
  constant_parameters.push_back(2);
  EXPECT_DEATH_IF_SUPPORTED(
      SubsetParameterization parameterization(2, constant_parameters),
      "less than the size");
}

TEST(SubsetParameterization, DuplicateParametersDeathTest) {
  std::vector<int> constant_parameters;
  constant_parameters.push_back(1);
  constant_parameters.push_back(1);
  EXPECT_DEATH_IF_SUPPORTED(
      SubsetParameterization parameterization(2, constant_parameters),
      "duplicates");
}

TEST(SubsetParameterization,
     ProductParameterizationWithZeroLocalSizeSubsetParameterization1) {
  std::vector<int> constant_parameters;
  constant_parameters.push_back(0);
  LocalParameterization* subset_param =
      new SubsetParameterization(1, constant_parameters);
  LocalParameterization* identity_param = new IdentityParameterization(2);
  ProductParameterization product_param(subset_param, identity_param);
  EXPECT_EQ(product_param.GlobalSize(), 3);
  EXPECT_EQ(product_param.LocalSize(), 2);
  double x[] = {1.0, 1.0, 1.0};
  double delta[] = {2.0, 3.0};
  double x_plus_delta[] = {0.0, 0.0, 0.0};
  EXPECT_TRUE(product_param.Plus(x, delta, x_plus_delta));
  EXPECT_EQ(x_plus_delta[0], x[0]);
  EXPECT_EQ(x_plus_delta[1], x[1] + delta[0]);
  EXPECT_EQ(x_plus_delta[2], x[2] + delta[1]);

  Matrix actual_jacobian(3, 2);
  EXPECT_TRUE(product_param.ComputeJacobian(x, actual_jacobian.data()));
}

TEST(SubsetParameterization,
     ProductParameterizationWithZeroLocalSizeSubsetParameterization2) {
  std::vector<int> constant_parameters;
  constant_parameters.push_back(0);
  LocalParameterization* subset_param =
      new SubsetParameterization(1, constant_parameters);
  LocalParameterization* identity_param = new IdentityParameterization(2);
  ProductParameterization product_param(identity_param, subset_param);
  EXPECT_EQ(product_param.GlobalSize(), 3);
  EXPECT_EQ(product_param.LocalSize(), 2);
  double x[] = {1.0, 1.0, 1.0};
  double delta[] = {2.0, 3.0};
  double x_plus_delta[] = {0.0, 0.0, 0.0};
  EXPECT_TRUE(product_param.Plus(x, delta, x_plus_delta));
  EXPECT_EQ(x_plus_delta[0], x[0] + delta[0]);
  EXPECT_EQ(x_plus_delta[1], x[1] + delta[1]);
  EXPECT_EQ(x_plus_delta[2], x[2]);

  Matrix actual_jacobian(3, 2);
  EXPECT_TRUE(product_param.ComputeJacobian(x, actual_jacobian.data()));
}

TEST(SubsetParameterization, NormalFunctionTest) {
  const int kGlobalSize = 4;
  const int kLocalSize = 3;

  double x[kGlobalSize] = {1.0, 2.0, 3.0, 4.0};
  for (int i = 0; i < kGlobalSize; ++i) {
    std::vector<int> constant_parameters;
    constant_parameters.push_back(i);
    SubsetParameterization parameterization(kGlobalSize, constant_parameters);
    double delta[kLocalSize] = {1.0, 2.0, 3.0};
    double x_plus_delta[kGlobalSize] = {0.0, 0.0, 0.0};

    parameterization.Plus(x, delta, x_plus_delta);
    int k = 0;
    for (int j = 0; j < kGlobalSize; ++j) {
      if (j == i)  {
        EXPECT_EQ(x_plus_delta[j], x[j]);
      } else {
        EXPECT_EQ(x_plus_delta[j], x[j] + delta[k++]);
      }
    }

    double jacobian[kGlobalSize * kLocalSize];
    parameterization.ComputeJacobian(x, jacobian);
    int delta_cursor = 0;
    int jacobian_cursor = 0;
    for (int j = 0; j < kGlobalSize; ++j) {
      if (j != i) {
        for (int k = 0; k < kLocalSize; ++k, jacobian_cursor++) {
          EXPECT_EQ(jacobian[jacobian_cursor], delta_cursor == k ? 1.0 : 0.0);
        }
        ++delta_cursor;
      } else {
        for (int k = 0; k < kLocalSize; ++k, jacobian_cursor++) {
          EXPECT_EQ(jacobian[jacobian_cursor], 0.0);
        }
      }
    }

    Matrix global_matrix = Matrix::Ones(10, kGlobalSize);
    for (int row = 0; row < kGlobalSize; ++row) {
      for (int col = 0; col < kGlobalSize; ++col) {
        global_matrix(row, col) = col;
      }
    }

    Matrix local_matrix = Matrix::Zero(10, kLocalSize);
    parameterization.MultiplyByJacobian(x,
                                        10,
                                        global_matrix.data(),
                                        local_matrix.data());
    Matrix expected_local_matrix =
        global_matrix * MatrixRef(jacobian, kGlobalSize, kLocalSize);
    EXPECT_EQ((local_matrix - expected_local_matrix).norm(), 0.0);
  }
}

// Functor needed to implement automatically differentiated Plus for
// quaternions.
struct QuaternionPlus {
  template<typename T>
  bool operator()(const T* x, const T* delta, T* x_plus_delta) const {
    const T squared_norm_delta =
        delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

    T q_delta[4];
    if (squared_norm_delta > T(0.0)) {
      T norm_delta = sqrt(squared_norm_delta);
      const T sin_delta_by_delta = sin(norm_delta) / norm_delta;
      q_delta[0] = cos(norm_delta);
      q_delta[1] = sin_delta_by_delta * delta[0];
      q_delta[2] = sin_delta_by_delta * delta[1];
      q_delta[3] = sin_delta_by_delta * delta[2];
    } else {
      // We do not just use q_delta = [1,0,0,0] here because that is a
      // constant and when used for automatic differentiation will
      // lead to a zero derivative. Instead we take a first order
      // approximation and evaluate it at zero.
      q_delta[0] = T(1.0);
      q_delta[1] = delta[0];
      q_delta[2] = delta[1];
      q_delta[3] = delta[2];
    }

    QuaternionProduct(q_delta, x, x_plus_delta);
    return true;
  }
};

template<typename Parameterization, typename Plus>
void QuaternionParameterizationTestHelper(
    const double* x, const double* delta,
    const double* x_plus_delta_ref) {
  const int kGlobalSize = 4;
  const int kLocalSize = 3;

  const double kTolerance = 1e-14;

  double x_plus_delta[kGlobalSize] = {0.0, 0.0, 0.0, 0.0};
  Parameterization parameterization;
  parameterization.Plus(x, delta, x_plus_delta);
  for (int i = 0; i < kGlobalSize; ++i) {
    EXPECT_NEAR(x_plus_delta[i], x_plus_delta[i], kTolerance);
  }

  const double x_plus_delta_norm =
      sqrt(x_plus_delta[0] * x_plus_delta[0] +
           x_plus_delta[1] * x_plus_delta[1] +
           x_plus_delta[2] * x_plus_delta[2] +
           x_plus_delta[3] * x_plus_delta[3]);

  EXPECT_NEAR(x_plus_delta_norm, 1.0, kTolerance);

  double jacobian_ref[12];
  double zero_delta[kLocalSize] = {0.0, 0.0, 0.0};
  const double* parameters[2] = {x, zero_delta};
  double* jacobian_array[2] = { NULL, jacobian_ref };

  // Autodiff jacobian at delta_x = 0.
  internal::AutoDiff<Plus,
                     double,
                     kGlobalSize,
                     kLocalSize>::Differentiate(Plus(),
                                                parameters,
                                                kGlobalSize,
                                                x_plus_delta,
                                                jacobian_array);

  double jacobian[12];
  parameterization.ComputeJacobian(x, jacobian);
  for (int i = 0; i < 12; ++i) {
    EXPECT_TRUE(IsFinite(jacobian[i]));
    EXPECT_NEAR(jacobian[i], jacobian_ref[i], kTolerance)
        << "Jacobian mismatch: i = " << i
        << "\n Expected \n"
        << ConstMatrixRef(jacobian_ref, kGlobalSize, kLocalSize)
        << "\n Actual \n"
        << ConstMatrixRef(jacobian, kGlobalSize, kLocalSize);
  }

  Matrix global_matrix = Matrix::Random(10, kGlobalSize);
  Matrix local_matrix = Matrix::Zero(10, kLocalSize);
  parameterization.MultiplyByJacobian(x,
                                      10,
                                      global_matrix.data(),
                                      local_matrix.data());
  Matrix expected_local_matrix =
      global_matrix * MatrixRef(jacobian, kGlobalSize, kLocalSize);
  EXPECT_NEAR((local_matrix - expected_local_matrix).norm(),
              0.0,
              10.0 * std::numeric_limits<double>::epsilon());
}

template <int N>
void Normalize(double* x) {
  VectorRef(x, N).normalize();
}

TEST(QuaternionParameterization, ZeroTest) {
  double x[4] = {0.5, 0.5, 0.5, 0.5};
  double delta[3] = {0.0, 0.0, 0.0};
  double q_delta[4] = {1.0, 0.0, 0.0, 0.0};
  double x_plus_delta[4] = {0.0, 0.0, 0.0, 0.0};
  QuaternionProduct(q_delta, x, x_plus_delta);
  QuaternionParameterizationTestHelper<QuaternionParameterization,
                                       QuaternionPlus>(x, delta, x_plus_delta);
}

TEST(QuaternionParameterization, NearZeroTest) {
  double x[4] = {0.52, 0.25, 0.15, 0.45};
  Normalize<4>(x);

  double delta[3] = {0.24, 0.15, 0.10};
  for (int i = 0; i < 3; ++i) {
    delta[i] = delta[i] * 1e-14;
  }

  double q_delta[4];
  q_delta[0] = 1.0;
  q_delta[1] = delta[0];
  q_delta[2] = delta[1];
  q_delta[3] = delta[2];

  double x_plus_delta[4] = {0.0, 0.0, 0.0, 0.0};
  QuaternionProduct(q_delta, x, x_plus_delta);
  QuaternionParameterizationTestHelper<QuaternionParameterization,
                                       QuaternionPlus>(x, delta, x_plus_delta);
}

TEST(QuaternionParameterization, AwayFromZeroTest) {
  double x[4] = {0.52, 0.25, 0.15, 0.45};
  Normalize<4>(x);

  double delta[3] = {0.24, 0.15, 0.10};
  const double delta_norm = sqrt(delta[0] * delta[0] +
                                 delta[1] * delta[1] +
                                 delta[2] * delta[2]);
  double q_delta[4];
  q_delta[0] = cos(delta_norm);
  q_delta[1] = sin(delta_norm) / delta_norm * delta[0];
  q_delta[2] = sin(delta_norm) / delta_norm * delta[1];
  q_delta[3] = sin(delta_norm) / delta_norm * delta[2];

  double x_plus_delta[4] = {0.0, 0.0, 0.0, 0.0};
  QuaternionProduct(q_delta, x, x_plus_delta);
  QuaternionParameterizationTestHelper<QuaternionParameterization,
                                       QuaternionPlus>(x, delta, x_plus_delta);
}

// Functor needed to implement automatically differentiated Plus for
// Eigen's quaternion.
struct EigenQuaternionPlus {
  template<typename T>
  bool operator()(const T* x, const T* delta, T* x_plus_delta) const {
    const T norm_delta =
        sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);

    Eigen::Quaternion<T> q_delta;
    if (norm_delta > T(0.0)) {
      const T sin_delta_by_delta = sin(norm_delta) / norm_delta;
      q_delta.coeffs() << sin_delta_by_delta * delta[0],
          sin_delta_by_delta * delta[1], sin_delta_by_delta * delta[2],
          cos(norm_delta);
    } else {
      // We do not just use q_delta = [0,0,0,1] here because that is a
      // constant and when used for automatic differentiation will
      // lead to a zero derivative. Instead we take a first order
      // approximation and evaluate it at zero.
      q_delta.coeffs() <<  delta[0], delta[1], delta[2], T(1.0);
    }

    Eigen::Map<Eigen::Quaternion<T> > x_plus_delta_ref(x_plus_delta);
    Eigen::Map<const Eigen::Quaternion<T> > x_ref(x);
    x_plus_delta_ref = q_delta * x_ref;
    return true;
  }
};

TEST(EigenQuaternionParameterization, ZeroTest) {
  Eigen::Quaterniond x(0.5, 0.5, 0.5, 0.5);
  double delta[3] = {0.0, 0.0, 0.0};
  Eigen::Quaterniond q_delta(1.0, 0.0, 0.0, 0.0);
  Eigen::Quaterniond x_plus_delta = q_delta * x;
  QuaternionParameterizationTestHelper<EigenQuaternionParameterization,
                                       EigenQuaternionPlus>(
      x.coeffs().data(), delta, x_plus_delta.coeffs().data());
}

TEST(EigenQuaternionParameterization, NearZeroTest) {
  Eigen::Quaterniond x(0.52, 0.25, 0.15, 0.45);
  x.normalize();

  double delta[3] = {0.24, 0.15, 0.10};
  for (int i = 0; i < 3; ++i) {
    delta[i] = delta[i] * 1e-14;
  }

  // Note: w is first in the constructor.
  Eigen::Quaterniond q_delta(1.0, delta[0], delta[1], delta[2]);

  Eigen::Quaterniond x_plus_delta = q_delta * x;
  QuaternionParameterizationTestHelper<EigenQuaternionParameterization,
                                       EigenQuaternionPlus>(
      x.coeffs().data(), delta, x_plus_delta.coeffs().data());
}

TEST(EigenQuaternionParameterization, AwayFromZeroTest) {
  Eigen::Quaterniond x(0.52, 0.25, 0.15, 0.45);
  x.normalize();

  double delta[3] = {0.24, 0.15, 0.10};
  const double delta_norm = sqrt(delta[0] * delta[0] +
                                 delta[1] * delta[1] +
                                 delta[2] * delta[2]);

  // Note: w is first in the constructor.
  Eigen::Quaterniond q_delta(cos(delta_norm),
                             sin(delta_norm) / delta_norm * delta[0],
                             sin(delta_norm) / delta_norm * delta[1],
                             sin(delta_norm) / delta_norm * delta[2]);

  Eigen::Quaterniond x_plus_delta = q_delta * x;
  QuaternionParameterizationTestHelper<EigenQuaternionParameterization,
                                       EigenQuaternionPlus>(
      x.coeffs().data(), delta, x_plus_delta.coeffs().data());
}

// Functor needed to implement automatically differentiated Plus for
// homogeneous vectors. Note this explicitly defined for vectors of size 4.
struct HomogeneousVectorParameterizationPlus {
  template<typename Scalar>
  bool operator()(const Scalar* p_x, const Scalar* p_delta,
                  Scalar* p_x_plus_delta) const {
    Eigen::Map<const Eigen::Matrix<Scalar, 4, 1> > x(p_x);
    Eigen::Map<const Eigen::Matrix<Scalar, 3, 1> > delta(p_delta);
    Eigen::Map<Eigen::Matrix<Scalar, 4, 1> > x_plus_delta(p_x_plus_delta);

    const Scalar squared_norm_delta =
        delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

    Eigen::Matrix<Scalar, 4, 1> y;
    Scalar one_half(0.5);
    if (squared_norm_delta > Scalar(0.0)) {
      Scalar norm_delta = sqrt(squared_norm_delta);
      Scalar norm_delta_div_2 = 0.5 * norm_delta;
      const Scalar sin_delta_by_delta = sin(norm_delta_div_2) /
          norm_delta_div_2;
      y[0] = sin_delta_by_delta * delta[0] * one_half;
      y[1] = sin_delta_by_delta * delta[1] * one_half;
      y[2] = sin_delta_by_delta * delta[2] * one_half;
      y[3] = cos(norm_delta_div_2);

    } else {
      // We do not just use y = [0,0,0,1] here because that is a
      // constant and when used for automatic differentiation will
      // lead to a zero derivative. Instead we take a first order
      // approximation and evaluate it at zero.
      y[0] = delta[0] * one_half;
      y[1] = delta[1] * one_half;
      y[2] = delta[2] * one_half;
      y[3] = Scalar(1.0);
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(4);
    Scalar beta;
    internal::ComputeHouseholderVector<Scalar>(x, &v, &beta);

    x_plus_delta = x.norm() * (y - v * (beta * v.dot(y)));

    return true;
  }
};

void HomogeneousVectorParameterizationHelper(const double* x,
                                             const double* delta) {
  const double kTolerance = 1e-14;

  HomogeneousVectorParameterization homogeneous_vector_parameterization(4);

  // Ensure the update maintains the norm.
  double x_plus_delta[4] = {0.0, 0.0, 0.0, 0.0};
  homogeneous_vector_parameterization.Plus(x, delta, x_plus_delta);

  const double x_plus_delta_norm =
      sqrt(x_plus_delta[0] * x_plus_delta[0] +
           x_plus_delta[1] * x_plus_delta[1] +
           x_plus_delta[2] * x_plus_delta[2] +
           x_plus_delta[3] * x_plus_delta[3]);

  const double x_norm = sqrt(x[0] * x[0] + x[1] * x[1] +
                             x[2] * x[2] + x[3] * x[3]);

  EXPECT_NEAR(x_plus_delta_norm, x_norm, kTolerance);

  // Autodiff jacobian at delta_x = 0.
  AutoDiffLocalParameterization<HomogeneousVectorParameterizationPlus, 4, 3>
      autodiff_jacobian;

  double jacobian_autodiff[12];
  double jacobian_analytic[12];

  homogeneous_vector_parameterization.ComputeJacobian(x, jacobian_analytic);
  autodiff_jacobian.ComputeJacobian(x, jacobian_autodiff);

  for (int i = 0; i < 12; ++i) {
    EXPECT_TRUE(ceres::IsFinite(jacobian_analytic[i]));
    EXPECT_NEAR(jacobian_analytic[i], jacobian_autodiff[i], kTolerance)
        << "Jacobian mismatch: i = " << i << ", " << jacobian_analytic[i] << " "
        << jacobian_autodiff[i];
  }
}

TEST(HomogeneousVectorParameterization, ZeroTest) {
  double x[4] = {0.0, 0.0, 0.0, 1.0};
  Normalize<4>(x);
  double delta[3] = {0.0, 0.0, 0.0};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, NearZeroTest1) {
  double x[4] = {1e-5, 1e-5, 1e-5, 1.0};
  Normalize<4>(x);
  double delta[3] = {0.0, 1.0, 0.0};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, NearZeroTest2) {
  double x[4] = {0.001, 0.0, 0.0, 0.0};
  double delta[3] = {0.0, 1.0, 0.0};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, AwayFromZeroTest1) {
  double x[4] = {0.52, 0.25, 0.15, 0.45};
  Normalize<4>(x);
  double delta[3] = {0.0, 1.0, -0.5};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, AwayFromZeroTest2) {
  double x[4] = {0.87, -0.25, -0.34, 0.45};
  Normalize<4>(x);
  double delta[3] = {0.0, 0.0, -0.5};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, AwayFromZeroTest3) {
  double x[4] = {0.0, 0.0, 0.0, 2.0};
  double delta[3] = {0.0, 0.0, 0};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, AwayFromZeroTest4) {
  double x[4] = {0.2, -1.0, 0.0, 2.0};
  double delta[3] = {1.4, 0.0, -0.5};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, AwayFromZeroTest5) {
  double x[4] = {2.0, 0.0, 0.0, 0.0};
  double delta[3] = {1.4, 0.0, -0.5};

  HomogeneousVectorParameterizationHelper(x, delta);
}

TEST(HomogeneousVectorParameterization, DeathTests) {
  EXPECT_DEATH_IF_SUPPORTED(HomogeneousVectorParameterization x(1), "size");
}


class ProductParameterizationTest : public ::testing::Test {
 protected :
  virtual void SetUp() {
    const int global_size1 = 5;
    std::vector<int> constant_parameters1;
    constant_parameters1.push_back(2);
    param1_.reset(new SubsetParameterization(global_size1,
                                             constant_parameters1));

    const int global_size2 = 3;
    std::vector<int> constant_parameters2;
    constant_parameters2.push_back(0);
    constant_parameters2.push_back(1);
    param2_.reset(new SubsetParameterization(global_size2,
                                             constant_parameters2));

    const int global_size3 = 4;
    std::vector<int> constant_parameters3;
    constant_parameters3.push_back(1);
    param3_.reset(new SubsetParameterization(global_size3,
                                             constant_parameters3));

    const int global_size4 = 2;
    std::vector<int> constant_parameters4;
    constant_parameters4.push_back(1);
    param4_.reset(new SubsetParameterization(global_size4,
                                             constant_parameters4));
  }

  scoped_ptr<LocalParameterization> param1_;
  scoped_ptr<LocalParameterization> param2_;
  scoped_ptr<LocalParameterization> param3_;
  scoped_ptr<LocalParameterization> param4_;
};

TEST_F(ProductParameterizationTest, LocalAndGlobalSize2) {
  LocalParameterization* param1 = param1_.release();
  LocalParameterization* param2 = param2_.release();

  ProductParameterization product_param(param1, param2);
  EXPECT_EQ(product_param.LocalSize(),
            param1->LocalSize() + param2->LocalSize());
  EXPECT_EQ(product_param.GlobalSize(),
            param1->GlobalSize() + param2->GlobalSize());
}


TEST_F(ProductParameterizationTest, LocalAndGlobalSize3) {
  LocalParameterization* param1 = param1_.release();
  LocalParameterization* param2 = param2_.release();
  LocalParameterization* param3 = param3_.release();

  ProductParameterization product_param(param1, param2, param3);
  EXPECT_EQ(product_param.LocalSize(),
            param1->LocalSize() + param2->LocalSize() + param3->LocalSize());
  EXPECT_EQ(product_param.GlobalSize(),
            param1->GlobalSize() + param2->GlobalSize() + param3->GlobalSize());
}

TEST_F(ProductParameterizationTest, LocalAndGlobalSize4) {
  LocalParameterization* param1 = param1_.release();
  LocalParameterization* param2 = param2_.release();
  LocalParameterization* param3 = param3_.release();
  LocalParameterization* param4 = param4_.release();

  ProductParameterization product_param(param1, param2, param3, param4);
  EXPECT_EQ(product_param.LocalSize(),
            param1->LocalSize() +
            param2->LocalSize() +
            param3->LocalSize() +
            param4->LocalSize());
  EXPECT_EQ(product_param.GlobalSize(),
            param1->GlobalSize() +
            param2->GlobalSize() +
            param3->GlobalSize() +
            param4->GlobalSize());
}

TEST_F(ProductParameterizationTest, Plus) {
  LocalParameterization* param1 = param1_.release();
  LocalParameterization* param2 = param2_.release();
  LocalParameterization* param3 = param3_.release();
  LocalParameterization* param4 = param4_.release();

  ProductParameterization product_param(param1, param2, param3, param4);
  std::vector<double> x(product_param.GlobalSize(), 0.0);
  std::vector<double> delta(product_param.LocalSize(), 0.0);
  std::vector<double> x_plus_delta_expected(product_param.GlobalSize(), 0.0);
  std::vector<double> x_plus_delta(product_param.GlobalSize(), 0.0);

  for (int i = 0; i < product_param.GlobalSize(); ++i) {
    x[i] = RandNormal();
  }

  for (int i = 0; i < product_param.LocalSize(); ++i) {
    delta[i] = RandNormal();
  }

  EXPECT_TRUE(product_param.Plus(&x[0], &delta[0], &x_plus_delta_expected[0]));
  int x_cursor = 0;
  int delta_cursor = 0;

  EXPECT_TRUE(param1->Plus(&x[x_cursor],
                           &delta[delta_cursor],
                           &x_plus_delta[x_cursor]));
  x_cursor += param1->GlobalSize();
  delta_cursor += param1->LocalSize();

  EXPECT_TRUE(param2->Plus(&x[x_cursor],
                           &delta[delta_cursor],
                           &x_plus_delta[x_cursor]));
  x_cursor += param2->GlobalSize();
  delta_cursor += param2->LocalSize();

  EXPECT_TRUE(param3->Plus(&x[x_cursor],
                           &delta[delta_cursor],
                           &x_plus_delta[x_cursor]));
  x_cursor += param3->GlobalSize();
  delta_cursor += param3->LocalSize();

  EXPECT_TRUE(param4->Plus(&x[x_cursor],
                           &delta[delta_cursor],
                           &x_plus_delta[x_cursor]));
  x_cursor += param4->GlobalSize();
  delta_cursor += param4->LocalSize();

  for (int i = 0; i < x.size(); ++i) {
    EXPECT_EQ(x_plus_delta[i], x_plus_delta_expected[i]);
  }
}

TEST_F(ProductParameterizationTest, ComputeJacobian) {
  LocalParameterization* param1 = param1_.release();
  LocalParameterization* param2 = param2_.release();
  LocalParameterization* param3 = param3_.release();
  LocalParameterization* param4 = param4_.release();

  ProductParameterization product_param(param1, param2, param3, param4);
  std::vector<double> x(product_param.GlobalSize(), 0.0);

  for (int i = 0; i < product_param.GlobalSize(); ++i) {
    x[i] = RandNormal();
  }

  Matrix jacobian = Matrix::Random(product_param.GlobalSize(),
                                   product_param.LocalSize());
  EXPECT_TRUE(product_param.ComputeJacobian(&x[0], jacobian.data()));
  int x_cursor = 0;
  int delta_cursor = 0;

  Matrix jacobian1(param1->GlobalSize(), param1->LocalSize());
  EXPECT_TRUE(param1->ComputeJacobian(&x[x_cursor], jacobian1.data()));
  jacobian.block(x_cursor, delta_cursor,
                 param1->GlobalSize(),
                 param1->LocalSize())
      -= jacobian1;
  x_cursor += param1->GlobalSize();
  delta_cursor += param1->LocalSize();

  Matrix jacobian2(param2->GlobalSize(), param2->LocalSize());
  EXPECT_TRUE(param2->ComputeJacobian(&x[x_cursor], jacobian2.data()));
  jacobian.block(x_cursor, delta_cursor,
                 param2->GlobalSize(),
                 param2->LocalSize())
      -= jacobian2;
  x_cursor += param2->GlobalSize();
  delta_cursor += param2->LocalSize();

  Matrix jacobian3(param3->GlobalSize(), param3->LocalSize());
  EXPECT_TRUE(param3->ComputeJacobian(&x[x_cursor], jacobian3.data()));
  jacobian.block(x_cursor, delta_cursor,
                 param3->GlobalSize(),
                 param3->LocalSize())
      -= jacobian3;
  x_cursor += param3->GlobalSize();
  delta_cursor += param3->LocalSize();

  Matrix jacobian4(param4->GlobalSize(), param4->LocalSize());
  EXPECT_TRUE(param4->ComputeJacobian(&x[x_cursor], jacobian4.data()));
  jacobian.block(x_cursor, delta_cursor,
                 param4->GlobalSize(),
                 param4->LocalSize())
      -= jacobian4;
  x_cursor += param4->GlobalSize();
  delta_cursor += param4->LocalSize();

  EXPECT_NEAR(jacobian.norm(), 0.0, std::numeric_limits<double>::epsilon());
}

}  // namespace internal
}  // namespace ceres
