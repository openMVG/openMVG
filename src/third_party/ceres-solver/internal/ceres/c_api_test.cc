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
// Author: mierle@gmail.com (Keir Mierle)

#include "ceres/c_api.h"

#include <cmath>

#include "glog/logging.h"
#include "gtest/gtest.h"

// Duplicated from curve_fitting.cc.
int num_observations = 67;
double data[] = {
  0.000000e+00, 1.133898e+00,
  7.500000e-02, 1.334902e+00,
  1.500000e-01, 1.213546e+00,
  2.250000e-01, 1.252016e+00,
  3.000000e-01, 1.392265e+00,
  3.750000e-01, 1.314458e+00,
  4.500000e-01, 1.472541e+00,
  5.250000e-01, 1.536218e+00,
  6.000000e-01, 1.355679e+00,
  6.750000e-01, 1.463566e+00,
  7.500000e-01, 1.490201e+00,
  8.250000e-01, 1.658699e+00,
  9.000000e-01, 1.067574e+00,
  9.750000e-01, 1.464629e+00,
  1.050000e+00, 1.402653e+00,
  1.125000e+00, 1.713141e+00,
  1.200000e+00, 1.527021e+00,
  1.275000e+00, 1.702632e+00,
  1.350000e+00, 1.423899e+00,
  1.425000e+00, 1.543078e+00,
  1.500000e+00, 1.664015e+00,
  1.575000e+00, 1.732484e+00,
  1.650000e+00, 1.543296e+00,
  1.725000e+00, 1.959523e+00,
  1.800000e+00, 1.685132e+00,
  1.875000e+00, 1.951791e+00,
  1.950000e+00, 2.095346e+00,
  2.025000e+00, 2.361460e+00,
  2.100000e+00, 2.169119e+00,
  2.175000e+00, 2.061745e+00,
  2.250000e+00, 2.178641e+00,
  2.325000e+00, 2.104346e+00,
  2.400000e+00, 2.584470e+00,
  2.475000e+00, 1.914158e+00,
  2.550000e+00, 2.368375e+00,
  2.625000e+00, 2.686125e+00,
  2.700000e+00, 2.712395e+00,
  2.775000e+00, 2.499511e+00,
  2.850000e+00, 2.558897e+00,
  2.925000e+00, 2.309154e+00,
  3.000000e+00, 2.869503e+00,
  3.075000e+00, 3.116645e+00,
  3.150000e+00, 3.094907e+00,
  3.225000e+00, 2.471759e+00,
  3.300000e+00, 3.017131e+00,
  3.375000e+00, 3.232381e+00,
  3.450000e+00, 2.944596e+00,
  3.525000e+00, 3.385343e+00,
  3.600000e+00, 3.199826e+00,
  3.675000e+00, 3.423039e+00,
  3.750000e+00, 3.621552e+00,
  3.825000e+00, 3.559255e+00,
  3.900000e+00, 3.530713e+00,
  3.975000e+00, 3.561766e+00,
  4.050000e+00, 3.544574e+00,
  4.125000e+00, 3.867945e+00,
  4.200000e+00, 4.049776e+00,
  4.275000e+00, 3.885601e+00,
  4.350000e+00, 4.110505e+00,
  4.425000e+00, 4.345320e+00,
  4.500000e+00, 4.161241e+00,
  4.575000e+00, 4.363407e+00,
  4.650000e+00, 4.161576e+00,
  4.725000e+00, 4.619728e+00,
  4.800000e+00, 4.737410e+00,
  4.875000e+00, 4.727863e+00,
  4.950000e+00, 4.669206e+00,
};

// A test cost function, similar to the one in curve_fitting.c.
int exponential_residual(void* user_data,
                         double** parameters,
                         double* residuals,
                         double** jacobians) {
  double* measurement = (double*) user_data;
  double x = measurement[0];
  double y = measurement[1];
  double m = parameters[0][0];
  double c = parameters[1][0];

  residuals[0] = y - exp(m * x + c);
  if (jacobians == NULL) {
    return 1;
  }
  if (jacobians[0] != NULL) {
    jacobians[0][0] = - x * exp(m * x + c);  // dr/dm
  }
  if (jacobians[1] != NULL) {
    jacobians[1][0] =     - exp(m * x + c);  // dr/dc
  }
  return 1;
}

namespace ceres {
namespace internal {

TEST(C_API, SimpleEndToEndTest) {
  double m = 0.0;
  double c = 0.0;
  double *parameter_pointers[] = { &m, &c };
  int parameter_sizes[] = { 1, 1 };

  ceres_problem_t* problem = ceres_create_problem();
  for (int i = 0; i < num_observations; ++i) {
    ceres_problem_add_residual_block(
        problem,
        exponential_residual,  // Cost function
        &data[2 * i],          // Points to the (x,y) measurement
        NULL,                  // Loss function
        NULL,                  // Loss function user data
        1,                     // Number of residuals
        2,                     // Number of parameter blocks
        parameter_sizes,
        parameter_pointers);
  }

  ceres_solve(problem);

  EXPECT_NEAR(0.3, m, 0.02);
  EXPECT_NEAR(0.1, c, 0.04);

  ceres_free_problem(problem);
}

template<typename T>
class ScopedSetValue {
 public:
  ScopedSetValue(T* variable, T new_value)
      : variable_(variable), old_value_(*variable) {
    *variable = new_value;
  }
  ~ScopedSetValue() {
    *variable_ = old_value_;
  }

 private:
  T* variable_;
  T old_value_;
};

TEST(C_API, LossFunctions) {
  double m = 0.2;
  double c = 0.03;
  double *parameter_pointers[] = { &m, &c };
  int parameter_sizes[] = { 1, 1 };

  // Create two outliers, but be careful to leave the data intact.
  ScopedSetValue<double> outlier1x(&data[12], 2.5);
  ScopedSetValue<double> outlier1y(&data[13], 1.0e3);
  ScopedSetValue<double> outlier2x(&data[14], 3.2);
  ScopedSetValue<double> outlier2y(&data[15], 30e3);

  // Create a cauchy cost function, and reuse it many times.
  void* cauchy_loss_data =
      ceres_create_cauchy_loss_function_data(5.0);

  ceres_problem_t* problem = ceres_create_problem();
  for (int i = 0; i < num_observations; ++i) {
    ceres_problem_add_residual_block(
        problem,
        exponential_residual,  // Cost function
        &data[2 * i],          // Points to the (x,y) measurement
        ceres_stock_loss_function,
        cauchy_loss_data,      // Loss function user data
        1,                     // Number of residuals
        2,                     // Number of parameter blocks
        parameter_sizes,
        parameter_pointers);
  }

  ceres_solve(problem);

  EXPECT_NEAR(0.3, m, 0.02);
  EXPECT_NEAR(0.1, c, 0.04);

  ceres_free_stock_loss_function_data(cauchy_loss_data);
  ceres_free_problem(problem);
}

}  // namespace internal
}  // namespace ceres
