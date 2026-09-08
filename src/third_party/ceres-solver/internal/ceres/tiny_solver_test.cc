
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
// Author: mierle@gmail.com (Keir Mierle)

#include "ceres/tiny_solver.h"

#include <algorithm>
#include <cmath>

#include "ceres/tiny_solver_test_util.h"
#include "gtest/gtest.h"

namespace ceres {

using Vec2 = Eigen::Matrix<double, 2, 1>;
using Vec3 = Eigen::Matrix<double, 3, 1>;
using VecX = Eigen::VectorXd;

class ExampleStatic {
 public:
  using Scalar = double;
  enum {
    // Can also be Eigen::Dynamic.
    NUM_RESIDUALS = 2,
    NUM_PARAMETERS = 3,
  };
  bool operator()(const double* parameters,
                  double* residuals,
                  double* jacobian) const {
    return EvaluateResidualsAndJacobians(parameters, residuals, jacobian);
  }
};

class ExampleParametersDynamic {
 public:
  using Scalar = double;
  enum {
    NUM_RESIDUALS = 2,
    NUM_PARAMETERS = Eigen::Dynamic,
  };

  int NumParameters() const { return 3; }

  bool operator()(const double* parameters,
                  double* residuals,
                  double* jacobian) const {
    return EvaluateResidualsAndJacobians(parameters, residuals, jacobian);
  }
};

class ExampleResidualsDynamic {
 public:
  using Scalar = double;
  enum {
    NUM_RESIDUALS = Eigen::Dynamic,
    NUM_PARAMETERS = 3,
  };

  int NumResiduals() const { return 2; }

  bool operator()(const double* parameters,
                  double* residuals,
                  double* jacobian) const {
    return EvaluateResidualsAndJacobians(parameters, residuals, jacobian);
  }
};

class ExampleAllDynamic {
 public:
  using Scalar = double;
  enum {
    NUM_RESIDUALS = Eigen::Dynamic,
    NUM_PARAMETERS = Eigen::Dynamic,
  };

  int NumResiduals() const { return 2; }

  int NumParameters() const { return 3; }

  bool operator()(const double* parameters,
                  double* residuals,
                  double* jacobian) const {
    return EvaluateResidualsAndJacobians(parameters, residuals, jacobian);
  }
};

template <typename Function, typename Vector>
void TestHelper(const Function& f, const Vector& x0) {
  Vector x = x0;
  Vec2 residuals;
  f(x.data(), residuals.data(), nullptr);
  EXPECT_GT(residuals.squaredNorm() / 2.0, 1e-10);

  TinySolver<Function> solver;
  solver.Solve(f, &x);
  EXPECT_NEAR(0.0, solver.summary.final_cost, 1e-10);
}

// A test case for when the cost function is statically sized.
TEST(TinySolver, SimpleExample) {
  Vec3 x0(0.76026643, -30.01799744, 0.55192142);
  ExampleStatic f;

  TestHelper(f, x0);
}

// A test case for when the number of parameters is dynamically sized.
TEST(TinySolver, ParametersDynamic) {
  VecX x0(3);
  x0 << 0.76026643, -30.01799744, 0.55192142;

  ExampleParametersDynamic f;

  TestHelper(f, x0);
}

// A test case for when the number of residuals is dynamically sized.
TEST(TinySolver, ResidualsDynamic) {
  Vec3 x0(0.76026643, -30.01799744, 0.55192142);

  ExampleResidualsDynamic f;

  TestHelper(f, x0);
}

// A test case for when the number of parameters and residuals is
// dynamically sized.
TEST(TinySolver, ParametersAndResidualsDynamic) {
  VecX x0(3);
  x0 << 0.76026643, -30.01799744, 0.55192142;

  ExampleAllDynamic f;

  TestHelper(f, x0);
}

}  // namespace ceres
