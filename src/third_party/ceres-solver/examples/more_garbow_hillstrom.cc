// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2014 Google Inc. All rights reserved.
// http://code.google.com/p/ceres-solver/
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
//
// Bounds constrained test problems from the paper
//
// Testing Unconstrained Optimization Software
// Jorge J. More, Burton S. Garbow and Kenneth E. Hillstrom
// ACM Transactions on Mathematical Software, 7(1), pp. 17-41, 1981
//
// A subset of these problems were augmented with bounds and used for
// testing bounds constrained optimization algorithms by
//
// A Trust Region Approach to Linearly Constrained Optimization
// David M. Gay
// Numerical Analysis (Griffiths, D.F., ed.), pp. 72-105
// Lecture Notes in Mathematics 1066, Springer Verlag, 1984.
//
// The latter paper is behind a paywall. We obtained the bounds on the
// variables and the function values at the global minimums from
//
// http://www.mat.univie.ac.at/~neum/glopt/bounds.html
//
// A problem is considered solved if of the log relative error of its
// objective function is at least 5.


#include <cmath>
#include <iostream>  // NOLINT
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

namespace ceres {
namespace examples {

const double kDoubleMax = std::numeric_limits<double>::max();

#define BEGIN_MGH_PROBLEM(name, num_parameters, num_residuals)          \
  struct name {                                                         \
    static const int kNumParameters = num_parameters;                   \
    static const double initial_x[kNumParameters];                      \
    static const double lower_bounds[kNumParameters];                   \
    static const double upper_bounds[kNumParameters];                   \
    static const double constrained_optimal_cost;                       \
    static const double unconstrained_optimal_cost;                     \
    static CostFunction* Create() {                                     \
      return new AutoDiffCostFunction<name,                             \
                                      num_residuals,                    \
                                      num_parameters>(new name);        \
    }                                                                   \
    template <typename T>                                               \
    bool operator()(const T* const x, T* residual) const {

#define END_MGH_PROBLEM return true; } };  // NOLINT

// Rosenbrock function.
BEGIN_MGH_PROBLEM(TestProblem1, 2, 2)
  const T x1 = x[0];
  const T x2 = x[1];
  residual[0] = T(10.0) * (x2 - x1 * x1);
  residual[1] = T(1.0) - x1;
END_MGH_PROBLEM;

const double TestProblem1::initial_x[] = {-1.2, 1.0};
const double TestProblem1::lower_bounds[] = {-kDoubleMax, -kDoubleMax};
const double TestProblem1::upper_bounds[] = {kDoubleMax, kDoubleMax};
const double TestProblem1::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem1::unconstrained_optimal_cost = 0.0;

// Freudenstein and Roth function.
BEGIN_MGH_PROBLEM(TestProblem2, 2, 2)
  const T x1 = x[0];
  const T x2 = x[1];
  residual[0] = T(-13.0) + x1 + ((T(5.0) - x2) * x2 - T(2.0)) * x2;
  residual[1] = T(-29.0) + x1 + ((x2 + T(1.0)) * x2 - T(14.0)) * x2;
END_MGH_PROBLEM;

const double TestProblem2::initial_x[] = {0.5, -2.0};
const double TestProblem2::lower_bounds[] = {-kDoubleMax, -kDoubleMax};
const double TestProblem2::upper_bounds[] = {kDoubleMax, kDoubleMax};
const double TestProblem2::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem2::unconstrained_optimal_cost = 0.0;

// Powell badly scaled function.
BEGIN_MGH_PROBLEM(TestProblem3, 2, 2)
  const T x1 = x[0];
  const T x2 = x[1];
  residual[0] = T(10000.0) * x1 * x2 - T(1.0);
  residual[1] = exp(-x1) + exp(-x2) - T(1.0001);
END_MGH_PROBLEM;

const double TestProblem3::initial_x[] = {0.0, 1.0};
const double TestProblem3::lower_bounds[] = {0.0, 1.0};
const double TestProblem3::upper_bounds[] = {1.0, 9.0};
const double TestProblem3::constrained_optimal_cost = 0.15125900e-9;
const double TestProblem3::unconstrained_optimal_cost = 0.0;

// Brown badly scaled function.
BEGIN_MGH_PROBLEM(TestProblem4, 2, 3)
  const T x1 = x[0];
  const T x2 = x[1];
  residual[0] = x1  - T(1000000.0);
  residual[1] = x2 - T(0.000002);
  residual[2] = x1 * x2 - T(2.0);
END_MGH_PROBLEM;

const double TestProblem4::initial_x[] = {1.0, 1.0};
const double TestProblem4::lower_bounds[] = {0.0, 0.00003};
const double TestProblem4::upper_bounds[] = {1000000.0, 100.0};
const double TestProblem4::constrained_optimal_cost = 0.78400000e3;
const double TestProblem4::unconstrained_optimal_cost = 0.0;

// Beale function.
BEGIN_MGH_PROBLEM(TestProblem5, 2, 3)
  const T x1 = x[0];
  const T x2 = x[1];
  residual[0] = T(1.5) - x1 * (T(1.0) - x2);
  residual[1] = T(2.25) - x1 * (T(1.0) - x2 * x2);
  residual[2] = T(2.625) - x1 * (T(1.0) - x2 * x2 * x2);
END_MGH_PROBLEM;

const double TestProblem5::initial_x[] = {1.0, 1.0};
const double TestProblem5::lower_bounds[] = {0.6, 0.5};
const double TestProblem5::upper_bounds[] = {10.0, 100.0};
const double TestProblem5::constrained_optimal_cost = 0.0;
const double TestProblem5::unconstrained_optimal_cost = 0.0;

// Jennrich and Sampson function.
BEGIN_MGH_PROBLEM(TestProblem6, 2, 10)
  const T x1 = x[0];
  const T x2 = x[1];
  for (int i = 1; i <= 10; ++i) {
    residual[i - 1] = T(2.0) + T(2.0 * i) -
        exp(T(static_cast<double>(i)) * x1) -
        exp(T(static_cast<double>(i) * x2));
  }
END_MGH_PROBLEM;

const double TestProblem6::initial_x[] = {1.0, 1.0};
const double TestProblem6::lower_bounds[] = {-kDoubleMax, -kDoubleMax};
const double TestProblem6::upper_bounds[] = {kDoubleMax, kDoubleMax};
const double TestProblem6::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem6::unconstrained_optimal_cost = 124.362;

// Helical valley function.
BEGIN_MGH_PROBLEM(TestProblem7, 3, 3)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T theta = T(0.5 / M_PI)  * atan(x2 / x1) + (x1 > 0.0 ? T(0.0) : T(0.5));

  residual[0] = T(10.0) * (x3 - T(10.0) * theta);
  residual[1] = T(10.0) * (sqrt(x1 * x1 + x2 * x2) - T(1.0));
  residual[2] = x3;
END_MGH_PROBLEM;

const double TestProblem7::initial_x[] = {-1.0, 0.0, 0.0};
const double TestProblem7::lower_bounds[] = {-100.0, -1.0, -1.0};
const double TestProblem7::upper_bounds[] = {0.8, 1.0, 1.0};
const double TestProblem7::constrained_optimal_cost = 0.99042212;
const double TestProblem7::unconstrained_optimal_cost = 0.0;

// Bard function
BEGIN_MGH_PROBLEM(TestProblem8, 3, 15)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];

  double y[] = {0.14, 0.18, 0.22, 0.25,
                0.29, 0.32, 0.35, 0.39, 0.37, 0.58,
                0.73, 0.96, 1.34, 2.10, 4.39};

  for (int i = 1; i <=15; ++i) {
    const T u = T(static_cast<double>(i));
    const T v = T(static_cast<double>(16 - i));
    const T w = T(static_cast<double>(std::min(i, 16 - i)));
    residual[i - 1] = T(y[i - 1]) - x1 + u / (v * x2 + w * x3);
  }
END_MGH_PROBLEM;

const double TestProblem8::initial_x[] = {1.0, 1.0, 1.0};
const double TestProblem8::lower_bounds[] = {
  -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem8::upper_bounds[] = {
  kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem8::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem8::unconstrained_optimal_cost = 8.21487e-3;

// Gaussian function.
BEGIN_MGH_PROBLEM(TestProblem9, 3, 15)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];

  const double y[] = {0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521,
                      0.3989,
                      0.3521, 0.2420, 0.1295, 0.0540, 0.0175, 0.0044, 0.0009};
  for (int i = 0; i < 15; ++i) {
    const T t_i = T((8.0 - i - 1.0) / 2.0);
    const T y_i = T(y[i]);
    residual[i] = x1 * exp(-x2 * (t_i - x3) * (t_i - x3) / T(2.0)) - y_i;
  }
END_MGH_PROBLEM;

const double TestProblem9::initial_x[] = {0.4, 1.0, 0.0};
const double TestProblem9::lower_bounds[] = {0.398, 1.0, -0.5};
const double TestProblem9::upper_bounds[] = {4.2, 2.0, 0.1};
const double TestProblem9::constrained_optimal_cost = 0.11279300e-7;
const double TestProblem9::unconstrained_optimal_cost = 0.112793e-7;

// Meyer function.
BEGIN_MGH_PROBLEM(TestProblem10, 3, 16)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];

  const double y[] = {34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744,
                      8261, 7030, 6005, 5147, 4427, 3820, 3307, 2872};

  for (int i = 0; i < 16; ++i) {
    T t = T(45 + 5.0 * (i + 1));
    residual[i] = x1 * exp(x2 / (t + x3)) - y[i];
  }
END_MGH_PROBLEM


const double TestProblem10::initial_x[] = {0.02, 4000, 250};
const double TestProblem10::lower_bounds[] ={
  -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem10::upper_bounds[] ={
  kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem10::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem10::unconstrained_optimal_cost = 87.9458;

#undef BEGIN_MGH_PROBLEM
#undef END_MGH_PROBLEM

template<typename TestProblem> string ConstrainedSolve() {
  double x[TestProblem::kNumParameters];
  std::copy(TestProblem::initial_x,
            TestProblem::initial_x + TestProblem::kNumParameters,
            x);

  Problem problem;
  problem.AddResidualBlock(TestProblem::Create(), NULL, x);
  for (int i = 0; i < TestProblem::kNumParameters; ++i) {
    problem.SetParameterLowerBound(x, i, TestProblem::lower_bounds[i]);
    problem.SetParameterUpperBound(x, i, TestProblem::upper_bounds[i]);
  }

  Solver::Options options;
  options.parameter_tolerance = 1e-18;
  options.function_tolerance = 1e-18;
  options.gradient_tolerance = 1e-18;
  options.max_num_iterations = 1000;
  options.linear_solver_type = DENSE_QR;
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  const double kMinLogRelativeError = 5.0;
  const double log_relative_error = -std::log10(
      std::abs(2.0 * summary.final_cost -
               TestProblem::constrained_optimal_cost) /
      (TestProblem::constrained_optimal_cost > 0.0
       ? TestProblem::constrained_optimal_cost
       : 1.0));

  return (log_relative_error >= kMinLogRelativeError
          ? "Success\n"
          : "Failure\n");
}

template<typename TestProblem> string UnconstrainedSolve() {
  double x[TestProblem::kNumParameters];
  std::copy(TestProblem::initial_x,
            TestProblem::initial_x + TestProblem::kNumParameters,
            x);

  Problem problem;
  problem.AddResidualBlock(TestProblem::Create(), NULL, x);

  Solver::Options options;
  options.parameter_tolerance = 1e-18;
  options.function_tolerance = 0.0;
  options.gradient_tolerance = 1e-18;
  options.max_num_iterations = 1000;
  options.linear_solver_type = DENSE_QR;
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  const double kMinLogRelativeError = 5.0;
  const double log_relative_error = -std::log10(
      std::abs(2.0 * summary.final_cost -
               TestProblem::unconstrained_optimal_cost) /
      (TestProblem::unconstrained_optimal_cost > 0.0
       ? TestProblem::unconstrained_optimal_cost
       : 1.0));

  return (log_relative_error >= kMinLogRelativeError
          ? "Success\n"
          : "Failure\n");
}

}  // namespace examples
}  // namespace ceres

int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  using ceres::examples::UnconstrainedSolve;
  using ceres::examples::ConstrainedSolve;

#define UNCONSTRAINED_SOLVE(n)                                          \
  std::cout << "Problem " << n << " : "                                 \
            << UnconstrainedSolve<ceres::examples::TestProblem##n>();

#define CONSTRAINED_SOLVE(n)                                            \
  std::cout << "Problem " << n << " : "                                 \
            << ConstrainedSolve<ceres::examples::TestProblem##n>();

  std::cout << "Unconstrained problems\n";
  UNCONSTRAINED_SOLVE(1);
  UNCONSTRAINED_SOLVE(2);
  UNCONSTRAINED_SOLVE(3);
  UNCONSTRAINED_SOLVE(4);
  UNCONSTRAINED_SOLVE(5);
  UNCONSTRAINED_SOLVE(6);
  UNCONSTRAINED_SOLVE(7);
  UNCONSTRAINED_SOLVE(8);
  UNCONSTRAINED_SOLVE(9);
  UNCONSTRAINED_SOLVE(10);

  std::cout << "\nConstrained problems\n";
  CONSTRAINED_SOLVE(3);
  CONSTRAINED_SOLVE(4);
  CONSTRAINED_SOLVE(5);
  CONSTRAINED_SOLVE(7);
  CONSTRAINED_SOLVE(9);

  return 0;
}
