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
//
// Test problems from the paper
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
// objective function is at least 4.


#include <cmath>
#include <iostream>  // NOLINT
#include <sstream>   // NOLINT
#include <string>
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

DEFINE_string(problem, "all", "Which problem to solve");
DEFINE_bool(use_numeric_diff, false,
            "Use numeric differentiation instead of automatic "
            "differentiation.");
DEFINE_string(numeric_diff_method, "ridders", "When using numeric "
              "differentiation, selects algorithm. Options are: central, "
              "forward, ridders.");
DEFINE_int32(ridders_extrapolations, 3, "Maximal number of extrapolations in "
             "Ridders' method.");

namespace ceres {
namespace examples {

const double kDoubleMax = std::numeric_limits<double>::max();

static void SetNumericDiffOptions(ceres::NumericDiffOptions* options) {
  options->max_num_ridders_extrapolations = FLAGS_ridders_extrapolations;
}

#define BEGIN_MGH_PROBLEM(name, num_parameters, num_residuals)            \
  struct name {                                                           \
    static const int kNumParameters = num_parameters;                     \
    static const double initial_x[kNumParameters];                        \
    static const double lower_bounds[kNumParameters];                     \
    static const double upper_bounds[kNumParameters];                     \
    static const double constrained_optimal_cost;                         \
    static const double unconstrained_optimal_cost;                       \
    static CostFunction* Create() {                                       \
      if (FLAGS_use_numeric_diff) {                                       \
        ceres::NumericDiffOptions options;                                \
        SetNumericDiffOptions(&options);                                  \
        if (FLAGS_numeric_diff_method == "central") {                     \
          return new NumericDiffCostFunction<name,                        \
                                             ceres::CENTRAL,              \
                                             num_residuals,               \
                                             num_parameters>(             \
              new name, ceres::TAKE_OWNERSHIP, num_residuals, options);   \
        } else if (FLAGS_numeric_diff_method == "forward") {              \
          return new NumericDiffCostFunction<name,                        \
                                             ceres::FORWARD,              \
                                             num_residuals,               \
                                             num_parameters>(             \
              new name, ceres::TAKE_OWNERSHIP, num_residuals, options);   \
        } else if (FLAGS_numeric_diff_method == "ridders") {              \
          return new NumericDiffCostFunction<name,                        \
                                             ceres::RIDDERS,              \
                                             num_residuals,               \
                                             num_parameters>(             \
              new name, ceres::TAKE_OWNERSHIP, num_residuals, options);   \
        } else {                                                          \
          LOG(ERROR) << "Invalid numeric diff method specified";          \
          return NULL;                                                    \
        }                                                                 \
      } else {                                                            \
        return new AutoDiffCostFunction<name,                             \
                                        num_residuals,                    \
                                        num_parameters>(new name);        \
      }                                                                   \
    }                                                                     \
    template <typename T>                                                 \
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
        (exp(T(static_cast<double>(i)) * x1) +
         exp(T(static_cast<double>(i) * x2)));
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
    residual[i - 1] = T(y[i - 1]) - (x1 + u / (v * x2 + w * x3));
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
    const T ti = T(45 + 5.0 * (i + 1));
    const T yi = T(y[i]);
    residual[i] = x1 * exp(x2 / (ti + x3)) - yi;
  }
END_MGH_PROBLEM

const double TestProblem10::initial_x[] = {0.02, 4000, 250};
const double TestProblem10::lower_bounds[] = {
  -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem10::upper_bounds[] = {
  kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem10::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem10::unconstrained_optimal_cost = 87.9458;

// Gulf research and development function
BEGIN_MGH_PROBLEM(TestProblem11, 3, 100)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  for (int i = 1; i <= 100; ++i) {
    const double ti = static_cast<double>(i) / 100.0;
    const double yi = 25.0 + pow(-50.0 * log(ti), 2.0 / 3.0);
    residual[i - 1] = exp(-pow(abs(T(yi * 100.0 * i) * x2), x3) / x1) - T(ti);
  }
END_MGH_PROBLEM

const double TestProblem11::initial_x[] = {5.0, 2.5, 0.15};
const double TestProblem11::lower_bounds[] = {1e-16, 0.0, 0.0};
const double TestProblem11::upper_bounds[] = {10.0, 10.0, 10.0};
const double TestProblem11::constrained_optimal_cost = 0.58281431e-4;
const double TestProblem11::unconstrained_optimal_cost = 0.0;

// Box three-dimensional function.
BEGIN_MGH_PROBLEM(TestProblem12, 3, 3)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];

  const T t1 = T(0.1);
  const T t2 = T(0.2);
  const T t3 = T(0.3);

  residual[0] = exp(-t1 * x1) - exp(-t1 * x2) - x3 * (exp(-t1) - exp(-T(10.0) * t1));
  residual[1] = exp(-t2 * x1) - exp(-t2 * x2) - x3 * (exp(-t2) - exp(-T(10.0) * t2));
  residual[2] = exp(-t3 * x1) - exp(-t3 * x2) - x3 * (exp(-t3) - exp(-T(10.0) * t3));
END_MGH_PROBLEM

const double TestProblem12::initial_x[] = {0.0, 10.0, 20.0};
const double TestProblem12::lower_bounds[] = {0.0, 5.0, 0.0};
const double TestProblem12::upper_bounds[] = {2.0, 9.5, 20.0};
const double TestProblem12::constrained_optimal_cost = 0.30998153e-5;
const double TestProblem12::unconstrained_optimal_cost = 0.0;

// Powell Singular function.
BEGIN_MGH_PROBLEM(TestProblem13, 4, 4)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];

  residual[0] = x1 + T(10.0) * x2;
  residual[1] = T(sqrt(5.0)) * (x3 - x4);
  residual[2] = (x2 - T(2.0) * x3) * (x2 - T(2.0) * x3);
  residual[3] = sqrt(10.0) * (x1 - x4) * (x1 - x4);
END_MGH_PROBLEM

const double TestProblem13::initial_x[] = {3.0, -1.0, 0.0, 1.0};
const double TestProblem13::lower_bounds[] = {
  -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem13::upper_bounds[] = {
  kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem13::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem13::unconstrained_optimal_cost = 0.0;

// Wood function.
BEGIN_MGH_PROBLEM(TestProblem14, 4, 6)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];

  residual[0] = T(10.0) * (x2 - x1 * x1);
  residual[1] = T(1.0) - x1;
  residual[2] = T(sqrt(90.0)) * (x4 - x3 * x3);
  residual[3] = T(1.0) - x3;
  residual[4] = T(sqrt(10.0)) * (x2 + x4 - T(2.0));
  residual[5] = T(1.0/sqrt(10.0)) * (x2 - x4);
END_MGH_PROBLEM;

const double TestProblem14::initial_x[] = {-3.0, -1.0, -3.0, -1.0};
const double TestProblem14::lower_bounds[] = {-100.0, -100.0, -100.0, -100.0};
const double TestProblem14::upper_bounds[] = {0.0, 10.0, 100.0, 100.0};
const double TestProblem14::constrained_optimal_cost = 0.15567008e1;
const double TestProblem14::unconstrained_optimal_cost = 0.0;

// Kowalik and Osborne function.
BEGIN_MGH_PROBLEM(TestProblem15, 4, 11)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];

  const double y[] = {0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627,
                      0.0456, 0.0342, 0.0323, 0.0235, 0.0246};
  const double u[] = {4.0, 2.0, 1.0, 0.5, 0.25, 0.167, 0.125, 0.1,
                      0.0833, 0.0714, 0.0625};

  for (int i = 0; i < 11; ++i) {
    const T yi = T(y[i]);
    const T ui = T(u[i]);
    residual[i]  = yi - x1 * (ui * ui + ui * x2) / (ui * ui  + ui * x3 + x4);
  }
END_MGH_PROBLEM;

const double TestProblem15::initial_x[] = {0.25, 0.39, 0.415, 0.39};
const double TestProblem15::lower_bounds[] = {
  -kDoubleMax, -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem15::upper_bounds[] = {
  kDoubleMax, kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem15::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem15::unconstrained_optimal_cost = 3.07505e-4;

// Brown and Dennis function.
BEGIN_MGH_PROBLEM(TestProblem16, 4, 20)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];

  for (int i = 0; i < 20; ++i) {
    const T ti = T(static_cast<double>(i + 1) / 5.0);
    residual[i] = (x1 + ti * x2 - exp(ti)) * (x1 + ti * x2 - exp(ti)) +
        (x3 + x4 * sin(ti) - cos(ti)) * (x3 + x4 * sin(ti) - cos(ti));
  }
END_MGH_PROBLEM;

const double TestProblem16::initial_x[] = {25.0, 5.0, -5.0, -1.0};
const double TestProblem16::lower_bounds[] = {-10.0, 0.0, -100.0, -20.0};
const double TestProblem16::upper_bounds[] = {100.0, 15.0, 0.0, 0.2};
const double TestProblem16::constrained_optimal_cost = 0.88860479e5;
const double TestProblem16::unconstrained_optimal_cost = 85822.2;

// Osborne 1 function.
BEGIN_MGH_PROBLEM(TestProblem17, 5, 33)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];
  const T x5 = x[4];

  const double y[] = {0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850, 0.818,
                      0.784, 0.751, 0.718, 0.685, 0.658, 0.628, 0.603, 0.580, 0.558,
                      0.538, 0.522, 0.506, 0.490, 0.478, 0.467, 0.457, 0.448, 0.438,
                      0.431, 0.424, 0.420, 0.414, 0.411, 0.406};

  for (int i = 0; i < 33; ++i) {
    const T yi = T(y[i]);
    const T ti = T(10.0 * i);
    residual[i] = yi - (x1 + x2 * exp(-ti * x4) + x3 * exp(-ti * x5));
  }
END_MGH_PROBLEM;

const double TestProblem17::initial_x[] = {0.5, 1.5, -1.0, 0.01, 0.02};
const double TestProblem17::lower_bounds[] = {
  -kDoubleMax, -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem17::upper_bounds[] = {
  kDoubleMax, kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem17::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem17::unconstrained_optimal_cost = 5.46489e-5;

// Biggs EXP6 function.
BEGIN_MGH_PROBLEM(TestProblem18, 6, 13)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];
  const T x5 = x[4];
  const T x6 = x[5];

  for (int i = 0; i < 13; ++i) {
    const double ti = 0.1 * (i + 1.0);
    const double yi = exp(-ti) - 5.0 * exp(-10.0 * ti) + 3.0 * exp(-4.0 * ti);
    const T si = T(ti);
    residual[i] =x3 * exp(-si * x1) - x4 * exp(-si * x2) + x6 * exp(-si * x5) - T(yi);
  }
END_MGH_PROBLEM

const double TestProblem18::initial_x[] = {1.0, 2.0, 1.0, 1.0, 1.0, 1.0};
const double TestProblem18::lower_bounds[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
const double TestProblem18::upper_bounds[] = {2.0, 8.0, 1.0, 7.0, 5.0, 5.0};
const double TestProblem18::constrained_optimal_cost = 0.53209865e-3;
const double TestProblem18::unconstrained_optimal_cost = 0.0;

// Osborne 2 function.
BEGIN_MGH_PROBLEM(TestProblem19, 11, 65)
  const T x1 = x[0];
  const T x2 = x[1];
  const T x3 = x[2];
  const T x4 = x[3];
  const T x5 = x[4];
  const T x6 = x[5];
  const T x7 = x[6];
  const T x8 = x[7];
  const T x9 = x[8];
  const T x10 = x[9];
  const T x11 = x[10];

  const double y[] = {1.366, 1.191, 1.112, 1.013, 0.991,
                      0.885, 0.831, 0.847, 0.786, 0.725,
                      0.746, 0.679, 0.608, 0.655, 0.616,
                      0.606, 0.602, 0.626, 0.651, 0.724,
                      0.649, 0.649, 0.694, 0.644, 0.624,
                      0.661, 0.612, 0.558, 0.533, 0.495,
                      0.500, 0.423, 0.395, 0.375, 0.372,
                      0.391, 0.396, 0.405, 0.428, 0.429,
                      0.523, 0.562, 0.607, 0.653, 0.672,
                      0.708, 0.633, 0.668, 0.645, 0.632,
                      0.591, 0.559, 0.597, 0.625, 0.739,
                      0.710, 0.729, 0.720, 0.636, 0.581,
                      0.428, 0.292, 0.162, 0.098, 0.054};

  for (int i = 0; i < 65; ++i) {
    const T ti = T(static_cast<double>(i) / 10.0);
    residual[i] = T(y[i]) - (x1 * exp(-(ti * x5)) +
                             x2 * exp(-(ti - x9)  * (ti - x9)  * x6) +
                             x3 * exp(-(ti - x10) * (ti - x10) * x7) +
                             x4 * exp(-(ti - x11) * (ti - x11) * x8));
  }
END_MGH_PROBLEM;

const double TestProblem19::initial_x[] = {1.3, 0.65, 0.65, 0.7, 0.6,
                                           3.0, 5.0, 7.0, 2.0, 4.5, 5.5};
const double TestProblem19::lower_bounds[] = {
  -kDoubleMax, -kDoubleMax, -kDoubleMax, -kDoubleMax};
const double TestProblem19::upper_bounds[] = {
  kDoubleMax, kDoubleMax, kDoubleMax, kDoubleMax};
const double TestProblem19::constrained_optimal_cost =
    std::numeric_limits<double>::quiet_NaN();
const double TestProblem19::unconstrained_optimal_cost = 4.01377e-2;


#undef BEGIN_MGH_PROBLEM
#undef END_MGH_PROBLEM

template<typename TestProblem> bool Solve(bool is_constrained, int trial) {
  double x[TestProblem::kNumParameters];
  for (int i = 0; i < TestProblem::kNumParameters; ++i) {
    x[i] = pow(10, trial) * TestProblem::initial_x[i];
  }

  Problem problem;
  problem.AddResidualBlock(TestProblem::Create(), NULL, x);
  double optimal_cost = TestProblem::unconstrained_optimal_cost;

  if (is_constrained) {
    for (int i = 0; i < TestProblem::kNumParameters; ++i) {
      problem.SetParameterLowerBound(x, i, TestProblem::lower_bounds[i]);
      problem.SetParameterUpperBound(x, i, TestProblem::upper_bounds[i]);
    }
    optimal_cost = TestProblem::constrained_optimal_cost;
  }

  Solver::Options options;
  options.parameter_tolerance = 1e-18;
  options.function_tolerance = 1e-18;
  options.gradient_tolerance = 1e-18;
  options.max_num_iterations = 1000;
  options.linear_solver_type = DENSE_QR;
  Solver::Summary summary;
  Solve(options, &problem, &summary);

  const double kMinLogRelativeError = 4.0;
  const double log_relative_error = -std::log10(
      std::abs(2.0 * summary.final_cost - optimal_cost) /
      (optimal_cost > 0.0 ? optimal_cost : 1.0));

  const bool success = log_relative_error >= kMinLogRelativeError;
  LOG(INFO) << "Expected : " <<  optimal_cost
            << " actual: " << 2.0 * summary.final_cost
            << " " << success
            << " in " << summary.total_time_in_seconds
            << " seconds";
  return success;
}

}  // namespace examples
}  // namespace ceres

int main(int argc, char** argv) {
  CERES_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  using ceres::examples::Solve;

  int unconstrained_problems = 0;
  int unconstrained_successes = 0;
  int constrained_problems = 0;
  int constrained_successes = 0;
  std::stringstream ss;

#define UNCONSTRAINED_SOLVE(n)                                          \
  ss << "Unconstrained Problem " << n << " : ";                          \
  if (FLAGS_problem == #n || FLAGS_problem == "all") {                  \
    unconstrained_problems += 3;                                        \
    if (Solve<ceres::examples::TestProblem##n>(false, 0)) {             \
      unconstrained_successes += 1;                                     \
      ss <<  "Yes ";                                                    \
    } else {                                                            \
      ss << "No  ";                                                     \
    }                                                                   \
    if (Solve<ceres::examples::TestProblem##n>(false, 1)) {             \
      unconstrained_successes += 1;                                     \
      ss << "Yes ";                                                     \
    } else {                                                            \
      ss << "No  ";                                                     \
    }                                                                   \
    if (Solve<ceres::examples::TestProblem##n>(false, 2)) {             \
      unconstrained_successes += 1;                                     \
      ss << "Yes ";                                                     \
    } else {                                                            \
      ss << "No  ";                                                     \
    }                                                                   \
  }                                                                     \
  ss << std::endl;

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
  UNCONSTRAINED_SOLVE(11);
  UNCONSTRAINED_SOLVE(12);
  UNCONSTRAINED_SOLVE(13);
  UNCONSTRAINED_SOLVE(14);
  UNCONSTRAINED_SOLVE(15);
  UNCONSTRAINED_SOLVE(16);
  UNCONSTRAINED_SOLVE(17);
  UNCONSTRAINED_SOLVE(18);
  UNCONSTRAINED_SOLVE(19);

  ss << "Unconstrained : "
     << unconstrained_successes
     << "/"
     << unconstrained_problems << std::endl;

#define CONSTRAINED_SOLVE(n)                                            \
  ss << "Constrained Problem " << n << " : ";                           \
  if (FLAGS_problem == #n || FLAGS_problem == "all") {                  \
    constrained_problems += 1;                                          \
    if (Solve<ceres::examples::TestProblem##n>(true, 0)) {              \
      constrained_successes += 1;                                       \
      ss << "Yes ";                                                     \
    } else {                                                            \
      ss << "No  ";                                                     \
    }                                                                   \
  }                                                                     \
  ss << std::endl;

  CONSTRAINED_SOLVE(3);
  CONSTRAINED_SOLVE(4);
  CONSTRAINED_SOLVE(5);
  CONSTRAINED_SOLVE(7);
  CONSTRAINED_SOLVE(9);
  CONSTRAINED_SOLVE(11);
  CONSTRAINED_SOLVE(12);
  CONSTRAINED_SOLVE(14);
  CONSTRAINED_SOLVE(16);
  CONSTRAINED_SOLVE(18);
  ss << "Constrained : "
     << constrained_successes
     << "/"
     << constrained_problems << std::endl;

  std::cout << ss.str();
  return 0;
}
