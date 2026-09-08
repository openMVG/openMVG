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
// Bicubic interpolation with analytic differentiation
//
// We will use estimation of 2d shift as a sample problem for bicubic
// interpolation.
//
// Let us define f(x, y) = x * x - y * x + y * y
// And optimize cost function sum_i [f(x_i + s_x, y_i + s_y) - v_i]^2
//
// Bicubic interpolation of f(x, y) will be exact, thus we can expect close to
// perfect convergence

#include <utility>

#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include "glog/logging.h"

using Grid = ceres::Grid2D<double>;
using Interpolator = ceres::BiCubicInterpolator<Grid>;

// Cost-function using analytic interface of BiCubicInterpolator
struct AnalyticBiCubicCost : public ceres::CostFunction {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  bool Evaluate(double const* const* parameters,
                double* residuals,
                double** jacobians) const override {
    Eigen::Map<const Eigen::Vector2d> shift(parameters[0]);

    const Eigen::Vector2d point = point_ + shift;

    double* f = residuals;
    double* dfdr = nullptr;
    double* dfdc = nullptr;
    if (jacobians && jacobians[0]) {
      dfdc = jacobians[0];
      dfdr = dfdc + 1;
    }

    interpolator_.Evaluate(point.y(), point.x(), f, dfdr, dfdc);

    if (residuals) {
      *f -= value_;
    }
    return true;
  }

  AnalyticBiCubicCost(const Interpolator& interpolator,
                      Eigen::Vector2d point,
                      double value)
      : point_(std::move(point)), value_(value), interpolator_(interpolator) {
    set_num_residuals(1);
    *mutable_parameter_block_sizes() = {2};
  }

  static ceres::CostFunction* Create(const Interpolator& interpolator,
                                     const Eigen::Vector2d& point,
                                     double value) {
    return new AnalyticBiCubicCost(interpolator, point, value);
  }

  const Eigen::Vector2d point_;
  const double value_;
  const Interpolator& interpolator_;
};

// Function for input data generation
static double f(const double& x, const double& y) {
  return x * x - y * x + y * y;
}

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  // Problem sizes
  const int kGridRowsHalf = 9;
  const int kGridColsHalf = 11;
  const int kGridRows = 2 * kGridRowsHalf + 1;
  const int kGridCols = 2 * kGridColsHalf + 1;
  const int kPoints = 4;

  const Eigen::Vector2d shift(1.234, 2.345);
  const std::array<Eigen::Vector2d, kPoints> points = {
      Eigen::Vector2d{-2., -3.},
      Eigen::Vector2d{-2., 3.},
      Eigen::Vector2d{2., 3.},
      Eigen::Vector2d{2., -3.}};

  // Data is a row-major array of kGridRows x kGridCols values of function
  // f(x, y) on the grid, with x in {-kGridColsHalf, ..., +kGridColsHalf},
  // and y in {-kGridRowsHalf, ..., +kGridRowsHalf}
  double data[kGridRows * kGridCols];
  for (int i = 0; i < kGridRows; ++i) {
    for (int j = 0; j < kGridCols; ++j) {
      // Using row-major order
      int index = i * kGridCols + j;
      double y = i - kGridRowsHalf;
      double x = j - kGridColsHalf;

      data[index] = f(x, y);
    }
  }
  const Grid grid(data,
                  -kGridRowsHalf,
                  kGridRowsHalf + 1,
                  -kGridColsHalf,
                  kGridColsHalf + 1);
  const Interpolator interpolator(grid);

  Eigen::Vector2d shift_estimate(3.1415, 1.337);

  ceres::Problem problem;
  problem.AddParameterBlock(shift_estimate.data(), 2);

  for (const auto& p : points) {
    const Eigen::Vector2d shifted = p + shift;

    const double v = f(shifted.x(), shifted.y());
    problem.AddResidualBlock(AnalyticBiCubicCost::Create(interpolator, p, v),
                             nullptr,
                             shift_estimate.data());
  }

  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << '\n';

  std::cout << "Bicubic interpolation with analytic derivatives:\n";
  std::cout << "Estimated shift: " << shift_estimate.transpose()
            << ", ground-truth: " << shift.transpose()
            << " (error: " << (shift_estimate - shift).transpose() << ")"
            << std::endl;

  CHECK_LT((shift_estimate - shift).norm(), 1e-9);
  return 0;
}
