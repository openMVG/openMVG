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
// A simple example of optimizing a sampled function by using cubic
// interpolation.

#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include "glog/logging.h"

using ceres::Grid1D;
using ceres::CubicInterpolator;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

// A simple cost functor that interfaces an interpolated table of
// values with automatic differentiation.
struct InterpolatedCostFunctor {
  explicit InterpolatedCostFunctor(
      const CubicInterpolator<Grid1D<double> >& interpolator)
      : interpolator_(interpolator) {
  }

  template<typename T> bool operator()(const T* x, T* residuals) const {
    interpolator_.Evaluate(*x, residuals);
    return true;
  }

  static CostFunction* Create(
      const CubicInterpolator<Grid1D<double> >& interpolator) {
    return new AutoDiffCostFunction<InterpolatedCostFunctor, 1, 1>(
        new InterpolatedCostFunctor(interpolator));
  }

 private:
  const CubicInterpolator<Grid1D<double> >& interpolator_;
};

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  // Evaluate the function f(x) = (x - 4.5)^2;
  const int kNumSamples = 10;
  double values[kNumSamples];
  for (int i = 0; i < kNumSamples; ++i) {
    values[i] = (i - 4.5) * (i - 4.5);
  }

  Grid1D<double> array(values, 0, kNumSamples);
  CubicInterpolator<Grid1D<double> > interpolator(array);

  double x = 1.0;
  Problem problem;
  CostFunction* cost_function = InterpolatedCostFunctor::Create(interpolator);
  problem.AddResidualBlock(cost_function, NULL, &x);

  Solver::Options options;
  options.minimizer_progress_to_stdout = true;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";
  std::cout << "Expected x: 4.5. Actual x : " << x << std::endl;
  return 0;
}
