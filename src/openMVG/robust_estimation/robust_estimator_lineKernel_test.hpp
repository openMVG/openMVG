
// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATION_LINEKERNEL_TEST_H_
#define OPENMVG_ROBUST_ESTIMATION_LINEKERNEL_TEST_H_

#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace robust{

struct LineSolver {
  enum { MINIMUM_SAMPLES = 2 };
  enum { MAX_MODELS = 1 };

  static void Solve(const Mat &x, std::vector<Vec2> *lines)
  {
    Mat X(x.cols(), 2);
    X.col(0).setOnes();
    X.col(1) = x.row(0).transpose();
    Mat A(X.transpose() * X);
    const Vec b(X.transpose() * x.row(1).transpose());
    Eigen::JacobiSVD<Mat> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    lines->push_back(svd.solve(b));
  }
};

struct pointToLineError {
  static double Error(const Vec2 &lineEq, const Vec2 &xs) {
    const double b = lineEq[0];
    const double a = lineEq[1];
    const double x = xs[0];
    const double y = xs[1];
    const double e = y - (a*x + b);
    return e*e;
  }
};

// Embed the basic solver to fit from sampled point set
struct LineKernel {
  typedef Vec2 Model;  // line parametrization: a, b;
  enum { MINIMUM_SAMPLES = 2 };

  LineKernel(const Mat2X &xs) : xs_(xs) {}

  size_t NumSamples() const { return static_cast<size_t> (xs_.cols()); }

  void Fit(const std::vector<size_t> &samples, std::vector<Vec2> *lines) const {
    assert(samples.size() >= (unsigned int)MINIMUM_SAMPLES);
    // Standard least squares solution.
    const Mat2X sampled_xs = ExtractColumns(xs_, samples);

    LineSolver::Solve(sampled_xs, lines);
  }

  double Error(size_t sample, const Vec2 &ba) const {
    return pointToLineError::Error(ba, xs_.col(sample));
  }

  const Mat2X &xs_;
};

} // namespace robust
} // namespace openMVG

#endif // OPENMVG_ROBUST_ESTIMATION_LINEKERNEL_TEST_H_
