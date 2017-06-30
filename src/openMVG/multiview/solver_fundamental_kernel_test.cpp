
// Copyright (c) 2010 libmv authors.
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

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/numeric/numeric.h"

#include "testing/testing.h"

#include <numeric>

using namespace openMVG;
using namespace std;

// Check that sin(angle(a, b)) < tolerance.
template<typename A, typename B>
bool Colinear(const A &a, const B &b, double tolerance) {
  const bool dims_match = (a.rows() == b.rows()) && (a.cols() == b.cols());
  if (!dims_match) {
    return false;
  }
  const double c = CosinusBetweenMatrices(a, b);
  if (c * c < 1) {
    const double s = sqrt(1 - c * c);
    return std::abs(s) < tolerance;
  }
  return false;
}

// Check the properties of a fundamental matrix:
//
//   1. The determinant is 0 (rank deficient)
//   2. The condition x'T*F*x = 0 is satisfied to precision.
//
bool ExpectFundamentalProperties(const Mat3 &F,
                                 const Mat &ptsA,
                                 const Mat &ptsB,
                                 double precision) {
  bool bOk = true;
  bOk &= F.determinant() < precision;
  assert(ptsA.cols() == ptsB.cols());
  Mat hptsA, hptsB;
  EuclideanToHomogeneous(ptsA, &hptsA);
  EuclideanToHomogeneous(ptsB, &hptsB);
  for (int i = 0; i < ptsA.cols(); ++i) {
    const double residual = hptsB.col(i).dot(F * hptsA.col(i));
    bOk &= residual < precision;
  }
  return bOk;
}

// Check the fundamental fitting:
//
//   1. Estimate the fundamental matrix
//   2. Check epipolar distance.
//
template <class Kernel>
bool ExpectKernelProperties(const Mat &x1,
                            const Mat &x2,
                            Mat3 *F_expected = nullptr) {
  bool bOk = true;
  Kernel kernel(x1, x2);
  vector<uint32_t> samples(x1.cols());
  std::iota(samples.begin(), samples.end(), 0);
  vector<Mat3> Fs;
  kernel.Fit(samples, &Fs);

  bOk &= (!Fs.empty());
  for (int i = 0; i < Fs.size(); ++i) {
    bOk &= ExpectFundamentalProperties(Fs[i], x1, x2, 1e-8);
    if (F_expected) {
      bOk &= Colinear(Fs[i], *F_expected, 1e-6);
    }
  }
  return bOk;
}

TEST(SevenPointTest, EasyCase) {
  Mat x1(2, 7), x2(2, 7);
  x1 << 0, 0, 0, 1, 1, 1, 2,
        0, 1, 2, 0, 1, 2, 0;
  x2 << 0, 0, 0, 1, 1, 1, 2,
        1, 2, 3, 1, 2, 3, 1;
  using Kernel = fundamental::kernel::SevenPointKernel;
  EXPECT_TRUE(ExpectKernelProperties<Kernel>(x1, x2));
}

TEST(SevenPointTest_Normalized, EasyCase) {
  Mat x1(2, 7), x2(2, 7);
  x1 << 0, 0, 0, 1, 1, 1, 2,
    0, 1, 2, 0, 1, 2, 0;
  x2 << 0, 0, 0, 1, 1, 1, 2,
    1, 2, 3, 1, 2, 3, 1;
  using Kernel = fundamental::kernel::NormalizedSevenPointKernel;
  EXPECT_TRUE(ExpectKernelProperties<Kernel>(x1, x2));
}

TEST(EightPointTest, EasyCase) {
  Mat x1(2, 8), x2(2, 8);
  x1 << 0, 0, 0, 1, 1, 1, 2, 2,
        0, 1, 2, 0, 1, 2, 0, 1;
  x2 << 0, 0, 0, 1, 1, 1, 2, 2,
        1, 2, 3, 1, 2, 3, 1, 2;
  using Kernel = fundamental::kernel::EightPointKernel;
  EXPECT_TRUE(ExpectKernelProperties<Kernel>(x1, x2));
}

TEST(EightPointTest_Normalized, EasyCase) {
  Mat x1(2, 8), x2(2, 8);
  x1 << 0, 0, 0, 1, 1, 1, 2, 2,
    0, 1, 2, 0, 1, 2, 0, 1;
  x2 << 0, 0, 0, 1, 1, 1, 2, 2,
    1, 2, 3, 1, 2, 3, 1, 2;
  using Kernel = fundamental::kernel::NormalizedEightPointKernel;
  EXPECT_TRUE(ExpectKernelProperties<Kernel>(x1, x2));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
