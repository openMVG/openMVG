
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

#include "openMVG/multiview/solver_homography_kernel.hpp"

#include "testing/testing.h"

#include <vector>

using namespace std;
using namespace openMVG;

TEST(HomographyKernelTest, Fitting_Unnormalized) {
  // Define 3 knows homographies (Use as GT).
  vector<Mat3> H_gt(3);

  H_gt[0] = Mat3::Identity();
  H_gt[1] << 1,  0, -4, // Affine homography motion
             0,  1,  5,
             0,  0,  1;
  H_gt[2] << 1, -2,  3,
             4,  5, -6,
            -7,  8,  1;

  // Define a set of points.
  Mat x(2, 9), xh;
  x << 0, 0, 0, 1, 1, 1, 2, 2, 2,
       0, 1, 2, 0, 1, 2, 0, 1, 2;
  EuclideanToHomogeneous(x, &xh);

  for (size_t i = 0; i < H_gt.size(); ++i) {
    // Transform points by the ground truth homography.
    Mat y, yh = H_gt[i] * xh;
    HomogeneousToEuclidean(yh, &y);

    homography::kernel::UnnormalizedKernel kernel(x, y);

    vector<uint32_t> samples = {0,1,2,3,4};
    for (
      Mat::Index j = 4;
      static_cast<Mat::Index>(samples.size()) < x.cols();
      samples.push_back(j++))
    {
      vector<Mat3> Hs;
      kernel.Fit(samples, &Hs);
      CHECK_EQUAL(1, Hs.size());
      // Check that found matrix is equal to the GT
      EXPECT_MATRIX_PROP(H_gt[i], Hs[0], 1e-6);
    }
  }
}

TEST(HomographyKernelTest, Fitting_Normalized) {
  // Define 3 knows homographies (Use as GT).
  vector<Mat3> H_gt(3);

  H_gt[0] = Mat3::Identity();
  H_gt[1] << 1,  0, -4, // Affine homography motion
             0,  1,  5,
             0,  0,  1;
  H_gt[2] << 1, -2,  3,
             4,  5, -6,
            -7,  8,  1;

  // Define a set of points.
  Mat x(2, 9), xh;
  x << 0, 0, 0, 1, 1, 1, 2, 2, 2,
       0, 1, 2, 0, 1, 2, 0, 1, 2;
  EuclideanToHomogeneous(x, &xh);

  for (size_t i = 0; i < H_gt.size(); ++i) {
    // Transform points by the ground truth homography.
    Mat y, yh = H_gt[i] * xh;
    HomogeneousToEuclidean(yh, &y);

    homography::kernel::Kernel kernel(x, y);

    vector<uint32_t> samples = {0,1,2,3,4};
    for (
      Mat::Index j = 4;
      static_cast<Mat::Index>(samples.size()) < x.cols();
      samples.push_back(j++))
    {
      vector<Mat3> Hs;
      kernel.Fit(samples, &Hs);
      CHECK_EQUAL(1, Hs.size());
      // Check that found matrix is equal to the GT
      EXPECT_MATRIX_PROP(H_gt[i], Hs[0], 1e-6);
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
