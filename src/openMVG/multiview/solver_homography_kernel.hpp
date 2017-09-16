
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

#ifndef OPENMVG_MULTIVIEW_SOLVER_HOMOGRAPHY_KERNEL_HPP
#define OPENMVG_MULTIVIEW_SOLVER_HOMOGRAPHY_KERNEL_HPP

#include <vector>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/two_view_kernel.hpp"

namespace openMVG {
namespace homography {
namespace kernel {

struct FourPointSolver {
  enum { MINIMUM_SAMPLES = 4 };
  enum { MAX_MODELS = 1 };
  /**
   * Computes the homography that transforms x to y with the Direct Linear
   * Transform (DLT).
   *
   * \param x  A 2xN matrix of column vectors.
   * \param y  A 2xN matrix of column vectors.
   * \param Hs A vector into which the computed homography is stored.
   *
   * The estimated homography should approximately hold the condition y = H x.
   */
  static void Solve(const Mat &x, const Mat &y, std::vector<Mat3> *Hs);
};

// Should be distributed as Chi-squared with k = 2.
struct AsymmetricError {
  static double Error(const Mat &H, const Vec2 &x, const Vec2 &y) {
    return (y - Vec3( H * x.homogeneous()).hnormalized() ).squaredNorm();
  }
};

// Kernel that works on original data point
using UnnormalizedKernel =
  two_view::kernel::Kernel<
    FourPointSolver,
    AsymmetricError,
    Mat3>;

// By default use the normalized version for increased robustness.
using Kernel =
  two_view::kernel::Kernel<
    two_view::kernel::NormalizedSolver<FourPointSolver, UnnormalizerI>,
    AsymmetricError,
    Mat3>;

}  // namespace kernel
}  // namespace homography
}  // namespace openMVG

#endif // OPENMVG_MULTIVIEW_SOLVER_HOMOGRAPHY_KERNEL_HPP
