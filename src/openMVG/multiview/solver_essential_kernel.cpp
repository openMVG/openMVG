
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

#include "openMVG/multiview/solver_essential_three_point.hpp"
#include "openMVG/multiview/solver_essential_five_point.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/numeric/nullspace.hpp"
#include "openMVG/numeric/poly.h"

#include <cassert>

namespace openMVG {
namespace essential {
namespace kernel {

void FivePointSolver::Solve(const Mat3X &x1, const Mat3X &x2, std::vector<Mat3> *E) {
  assert(5 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  FivePointsRelativePose(x1, x2, E);
}

void ThreePointSolver::Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *E) {
  assert(2 == x1.rows());
  assert(3 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  ThreePointsRelativePose(x1, x2, E);
}

}  // namespace kernel
}  // namespace essential
}  // namespace openMVG
