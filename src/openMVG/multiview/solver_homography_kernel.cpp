
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


// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_homography_kernel.hpp"

namespace openMVG {
namespace homography {
namespace kernel {

/// Setup the Direct Linear Transform.
///  Use template in order to support fixed or dynamic sized matrix.
/// Allow solve H as homogeneous(x2) = H homogeneous(x1)
template<typename Matrix >
void BuildActionMatrix(Matrix & L, const Mat &x, const Mat &y)  {

  const Mat::Index n = x.cols();
  for (Mat::Index i = 0; i < n; ++i) {
    Mat::Index j = 2 * i;
    L(j, 0) = x(0, i);
    L(j, 1) = x(1, i);
    L(j, 2) = 1.0;
    L(j, 6) = -y(0, i) * x(0, i);
    L(j, 7) = -y(0, i) * x(1, i);
    L(j, 8) = -y(0, i);

    ++j;
    L(j, 3) = x(0, i);
    L(j, 4) = x(1, i);
    L(j, 5) = 1.0;
    L(j, 6) = -y(1, i) * x(0, i);
    L(j, 7) = -y(1, i) * x(1, i);
    L(j, 8) = -y(1, i);
  }
}

void FourPointSolver::Solve(const Mat &x, const Mat &y, vector<Mat3> *Hs) {
  assert(2 == x.rows());
  assert(4 <= x.cols());
  assert(x.rows() == y.rows());
  assert(x.cols() == y.cols());

  Mat::Index n = x.cols();

  Vec9 h;
  if (n == 4)  {
    // In the case of minimal configuration we use fixed sized matrix to let
    //  Eigen and the compiler doing the maximum of optimization.
    typedef Eigen::Matrix<double, 16, 9> Mat16_9;
    Mat16_9 L = Mat::Zero(16, 9);
    BuildActionMatrix(L, x, y);
    Nullspace(&L, &h);
  }
  else {
    MatX9 L = Mat::Zero(n * 2, 9);
    BuildActionMatrix(L, x, y);
    Nullspace(&L, &h);
  }
  Mat3 H = Map<RMat3>(h.data()); // map the linear vector as the H matrix
  Hs->push_back(H);
}

}  // namespace kernel
}  // namespace homography
}  // namespace openMVG
