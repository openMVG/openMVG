
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

#include "openMVG/multiview/solver_affine.hpp"

namespace openMVG {

// Parametrization
// a b x
// c d y
// 0 0 1

// It gives the following system A x = B :
// | X1 Y1 0  0  1 0 |  | a |   | X2 |
// | 0  0  X1 Y1 0 1 |  | b | = | Y2 |
//         ...          | c |     ..
//                      | d |
//                      | x |
//                      | y |
bool Affine2DFromCorrespondencesLinear(const Mat &x1, const Mat &x2,
                                       Mat3 *M,
                                       double expected_precision) {
  assert(2 == x1.rows());
  assert(3 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  const Mat::Index n = x1.cols();
  Mat A = Mat::Zero(2*n, 6);
  Mat b = Mat::Zero(2*n, 1);
  for (Mat::Index i = 0; i < n; ++i) {
    const Mat::Index j= i * 2;
    A(j,0) =  x1(0,i);
    A(j,1) =  x1(1,i);
    A(j,4) =  1.0;

    A(j+1,2) = x1(0,i);
    A(j+1,3) = x1(1,i);
    A(j+1,5) = 1.0;

    b(j,0)   = x2(0,i);
    b(j+1,0) = x2(1,i);
  }
  // Solve A x = B
  const Vec x = A.fullPivLu().solve(b);
  if ((A * x).isApprox(b, expected_precision))  {
    *M << x(0), x(1), x(4),
          x(2), x(3), x(5),
          0.0,  0.0,  1.0;
    return true;
  }
  return false;
}

// Parametrization
// a b c x
// d e f y
// g h i z
// 0 0 0 1

// It gives the following system A x = B :
// | X1 Y1 Z1 0  0  0  0  0  0  1 0 0|  | a |   | X2 |
// | 0  0  0  X1 Y1 Z1 0  0  0  0 1 0|  | b | = | Y2 |
// | 0  0  0  0  0  0  X1 Y1 Z1 0 0 1|  | c |   | Z2 |
//                   ...                | d |     ..
//                                      | e |
//                                      | f |
//                                      | g |
//                                      | h |
//                                      | i |
//                                      | x |
//                                      | y |
//                                      | z |
bool Affine3DFromCorrespondencesLinear(const Mat &x1,
                                       const Mat &x2,
                                       Mat4 *M,
                                       double expected_precision) {
  assert(3 == x1.rows());
  assert(4 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  const Mat::Index n = x1.cols();
  Mat A = Mat::Zero(3*n, 12);
  Mat b = Mat::Zero(3*n, 1);
  for (Mat::Index i = 0; i < n; ++i) {
    const Mat::Index j= i * 3;
    const Mat::Index j1= j + 1;
    const Mat::Index j2= j + 2;
    A(j,0) =  x1(0,i);
    A(j,1) =  x1(1,i);
    A(j,2) =  x1(2,i);
    A(j,9) =  1.0;

    A(j1,3) =  x1(0,i);
    A(j1,4) =  x1(1,i);
    A(j1,5) =  x1(2,i);
    A(j1,10)=  1.0;

    A(j2,6) =  x1(0,i);
    A(j2,7) =  x1(1,i);
    A(j2,8) =  x1(2,i);
    A(j2,11)=  1.0;

    b(j,0)  = x2(0,i);
    b(j1,0) = x2(1,i);
    b(j2,0) = x2(2,i);
  }
  // Solve A x = B
  const Vec x = A.fullPivLu().solve(b);
  if ((A * x).isApprox(b, expected_precision))  {
    *M << x(0), x(1), x(2), x(9), // a b c x
          x(3), x(4), x(5), x(10), // d e f y
          x(6), x(7), x(8), x(11), // g h i z
          0.0,   0.0,  0.0, 1.0;
    return true;
  }
  return false;
}
} // namespace openMVG
