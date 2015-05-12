
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

#ifndef OPENMVG_MULTIVIEW_SOLVER_FUNDAMENTAL_KERNEL_H_
#define OPENMVG_MULTIVIEW_SOLVER_FUNDAMENTAL_KERNEL_H_

#include <vector>
#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace fundamental {
namespace kernel {

using namespace std;

/**
 * Seven-point algorithm for solving for the fundamental matrix from point
 * correspondences. See page 281 in HZ, though oddly they use a different
 * equation: \f$det(\alpha F_1 + (1-\alpha)F_2) = 0\f$. Since \f$F_1\f$ and
 * \f$F2\f$ are projective, there's no need to balance the relative scale.
 * Instead, here, the simpler equation is solved: \f$det(F_1 + \alpha F_2) =
 * 0\f$.
 *
 * \see http://www.cs.unc.edu/~marc/tutorial/node55.html
 */
struct SevenPointSolver {
  enum { MINIMUM_SAMPLES = 7 };
  enum { MAX_MODELS = 3 };
  static void Solve(const Mat &x1, const Mat &x2, vector<Mat3> *F);
};

struct EightPointSolver {
  enum { MINIMUM_SAMPLES = 8 };
  enum { MAX_MODELS = 1 };
  static void Solve(const Mat &x1, const Mat &x2, vector<Mat3> *Fs);
};

/**
 * Build a 9 x n matrix from point matches, where each row is equivalent to the
 * equation x'T*F*x = 0 for a single correspondence pair (x', x). The domain of
 * the matrix is a 9 element vector corresponding to F. In other words, set up
 * the linear system
 *
 *   Af = 0,
 *
 * where f is the F matrix as a 9-vector rather than a 3x3 matrix (row
 * major). If the points are well conditioned and there are 8 or more, then
 * the nullspace should be rank one. If the nullspace is two dimensional,
 * then the rank 2 constraint must be enforced to identify the appropriate F
 * matrix.
 *
 * Note that this does not resize the matrix A; it is expected to have the
 * appropriate size already.
 */
template<typename TMatX, typename TMatA>
inline void EncodeEpipolarEquation(const TMatX &x1, const TMatX &x2, TMatA *A) {
  for (typename TMatX::Index i = 0; i < x1.cols(); ++i) {
    const Vec2 xx1 = x1.col(i);
    const Vec2 xx2 = x2.col(i);
    A->row(i) <<
      xx2(0) * xx1(0),  // 0 represents x coords,
      xx2(0) * xx1(1),  // 1 represents y coords.
      xx2(0),
      xx2(1) * xx1(0),
      xx2(1) * xx1(1),
      xx2(1),
      xx1(0),
      xx1(1),
      1.0;
  }
}

/// Compute SampsonError related to the Fundamental matrix and 2 correspondences
struct SampsonError {
  static double Error(const Mat3 &F, const Vec2 &x1, const Vec2 &x2) {
    Vec3 x(x1(0), x1(1), 1.0);
    Vec3 y(x2(0), x2(1), 1.0);
    // See page 287 equation (11.9) of HZ.
    Vec3 F_x = F * x;
    Vec3 Ft_y = F.transpose() * y;
    return Square(y.dot(F_x)) / (  F_x.head<2>().squaredNorm()
                                + Ft_y.head<2>().squaredNorm());
  }
};

struct SymmetricEpipolarDistanceError {
  static double Error(const Mat3 &F, const Vec2 &x1, const Vec2 &x2) {
    Vec3 x(x1(0), x1(1), 1.0);
    Vec3 y(x2(0), x2(1), 1.0);
    // See page 288 equation (11.10) of HZ.
    Vec3 F_x = F * x;
    Vec3 Ft_y = F.transpose() * y;
    return Square(y.dot(F_x)) * ( 1.0 / F_x.head<2>().squaredNorm()
                                + 1.0 / Ft_y.head<2>().squaredNorm())
      / 4.0;  // The divide by 4 is to make this match the Sampson distance.
  }
};

struct EpipolarDistanceError {
  static double Error(const Mat3 &F, const Vec2 &x1, const Vec2 &x2) {
    // Transfer error in image 2
    // See page 287 equation (11.9) of HZ.
    Vec3 x(x1(0), x1(1), 1.0);
    Vec3 y(x2(0), x2(1), 1.0);
    Vec3 F_x = F * x;
    return Square(F_x.dot(y)) /  F_x.head<2>().squaredNorm();
  }
};
typedef EpipolarDistanceError SimpleError;

//-- Kernel solver for the 8pt Fundamental Matrix Estimation
typedef two_view::kernel::Kernel<SevenPointSolver, SampsonError, Mat3>
  SevenPointKernel;

//-- Kernel solver for the 8pt Fundamental Matrix Estimation
typedef two_view::kernel::Kernel<EightPointSolver, SampsonError, Mat3>
  EightPointKernel;

//-- Normalized 7pt kernel -> conditioning from HZ (Algo 11.1) pag 282
typedef two_view::kernel::Kernel<
  two_view::kernel::NormalizedSolver<SevenPointSolver, UnnormalizerT>,
  SampsonError,
  Mat3>
  NormalizedSevenPointKernel;

//-- Normalized 8pt kernel -> conditioning from HZ (Algo 11.1) pag 282
typedef two_view::kernel::Kernel<
  two_view::kernel::NormalizedSolver<EightPointSolver, UnnormalizerT>,
  SampsonError,
  Mat3>
  NormalizedEightPointKernel;

}  // namespace kernel
}  // namespace fundamental
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_FUNDAMENTAL_KERNEL_H_
