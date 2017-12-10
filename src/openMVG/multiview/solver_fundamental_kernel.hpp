
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

#ifndef OPENMVG_MULTIVIEW_SOLVER_FUNDAMENTAL_KERNEL_HPP
#define OPENMVG_MULTIVIEW_SOLVER_FUNDAMENTAL_KERNEL_HPP

#include <vector>

#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace fundamental {
namespace kernel {


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
  static void Solve(const Mat2X &x1, const Mat2X &x2, std::vector<Mat3> *F);
};

struct EightPointSolver {
  enum { MINIMUM_SAMPLES = 8 };
  enum { MAX_MODELS = 1 };
  static void Solve(const Mat2X &x1, const Mat2X &x2, std::vector<Mat3> *Fs);
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
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());
  assert(x1.rows() == 3);
  for (typename TMatX::Index i = 0; i < x1.cols(); ++i) {
    const auto xx1 = x1.col(i).homogeneous().transpose();
    const auto xx2 = x2.col(i).homogeneous().transpose();
    A->row(i) <<
      x2(0, i) * x1.col(i).transpose(),
      x2(1, i) * x1.col(i).transpose(),
      x2(2, i) * x1.col(i).transpose();
  }
}

/// Compute SampsonError related to the Fundamental matrix and 2 correspondences
struct SampsonError {
  static double Error(const Mat3 &F, const Vec2 &x, const Vec2 &y);
};

struct SymmetricEpipolarDistanceError {
  static double Error(const Mat3 &F, const Vec2 &x, const Vec2 &y);
};

struct EpipolarDistanceError {
  static double Error(const Mat3 &F, const Vec2 &x, const Vec2 &y);
};

//-- Kernel solver for the 8pt Fundamental Matrix Estimation
using SevenPointKernel = two_view::kernel::Kernel<SevenPointSolver, SampsonError, Mat3>;

//-- Kernel solver for the 8pt Fundamental Matrix Estimation
using EightPointKernel = two_view::kernel::Kernel<EightPointSolver, SampsonError, Mat3>;

//-- Normalized 7pt kernel -> conditioning from HZ (Algo 11.1) pag 282
using NormalizedSevenPointKernel =
  two_view::kernel::Kernel<
    two_view::kernel::NormalizedSolver<SevenPointSolver, UnnormalizerT>,
    SampsonError,
    Mat3>;

//-- Normalized 8pt kernel -> conditioning from HZ (Algo 11.1) pag 282
using NormalizedEightPointKernel =
  two_view::kernel::Kernel<
    two_view::kernel::NormalizedSolver<EightPointSolver, UnnormalizerT>,
    SampsonError,
    Mat3>;

}  // namespace kernel
}  // namespace fundamental
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_FUNDAMENTAL_KERNEL_HPP
