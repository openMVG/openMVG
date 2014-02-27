// Copyright (c) 2014 Matthieu Tourne

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GOLD_STANDARD_AFFINE_H_
#define OPENMVG_GOLD_STANDARD_AFFINE_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/two_view_kernel.hpp"

namespace openMVG {
namespace fund_affine {
namespace kernel {

struct GSAffineSolver {
    enum { MINIMUM_SAMPLES = 4 };
    enum { MAX_MODELS = 1 };

    static void Solve(const Mat2X &x1, const Mat2X &x2, std::vector<Mat3> *Fa);
};

struct MinimalAffineSolver {
  enum { MINIMUM_SAMPLES = 4 };
  enum { MAX_MODELS = 1 };

  static void Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *Fa);
};

struct AffineSampsonError {
    // Sampson error
    // HZ 14.4
    // Shape of the Fa Matrix
    //  0 0 a
    //  0 0 b
    //  c d e
    static double Error(const Mat3 &Fa, const Vec2& x1, const Vec2& x2) {
        // stack vertically
        Vec4 X;
        X <<    x2,
                x1;

        // a b c d
        Vec4 affin_lin(Fa(0,2), Fa(1,2), Fa(2,0), Fa(2,1));
        double e = Fa(2, 2);

        // (ax2 + by2 + cx1 + dy1 + e) / (a2 + b2 + c2 + d2)
        double d = ((affin_lin.cwiseProduct(X)).sum() + e) / affin_lin.squaredNorm();
        Vec4 H = X - d * affin_lin;

        return (X - H).squaredNorm();
    }
};

typedef two_view::kernel::Kernel<MinimalAffineSolver, AffineSampsonError, Mat3>
  MinimalAffineKernel;

typedef two_view::kernel::Kernel<GSAffineSolver, AffineSampsonError, Mat3>
  GSAffineKernel;

} // namespace kernel
} // namespace fund_affine
} // namespace openMVG

#endif // OPENMVG_GOLD_STANDARD_AFFINE_H_
