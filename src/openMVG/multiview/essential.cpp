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

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/numeric.h"

#include <array>

namespace openMVG {

// HZ 9.6 page 257 (formula 9.12)
void EssentialFromFundamental(const Mat3 &F,
                              const Mat3 &K1,
                              const Mat3 &K2,
                              Mat3 *E) {
  *E = K2.transpose() * F * K1;
}

// HZ 9.6 page 257 (formula 9.12)
// Or http://ai.stanford.edu/~birch/projective/node20.html
void FundamentalFromEssential(const Mat3 &E,
                              const Mat3 &K1,
                              const Mat3 &K2,
                              Mat3 *F)  {
  *F = K2.inverse().transpose() * E * K1.inverse();
}

void RelativeCameraMotion(const Mat3 &R1,
                          const Vec3 &t1,
                          const Mat3 &R2,
                          const Vec3 &t2,
                          Mat3 *R,
                          Vec3 *t) {
  *R = R2 * R1.transpose();
  *t = t2 - (*R) * t1;
}

// HZ 9.6 page 257
void EssentialFromRt(const Mat3 &R1,
                     const Vec3 &t1,
                     const Mat3 &R2,
                     const Vec3 &t2,
                     Mat3 *E) {
  Mat3 R;
  Vec3 t;
  RelativeCameraMotion(R1, t1, R2, t2, &R, &t);
  const Mat3 Tx = CrossProductMatrix(t);
  *E = Tx * R;
}

// HZ 9.7 page 259 (Result 9.19)
void MotionFromEssential(const Mat3 &E,
                         std::vector<geometry::Pose3> * relative_poses) {
  Eigen::JacobiSVD<Mat3> USV(E, Eigen::ComputeFullU|Eigen::ComputeFullV);
  Mat3 U =  USV.matrixU();
  Mat3 Vt = USV.matrixV().transpose();

  // Last column of U is undetermined since d = (a a 0).
  if (U.determinant() < 0) {
    U.col(2) *= -1;
  }
  // Last row of Vt is undetermined since d = (a a 0).
  if (Vt.determinant() < 0) {
    Vt.row(2) *= -1;
  }

  Mat3 W;
  W << 0, -1,  0,
       1,  0,  0,
       0,  0,  1;

  const Mat3 U_W_Vt = U * W * Vt;
  const Mat3 U_Wt_Vt = U * W.transpose() * Vt;

  const std::array<Mat3, 2> R{{U_W_Vt, U_Wt_Vt}};
  const std::array<Vec3, 2> t{{U.col(2), -U.col(2)}};
  if (relative_poses)
  {
    relative_poses->clear();
    relative_poses->reserve(4);
    relative_poses->emplace_back(R[0], -R[0].transpose() * t[0]);
    relative_poses->emplace_back(R[1], -R[1].transpose() * t[1]);
    relative_poses->emplace_back(R[0], -R[0].transpose() * t[1]);
    relative_poses->emplace_back(R[1], -R[1].transpose() * t[0]);
  }
}

}  // namespace openMVG
