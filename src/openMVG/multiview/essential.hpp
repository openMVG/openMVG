
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

#ifndef OPENMVG_MULTIVIEW_ESSENTIAL_H_
#define OPENMVG_MULTIVIEW_ESSENTIAL_H_

#include <vector>

#include "openMVG/numeric/numeric.h"

namespace openMVG {

// Compute the relative camera motion between two cameras.
// Given the motion parameters of two cameras, computes the motion parameters
// of the second one assuming the first one to be at the origin.
// If T1 and T2 are the camera motions, the computed relative motion is
//
//      T = T2 T1^{-1}
//
void RelativeCameraMotion(const Mat3 &R1,
                          const Vec3 &t1,
                          const Mat3 &R2,
                          const Vec3 &t2,
                          Mat3 *R,
                          Vec3 *t);

/// Given F, Left/Right K matrix it compute the Essential matrix
void EssentialFromFundamental(const Mat3 &F,
                              const Mat3 &K1,
                              const Mat3 &K2,
                              Mat3 *E);

/// Compute E as E = [t12]x R12.
void EssentialFromRt(const Mat3 &R1,
                     const Vec3 &t1,
                     const Mat3 &R2,
                     const Vec3 &t2,
                     Mat3 *E);

/// Given E, Left/Right K matrix it compute the Fundamental matrix
void FundamentalFromEssential(const Mat3 &E,
                              const Mat3 &K1,
                              const Mat3 &K2,
                              Mat3 *F);

/// Test the possible R|t configuration to have point in front of the cameras
/// Return false if no possible configuration
bool MotionFromEssentialAndCorrespondence(const Mat3 &E,
                                          const Mat3 &K1,
                                          const Vec2 &x1,
                                          const Mat3 &K2,
                                          const Vec2 &x2,
                                          Mat3 *R,
                                          Vec3 *t);

/// Choose one of the four possible motion solutions from an essential matrix.
/// Decides the right solution by checking that the triangulation of a match
/// x1--x2 lies in front of the cameras.
/// Return the index of the right solution or -1 if no solution.
int MotionFromEssentialChooseSolution(const std::vector<Mat3> &Rs,
                                      const std::vector<Vec3> &ts,
                                      const Mat3 &K1,
                                      const Vec2 &x1,
                                      const Mat3 &K2,
                                      const Vec2 &x2);

// HZ 9.7 page 259 (Result 9.19)
void MotionFromEssential(const Mat3 &E,
  std::vector<Mat3> *Rs,
  std::vector<Vec3> *ts);


} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_ESSENTIAL_H_
