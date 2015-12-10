
/*
 * Copyright (c) 2011, Laurent Kneip, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ETH Zurich nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_P3P_H_
#define OPENMVG_MULTIVIEW_RESECTION_P3P_H_

#include <iostream>
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace euclidean_resection {

typedef Eigen::Matrix<double, 5, 1> Vec5;

static void solveQuartic(const Vec5 & factors, Vec4 & realRoots);

/*
 *      Author: Laurent Kneip, adapted to openMVG by Pierre Moulon
 * Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
 *   Reference: [1] A Novel Parametrization of the P3P-Problem for a Direct Computation of
 *              Absolute Camera Position and Orientation
 *              Kneip, L.; Scaramuzza, D. ; Siegwart, R.
 *              CVPR 2011
 *
 *       Input: featureVectors: 3x3 matrix with UNITARY feature vectors (each column is a vector)
 *              worldPoints: 3x3 matrix with corresponding 3D world points (each column is a point)
 *              solutions: 3x16 matrix that will contain the solutions
 *                         form: [ C1,R1, C2,R2 ... ]
 *                         the obtained orientation matrices are defined as transforming points from the cam to the world frame
 *      Output: bool: true if correct execution
 *                    false if world points aligned
 */

static bool compute_P3P_Poses(const Mat3 & featureVectors, const Mat3 & worldPoints, Mat & solutions);

struct P3PSolver
{

  enum
  {
    MINIMUM_SAMPLES = 3
  };

  enum
  {
    MAX_MODELS = 4
  };
  // Solve the problem of camera pose.
  static void Solve(const Mat &pt2D, const Mat &pt3D, std::vector<Mat34> *models);

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D);
};

class P3P_ResectionKernel_K
{
public:
  typedef Mat34 Model;

  enum
  {
    MINIMUM_SAMPLES = 3
  };

  P3P_ResectionKernel_K(const Mat2X &x_camera, const Mat3X &X, const Mat3 &K = Mat3::Identity());

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const;

  double Error(size_t sample, const Model &model) const;

  size_t NumSamples() const;

private:
  Mat2X x_image_; // camera coordinates
  Mat3X x_camera_; // camera coordinates (normalized)
  Mat3X X_; // 3D points
  Mat3 K_;
};

} // namespace euclidean_resection
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P3P_H_

