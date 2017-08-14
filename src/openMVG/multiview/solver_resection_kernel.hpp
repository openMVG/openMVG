
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

#ifndef OPENMVG_MULTIVIEW_RESECTION_KERNEL_HPP
#define OPENMVG_MULTIVIEW_RESECTION_KERNEL_HPP

#include <vector>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG {
namespace resection {
namespace kernel {

/**
 * Six-point resection
 * P Matrix estimation (Pose estimation)
 * Compute a projection matrix using linear least squares.
 * Rely on Linear Resection algorithm.
 * Work from 6 to N points.
 */
struct SixPointResectionSolver {
  enum { MINIMUM_SAMPLES = 6 };
  enum { MAX_MODELS = 1 };
  // Solve the problem of camera pose.
  // First 3d point will be translated in order to have X0 = (0,0,0,1)
  static void Solve(const Mat &pt2D, const Mat &pt3D, std::vector<Mat34> *P, bool bcheck = true);

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D){
    return (pt2D - Project(P, pt3D)).norm();
  }
};

//-- Generic Solver for the 6pt Resection algorithm using linear least squares.
template<typename SolverArg,
  typename ErrorArg,
  typename ModelArg = Mat34>
class ResectionKernel :
   public two_view::kernel::Kernel<SolverArg,ErrorArg, ModelArg>
{
public:
  // 2D / 3D points
  ResectionKernel(const Mat &pt2D, const Mat &pt3D):
  two_view::kernel::Kernel<SolverArg,ErrorArg, ModelArg>(pt2D,pt3D){}

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<ModelArg> *models
  )
  const
  {
    const Mat pt2d = ExtractColumns(this->x1_, samples);
    const Mat pt3d = ExtractColumns(this->x2_, samples);

    assert(2 == pt2d.rows());
    assert(3 == pt3d.rows());
    assert(SolverArg::MINIMUM_SAMPLES <= pt2d.cols());
    assert(pt2d.cols() == pt3d.cols());

    SolverArg::Solve(pt2d, pt3d, models);
  }

  // Error : re-projection error of the sample
  double Error
  (
    size_t sample,
    const ModelArg &model
  )
  const
  {
    return ErrorArg::Error(model, this->x1_.col(sample), this->x2_.col(sample));
  }
};

//-- Usable solver for the 6pt Resection estimation
using PoseResectionKernel =
  two_view::kernel::Kernel<
    SixPointResectionSolver, // Model estimator
    SixPointResectionSolver, // Error metric
    Mat34>;

}  // namespace kernel
}  // namespace resection
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_RESECTION_KERNEL_HPP
