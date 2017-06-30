
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

#ifndef OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_KERNEL_HPP
#define OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_KERNEL_HPP

#include <vector>

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/multiview/two_view_kernel.hpp"

namespace openMVG {
namespace essential {
namespace kernel {

/**
 * Eight-point algorithm for solving for the essential matrix from normalized
 * image coordinates of point correspondences.
 * See page 294 in HZ Result 11.1.
 *
 */
struct EightPointRelativePoseSolver {
  enum { MINIMUM_SAMPLES = 8 };
  enum { MAX_MODELS = 1 };
  static void Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *E);
};

/**
 * Five-point algorithm to solve the Essential matrix from 5 points
 * correspondences. It solves the relative pose problem.
 * Input point must be normalized one.
 */
struct FivePointSolver {
  enum { MINIMUM_SAMPLES = 5 };
  enum { MAX_MODELS = 10 };
  static void Solve(const Mat &x1, const Mat &x2, std::vector<Mat3> *E);
};

//-- Generic Solver for the 5pt Essential Matrix Estimation.
//-- Need a new Class that inherit of two_view::kernel::kernel.
//    Error must be overwrite in order to compute F from E and K's.
//-- Fitting must normalize image values to camera values.
template<typename SolverArg,
  typename ErrorArg,
  typename ModelArg = Mat3>
class EssentialKernel :
   public two_view::kernel::Kernel<SolverArg,ErrorArg, ModelArg>
{
public:
  EssentialKernel
  (
    const Mat &x1,
    const Mat &x2,
    const Mat3 &K1,
    const Mat3 &K2
   ):
    two_view::kernel::Kernel<SolverArg,ErrorArg, ModelArg>(x1,x2),
    K1_(K1),
    K2_(K2)
  {}

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<ModelArg> *models
  )
  const
  {
    const Mat x1 = ExtractColumns(this->x1_, samples);
    const Mat x2 = ExtractColumns(this->x2_, samples);

    assert(2 == x1.rows());
    assert(SolverArg::MINIMUM_SAMPLES <= x1.cols());
    assert(x1.rows() == x2.rows());
    assert(x1.cols() == x2.cols());

    // Normalize the data (image coords to camera coords).
    const Mat x1_normalized = (K1_.inverse() * x1.colwise().homogeneous()).colwise().hnormalized();
    const Mat x2_normalized = (K2_.inverse() * x2.colwise().homogeneous()).colwise().hnormalized();
    SolverArg::Solve(x1_normalized, x2_normalized, models);
  }

  double Error
  (
    size_t sample,
    const ModelArg &model
  )
  const
  {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    return ErrorArg::Error(F, this->x1_.col(sample), this->x2_.col(sample));
  }
protected:
  Mat3 K1_, K2_; // The two calibrated camera matrices
};

//-- Solver kernel for the 8pt Essential Matrix Estimation
using EightPointKernel = essential::kernel::EssentialKernel<EightPointRelativePoseSolver,
  fundamental::kernel::SampsonError, Mat3>;


//-- Solver kernel for the 5pt Essential Matrix Estimation
using FivePointKernel = essential::kernel::EssentialKernel<FivePointSolver,
  fundamental::kernel::SampsonError, Mat3>;


}  // namespace kernel
}  // namespace essential
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_ESSENTIAL_KERNEL_HPP
