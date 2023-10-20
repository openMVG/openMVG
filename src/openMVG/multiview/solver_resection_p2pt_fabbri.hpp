// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Ricardo Fabbri, Ariel Kovaljski
//
// Author: Ariel Kovaljski and Ricardo Fabbri
// Rio de Janeiro State University

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPENMVG_MULTIVIEW_RESECTION_P2Pt_FABBRI_HPP
#define OPENMVG_MULTIVIEW_RESECTION_P2Pt_FABBRI_HPP

#include "openMVG/multiview/two_view_kernel.hpp"
#include "openMVG/multiview/solver_resection_metrics.hpp"

namespace openMVG
{
namespace euclidean_resection
{

struct P2PtSolver_Fabbri
{
  enum { MINIMUM_SAMPLES = 2 };
  enum { MAX_MODELS = 8 };

  // Solve for absolute camera pose from 2 points with (SIFT) orientation
  // Paper: Camera Pose Estimation Using Curve Differential Geometry, 
  // IEEE Transactions on Pattern Analysis and Machine Intelligence, 2020, 
  // Ricardo Fabbri, Peter J. Giblin and Benjamin Kimia
  static void Solve(
      const Mat &bearing_vectors,
      const Mat &tangent_vectors,
      const Mat &X, // 3D points
      const Mat &T, // 3D tangents
      std::vector<Mat34> *models);

  static void Solve(
      const Mat &info_2d, // 2D points and tangents x y z tx ty tz where z=1 tz  =0
      const Mat &info_3d, // 3D points and tangents
      std::vector<Mat34> *models)
  {
    constexpr double eps = 1e-8;
    assert(info_2d.cols() == info_3d.cols());
    assert(info_2d.cols() == 2);
    assert(info_2d.rows() == 6); 
    assert(info_3d.rows() == 6);
    assert(fabs(info_2d(5,1)) < eps  && fabs(info_2d(5,2)) < eps);
    assert(fabs(info_2d(2,1) - 1) < eps && fabs(info_2d(2,2) - 1) < eps);

    // TODO h-normalize vectors to comput3 angular error

    Solve(info_2d.block<0,0>(2,2), info_2d.block<3,0>(2,2),
          info_3d.block<0,0>(3,2), info_3d.block<3,0>(3,2), models);
  }
};

//-- Usable solver for robust estimation framework
using PoseResectionKernel_P2Pt_Fabbri =
    two_view::kernel::Kernel<
       P2PtSolver_Fabbri, // Model estimator
        resection::AngularReprojectionError, // Error metric
        Mat34>;

} // namespace euclidean_resection
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P2Pt_FABBRI_HPP
