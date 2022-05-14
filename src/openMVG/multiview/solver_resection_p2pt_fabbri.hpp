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
  enum
  {
    MINIMUM_SAMPLES = 2
  };
  enum
  {
    MAX_MODELS = 8
  };

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
};

} // namespace euclidean_resection
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_RESECTION_P2Pt_FABBRI_HPP
