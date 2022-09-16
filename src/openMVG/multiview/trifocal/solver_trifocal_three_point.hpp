// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 11:55:58 -03 2021
//\author Pierre MOULON
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP
#define OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP

#include <vector>
#include "openMVG/multiview/trifocal/trifocal_model.hpp"

namespace openMVG {
namespace trifocal {
  
struct Trifocal3PointPositionTangentialSolver {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 312 }; // Usually returns around 40 real solutions

  // Solve for 3 cameras given 3 corresponding SIFT features across 3 views
  // Or 3 points with 3 lines passing through the points, across 3 views
  // Or 3 points on curves (points and curve tangents)
  //
  // Solver actually uses only tangent orientation on the first two points
  //
  // 
  // Reference:
  // Trifocal Relative Pose from Lines at Points and its Efficient Solution,
  // CVPR 2020 Ricardo Fabbri, Timothy Duff, Hongyi Fan, Margaret Regan, David
  // de Pinho, Elias Tsigaridas, Charles Wampler, Jonathan Hauenstein, Peter
  // Giblin, Benjamin Kimia, Anton Leykin and Tomas Pajdla
  static void Solve(
      const Mat &datum_0, // x, y, cos(feature orientation), sin(feature orientation) for 1st view
      const Mat &datum_1, // x, y, cos(feature orientation), sin(feature orientation) for 2nd view
      const Mat &datum_2, // x, y, cos(feature orientation), sin(feature orientation) for 3rd view
      std::vector<trifocal_model_t> *trifocal_tensor); // vector of solutions (each solution = 3 cameras)
};

} // namespace trifocal
} // namespace OpenMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP
