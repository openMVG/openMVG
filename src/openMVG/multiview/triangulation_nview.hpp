// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_HPP
#define OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_HPP

#include <utility>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

  /// Compute a 3D position of a point from several images of it. In particular,
  ///  compute the projective point X in R^4 such that x = PX.
  /// Algorithm is the standard DLT; for derivation see appendix of Keir's thesis.
  void TriangulateNView
  (
    const Mat3X &x, // x's are landmark bearing vectors in each camera
    const std::vector<Mat34> &Ps, // Ps are projective cameras
    Vec4 *X
  );

  // This method uses the algebraic distance approximation.
  // Note that this method works better when the 2D points are normalized
  // with an isotopic normalization.
  bool TriangulateNViewAlgebraic
  (
    const Mat3X &x, // x's are landmark bearing vectors in each camera
    const std::vector<Mat34> &Ps, // Ps are projective cameras.
    Vec4 *X
  );

}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_HPP
