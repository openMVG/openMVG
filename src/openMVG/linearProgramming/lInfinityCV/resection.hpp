// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_RESECTION_HPP
#define OPENMVG_LINFINITY_COMPUTER_VISION_RESECTION_HPP

#include <vector>

#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

//--
//- Implementation of algorithm from Paper titled :
//- [1] "Multiple-View Geometry under the L_\infty Norm."
//- Author : Fredrik Kahl, Richard Hartley.
//- Date : 9 sept 2008.

//- [2] "Multiple View Geometry and the L_\infty -norm"
//- Author : Fredrik Kahl
//- ICCV 2005.
//--

namespace openMVG   {
namespace lInfinityCV  {


//-- Camera Resection
//    - Estimation of [Ri|Ti] from xij and Xi
// [1] -> 5.4 Camera Resectioning
//
//  The projection is parametrized as :
//      | p11 p12 p13 p14 |
//   P= | p21 p22 p23 p24 |
//      | p31 p32 p33 1   |
//
//    - This implementation Use L1 norm instead of the L2 norm of
//      the paper, it allows to use standard standard LP
//      (simplex) instead of using SOCP (second order cone programming).
//      Implementation by Pierre Moulon
//
void EncodeResection
(
  const Mat2X & Pt2D,
  const Mat3X & Pt3D,
  double gamma, // Start upper bound
  sRMat & A, Vec & C
);

/// Kernel that set Linear constraints for the
///   - Translation Registration and Structure Problem.
///  Designed to be used with bisectionLP and LP_Solver interface.
///
/// Implementation of camera Resection
///    - Estimation of [Ri|Ti] from xij and Xi
/// [1] -> 5.4 Camera Resectioning
///
struct Resection_L1_ConstraintBuilder
{
  Resection_L1_ConstraintBuilder
  (
    const Mat2X & Pt2D,
    const Mat3X & Pt3D
  );

  /// Setup constraints for the Resection problem,
  ///  in the LP_Constraints object.
  bool Build
  (
    double gamma,
    linearProgramming::LP_Constraints_Sparse & constraint
  );

  Mat2X pt_2d_;
  Mat3X pt_3d_;
};

} // namespace lInfinityCV
} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_RESECTION_HPP
