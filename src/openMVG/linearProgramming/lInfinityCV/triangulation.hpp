// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_TRIANGULATION_HPP
#define OPENMVG_LINFINITY_COMPUTER_VISION_TRIANGULATION_HPP

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

//-- Triangulation
//    - Estimation of X from Pi and xij
// [1] -> 5.1 The triangulation problem
//
//    - This implementation Use L1 norm instead of the L2 norm of
//      the paper, it allows to use standard standard LP
//      (simplex) instead of using SOCP (second order cone programming).
//

void EncodeTriangulation
(
  const std::vector<Mat34> & Pi, // Projection matrices
  const Mat2X & x_ij, // corresponding observations
  double gamma, // Start upper bound
  Mat & A,
  Vec & C
);


/// Kernel that set Linear constraints for the Triangulation Problem.
///  Designed to be used with bisectionLP and LP_Solver interface.
///
/// Triangulation :
///    - Estimation of Xi from Pj and xij
/// Implementation of problem of [1] -> 5.1 The triangulation problem
///  under a Linear Program form.
struct Triangulation_L1_ConstraintBuilder
{
  Triangulation_L1_ConstraintBuilder
  (
    const std::vector<Mat34> & vec_Pi,
    const Mat2X & x_ij
  );

  /// Setup constraints of the triangulation problem as a Linear program
  bool Build
  (
    double gamma,
    linearProgramming::LP_Constraints & constraint
  );

  //-- Data required for triangulation :
  std::vector<Mat34> vec_Pi_;  // Projection matrix
  Mat2X x_ij_;                 // 2d projection : xij = Pj*Xi
};


} // namespace lInfinityCV
} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_TRIANGULATION_HPP
