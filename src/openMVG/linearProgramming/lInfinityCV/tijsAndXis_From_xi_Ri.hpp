// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_TRANSLATIONANDSTRUCTUREFrom_xi_RI_HPP
#define OPENMVG_LINFINITY_COMPUTER_VISION_TRANSLATIONANDSTRUCTUREFrom_xi_RI_HPP

#include <utility>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"

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

using namespace linearProgramming;

//-- Estimate the translation and the structure
//    from image points coordinates and camera rotations.
//    - Estimation of Ci from Ri and xij
// [1] -> 6.1 Cameras with Known Rotation
//
//    - This implementation Use L1 norm instead of the L2 norm of
//      the paper, it allows to use standard standard LP
//      (simplex) instead of using SOCP (second order cone programming).
//      Implementation by Pierre Moulon
//

/// Encode translation and structure linear program
void EncodeTiXi
(
  const Mat & M, //Scene representation
  const std::vector<Mat3> & Ri,
  double sigma, // Start upper bound
  sRMat & A, Vec & C,
  std::vector<LP_Constraints::eLP_SIGN> & vec_sign,
  std::vector<double> & vec_costs,
  std::vector< std::pair<double,double> > & vec_bounds
);



/// Kernel that set Linear constraints for the
///   - Translation Registration and Structure Problem.
///  Designed to be used with bisectionLP and LP_Solver interface.
///
/// Implementation of problem of [1] -> 6.1 Cameras with known rotation
///  under a Linear Program form. (With SPARSE constraint matrix).

struct Translation_Structure_L1_ConstraintBuilder
{
  Translation_Structure_L1_ConstraintBuilder
  (
    const std::vector<Mat3> & vec_Ri,
    const Mat & M
  );

  /// Setup constraints for the translation and structure problem,
  ///  in the LP_Constraints object.
  bool Build
  (
    double gamma,
    LP_Constraints_Sparse & constraint
  );

  std::vector<Mat3> vec_Ri_;  // Rotation matrix
  Mat M_; // M contains (X,Y,index3dPoint, indexCam)^T
};

} // namespace lInfinityCV
} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_TRANSLATIONANDSTRUCTUREFrom_xi_RI_HPP
