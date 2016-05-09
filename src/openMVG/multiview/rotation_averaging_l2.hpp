
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L2_H_
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L2_H_

#include "openMVG/multiview/rotation_averaging_common.hpp"

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

//--
//-- Implementation related to rotation averaging.
// . Compute global rotation from a list of relative estimates.
//
//- Implementation of algorithm from Thesis titled:
//- [1] "Robust Multiview Reconstruction."
//- Author : Daniel Martinec.
//- Date : July 2, 2008.
//--
namespace openMVG   {
namespace rotation_averaging  {
namespace l2  {

// [1] 6.7.2 Consistent Rotation page 89
// Closest Rotation Estimation R = U*transpose(V)
//  approximate rotation in the Frobenius norm using SVD
Mat3 ClosestSVDRotationMatrix(const Mat3 & rotMat);

//-- Solve the Global Rotation matrix registration for each camera given a list
//    of relative orientation using matrix parametrization
//    [1] formula 6.62 page 100. Dense formulation.
//- nCamera:               The number of camera to solve
//- vec_rotationEstimate:  The relative rotation i->j
//- vec_ApprRotMatrix:     The output global rotation

// Minimization of the norm of:
// => || wij * (rj - Rij * ri) ||= 0
// With rj et rj the global rotation and Rij the relative rotation from i to j.
//
// Example:
// 0_______2
//  \     /
//   \   /
//    \ /
//     1
//
// nCamera = 3
// vector.add( RelativeRotation(0,1, R01) );
// vector.add( RelativeRotation(1,2, R12) );
// vector.add( RelativeRotation(0,2, R02) );
//
bool L2RotationAveraging( size_t nCamera,
  const RelativeRotations& vec_relativeRot,
  // Output
  std::vector<Mat3> & vec_ApprRotMatrix);

// None linear refinement of the rotation using an angle-axis representation
bool L2RotationAveraging_Refine(
  const RelativeRotations & vec_relativeRot,
  std::vector<openMVG::Mat3> & vec_ApprRotMatrix);

} // namespace l2
} // namespace rotation_averaging
} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L2_H_

