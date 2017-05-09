// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 cDc and Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L1_HPP
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L1_HPP

#include <vector>
#include "openMVG/multiview/rotation_averaging_common.hpp"
#include "openMVG/numeric/numeric.h"

//------------------
//-- Bibliography --
//------------------
//- [1] "Efficient and Robust Large-Scale Rotation Averaging"
//- Authors: Avishek Chatterjee and Venu Madhav Govindu
//- Date: December 2013.
//- Conference: ICCV.

namespace openMVG   {
namespace rotation_averaging  {
namespace l1  {

// D E F I N E S ///////////////////////////////////////////////////
using Matrix3x3Arr = std::vector<openMVG::Mat3>;

/**
 * @brief Compute an initial estimation of global rotation (chain rotations along a MST).
 *
 * @param[in] RelRs Relative weighted rotation matrices
 * @param[out] Rs output global rotation matrices
 * @param[in] nMainViewID Id of the image considered as Identity (unit rotation)
 */
void InitRotationsMST
(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const uint32_t nMainViewID
);

/**
 * @brief Compute an initial estimation of global rotation and refines them under the L1 norm, [1].
 *
 * @param[in] RelRs Relative weighted rotation matrices
 * @param[out] Rs output global rotation matrices
 * @param[in] nMainViewID Id of the image considered as Identity (unit rotation)
 * @param[in] threshold (optionnal) threshold
 * @param[out] vec_inliers rotation labelled as inliers or outliers
 */
bool GlobalRotationsRobust(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const uint32_t nMainViewID,
  float threshold = 0.f,
  std::vector<bool> * vec_inliers = nullptr );

/**
 * @brief Implementation of Iteratively Reweighted Least Squares (IRLS) [1].
 *
 * @param[in] RelRs Relative weighted rotation matrices
 * @param[out] Rs output global rotation matrices
 * @param[in] nMainViewID Id of the image considered as Identity (unit rotation)
 * @param[in] sigma factor
 */
bool RefineRotationsAvgL1IRLS(
  const RelativeRotations& RelRs,
  Matrix3x3Arr& Rs,
  const uint32_t nMainViewID,
  const double sigma = openMVG::D2R(5));

/**
 * @brief Sort relative rotation as inlier, outlier rotations.
 *
 * @param[in] RelRs Relative weighted rotation matrices
 * @param[out] Rs output global rotation matrices
 * @param[in] threshold used to label rotations as inlier, or outlier (if 0, threshold is computed with the X84 law)
 * @param[in] vec_inliers inlier, outlier labels
 */
unsigned int FilterRelativeRotations(
  const RelativeRotations& RelRs,
  const Matrix3x3Arr& Rs,
  float threshold = 0.f,
  std::vector<bool> * vec_inliers = nullptr );

} // namespace l1
} // namespace rotation_averaging
} // namespace openMVG

#endif // OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L1_HPP
