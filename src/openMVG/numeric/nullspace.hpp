// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_NUMERIC_NULLSPACE_HPP
#define OPENMVG_NUMERIC_NULLSPACE_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG
{

/**
* @brief Solve linear system
*
* Linear system is given by : \n
* \f$ A x = 0 \f$
* Solution is found using the constraint on x : \f$ \| x \| = 1 \f$
*
* @param[in,out] A Input matrix storing the system to solve
* @param[out] nullspace result vector containing the solution of the system
* @return Singular value corresponding to the solution of the system
*
* @note Computation is made using SVD decomposition of input matrix
* @note Input matrix A content may be modified during computation
* @note Input vector nullspace may be resized to store the full result
*/
double Nullspace
(
  const Eigen::Ref<const Mat> & A,
  Eigen::Ref<Vec> nullspace
);

/**
* @brief Solve linear system and gives the two best solutions
*
* Linear system is given by : \n
* \f$ A x = 0 \f$
* Solution is found using the constraint on x : \f$ \| x \| = 1 \f$
*
* @param[in,out] A Input matrix storing the system to solve
* @param[out] x1 result vector containing the best solution of the system
* @param[out] x2 result vector containing the second best solution of the system
* @return Singular value corresponding to the best solution of the system
*
* @note Computation is made using SVD decomposition of input matrix
* @note Input matrix A content may be modified during computation
* @note Input vector nullspace may be resized to store the full result
*/
double Nullspace2
(
  const Eigen::Ref<const Mat> & A,
  Eigen::Ref<Vec> x1,
  Eigen::Ref<Vec> x2
);

} // namespace openMVG

#endif  // OPENMVG_NUMERIC_NULLSPACE_HPP
