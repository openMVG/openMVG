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

} // namespace openMVG

#endif  // OPENMVG_NUMERIC_NULLSPACE_HPP
