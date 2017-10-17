// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/lInfinityCV/triangulation.hpp"

#include <limits>
#include <utility>
#include <vector>

#include "openMVG/numeric/numeric.h"

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

using namespace openMVG::linearProgramming;

void EncodeTriangulation
(
  const std::vector<Mat34> & Pi, // Projection matrices
  const Mat2X & x_ij, // corresponding observations
  double gamma, // Start upper bound
  Mat & A,
  Vec & C
)
{
  // Build A, C matrix.

  const size_t nbCamera = Pi.size();
  A.resize(5*nbCamera,3);
  C.resize(5*nbCamera,1);

  int cpt = 0;
  for (size_t i = 0; i < nbCamera; ++i)
  {
    const Mat3 R = Pi[i].block<3,3>(0,0);
    const Vec3 t = Pi[i].block<3,1>(0,3);
    const Mat2X pt = x_ij.col(i);

    // A (Rotational part):
    A.block<1,3>(cpt,0)   = R.row(0) - pt(0) * R.row(2) - gamma * R.row(2);
    A.block<1,3>(cpt+1,0) = R.row(1) - pt(1) * R.row(2) - gamma * R.row(2);
    A.block<1,3>(cpt+2,0) = - R.row(2);
    A.block<1,3>(cpt+3,0) = - R.row(0) + pt(0) * R.row(2) - gamma * R.row(2);
    A.block<1,3>(cpt+4,0) = - R.row(1) + pt(1) * R.row(2) - gamma * R.row(2);

    // C (translation part):
    C(cpt)   = - t(0) + pt(0) * t(2) + gamma * t(2);
    C(cpt+1) = - t(1) + pt(1) * t(2) + gamma * t(2);
    C(cpt+2) = t(2);
    C(cpt+3) = t(0) - pt(0) * t(2) + gamma * t(2);
    C(cpt+4) = t(1) - pt(1) * t(2) + gamma * t(2);

    //- Next entry
    cpt += 5;
  }
}

Triangulation_L1_ConstraintBuilder::Triangulation_L1_ConstraintBuilder
(
  const std::vector<Mat34> & vec_Pi,
  const Mat2X & x_ij
)
: vec_Pi_(vec_Pi),
  x_ij_(x_ij)
{

}

/// Setup constraints of the triangulation problem as a Linear program
bool Triangulation_L1_ConstraintBuilder::Build
(
  double gamma,
  LP_Constraints & constraint
)
{
  EncodeTriangulation(vec_Pi_, x_ij_,
    gamma,
    constraint.constraint_mat_,
    constraint.constraint_objective_);
  //-- Setup additional information about the Linear Program constraint
  // We look for 3 variables [X,Y,Z] with no bounds.
  // Constraint sign are all less or equal (<=)
  constraint.nbParams_ = 3;
  constraint.vec_bounds_ = {{std::numeric_limits<double>::lowest(),
                             std::numeric_limits<double>::max()}};
  // Setup constraint sign
  constraint.vec_sign_.resize(constraint.constraint_mat_.rows());
  fill(constraint.vec_sign_.begin(), constraint.vec_sign_.end(),
    LP_Constraints::LP_LESS_OR_EQUAL);

  return true;
}

} // namespace lInfinityCV
} // namespace openMVG
