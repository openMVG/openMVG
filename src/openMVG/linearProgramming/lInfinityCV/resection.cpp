// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/lInfinityCV/resection.hpp"

#include <limits>

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
)
{
  // Build Constraint matrix.

  const int Nobs = Pt2D.cols();

  assert(Pt2D.cols() == Pt3D.cols());
  assert(Pt2D.cols() >= 6 && "The problem requires at least 6 points");

  A.resize(Nobs * 5, 11);

  C.resize(Nobs * 5, 1);
  C.fill(0.0);

  for (int p=0; p <Nobs; ++p)
  {
    const Vec2 pt2d = Pt2D.col(p);
    const Vec3 pt3d = Pt3D.col(p);

    // Compute and setup constraint:

    // Cheirality
    // R^3 X + t >0
    // - R^3 X - t < 0 => // - R^3 X < 1.0 because t3 = 1.0 by default
    int cpti = 8;
    int cptj = 5 * p;
    A.coeffRef(cptj, cpti) = -pt3d(0);
    A.coeffRef(cptj, cpti+1) = -pt3d(1);
    A.coeffRef(cptj, cpti+2) = -pt3d(2);
    C(cptj) = 1.0;

    // First constraint <= s
    // R1 * X + t1 + R3 * X *(-xij -s) <= xij t3 + s t3
    cpti = 0;
    cptj += 1;
    A.coeffRef(cptj, cpti) = pt3d(0);
    A.coeffRef(cptj, cpti+1) = pt3d(1);
    A.coeffRef(cptj, cpti+2) = pt3d(2);
    A.coeffRef(cptj, cpti+3) = 1.0;
    cpti = 4;
    A.coeffRef(cptj+1, cpti) = pt3d(0);
    A.coeffRef(cptj+1, cpti+1) = pt3d(1);
    A.coeffRef(cptj+1, cpti+2) = pt3d(2);
    A.coeffRef(cptj+1, cpti+3) = 1.0;

    cpti = 8;
    Mat temp;
    temp = Vec2((-pt2d).array() - gamma) * pt3d.transpose();
    A.coeffRef(cptj, cpti) = temp(0,0);
    A.coeffRef(cptj, cpti+1) = temp(0,1);
    A.coeffRef(cptj, cpti+2) = temp(0,2);

    A.coeffRef(cptj+1, cpti) = temp(1,0);
    A.coeffRef(cptj+1, cpti+1) = temp(1,1);
    A.coeffRef(cptj+1, cpti+2) = temp(1,2);

    C(cptj) = gamma + pt2d(0);
    C(cptj+1) = gamma + pt2d(1);


    // Second constraint >= s
    // -R1 * X - t1 + R3 * X *(xij -s) <= - xij t3 + s t3
    cpti = 0;
    cptj += 2;
    A.coeffRef(cptj, cpti) = - pt3d(0);
    A.coeffRef(cptj, cpti+1) = - pt3d(1);
    A.coeffRef(cptj, cpti+2) = - pt3d(2);
    A.coeffRef(cptj, cpti+3) = - 1.0;
    cpti = 4;
    A.coeffRef(cptj+1, cpti) = - pt3d(0);
    A.coeffRef(cptj+1, cpti+1) = - pt3d(1);
    A.coeffRef(cptj+1, cpti+2) = - pt3d(2);
    A.coeffRef(cptj+1, cpti+3) = - 1.0;

    cpti = 8;
    temp = Vec2(pt2d.array() - gamma) * pt3d.transpose();
    A.coeffRef(cptj, cpti) = temp(0,0);
    A.coeffRef(cptj, cpti+1) = temp(0,1);
    A.coeffRef(cptj, cpti+2) = temp(0,2);

    A.coeffRef(cptj+1, cpti) = temp(1,0);
    A.coeffRef(cptj+1, cpti+1) = temp(1,1);
    A.coeffRef(cptj+1, cpti+2) = temp(1,2);

    C(cptj) = gamma - pt2d(0);
    C(cptj+1) = gamma - pt2d(1);
  }
}

Resection_L1_ConstraintBuilder::Resection_L1_ConstraintBuilder
(
  const Mat2X & Pt2D,
  const Mat3X & Pt3D
)
: pt_2d_ (Pt2D),
  pt_3d_ (Pt3D)
{

}

/// Setup constraints for the Resection problem,
///  in the LP_Constraints object.
bool Resection_L1_ConstraintBuilder::Build
(
  double gamma,
  linearProgramming::LP_Constraints_Sparse & constraint
)
{
  EncodeResection(pt_2d_, pt_3d_,
    gamma,
    constraint.constraint_mat_,
    constraint.constraint_objective_);

  //-- Setup additional information about the Linear Program constraint
  // We look for:
  // P = [P00 P01 P02 P03;
  //      P10 P11 P12 P13;
  //      P20 P21 P22 1.0];
  const int NParams = 4 * 2 + 3;

  constraint.nbParams_ = NParams;
  constraint.vec_bounds_ = { {std::numeric_limits<double>::lowest(),
                              std::numeric_limits<double>::max()} };
  // Constraint sign are all LESS or equal (<=)
  constraint.vec_sign_.resize(constraint.constraint_mat_.rows());
  std::fill(constraint.vec_sign_.begin(), constraint.vec_sign_.end(),
    linearProgramming::LP_Constraints::LP_LESS_OR_EQUAL);

  return true;
}

} // namespace lInfinityCV
} // namespace openMVG

