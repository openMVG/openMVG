// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri_noise.hpp"

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
//
//  This implementation handle noisy measurement by adding a slack
//   variables for each x,y,z residual.
//  Based on idea expressed in (See Algorithm 2.0 of [3]):
// [3] Outlier Removal Using Duality.
//      Carl Olsson, Anders Eriksson and Richard Hartley, Richard.
//      CVPR 2010.

/// Encode translation and structure linear program with slack variables
///  in order to handle noisy measurements.
void EncodeTiXi_withNoise
(
  const Mat & M, //Scene representation
  const std::vector<Mat3> & Ri,
  double sigma, // Start upper bound
  sRMat & A, Vec & C,
  std::vector<LP_Constraints::eLP_SIGN> & vec_sign,
  std::vector<double> & vec_costs,
  std::vector< std::pair<double,double> > & vec_bounds
)
{
  // Build Constraint matrix.

  const size_t Ncam = (size_t) M.row(3).maxCoeff()+1;
  const size_t N3D  = (size_t) M.row(2).maxCoeff()+1;
  const size_t Nobs = M.cols();

  assert(Ncam == Ri.size());

  A.resize(5*Nobs+3, 3 * (Ncam + N3D + Nobs));

  C.resize(5 * Nobs + 3, 1);
  C.fill(0.0);
  vec_sign.resize(5*Nobs+3);

  const size_t transStart  = 0;
  const size_t pointStart  = transStart + 3*Ncam;
  const size_t offsetStart = pointStart + 3*N3D;
# define TVAR(i, el) (0 + 3*(i) + (el))
# define XVAR(j, el) (pointStart + 3*(j) + (el))
# define OFSVAR(k, el) (offsetStart + 3*(k) + (el))

  // Set objective coefficients.
  vec_costs = std::vector<double>(3 * (Ncam + N3D + Nobs), 0.0);
  for (size_t k = 0; k < Nobs; ++k)
  {
     vec_costs[OFSVAR(k, 0)] = 1.0;
     vec_costs[OFSVAR(k, 1)] = 1.0;
     vec_costs[OFSVAR(k, 2)] = 1.0;
  }

  // By default set free variable:
  vec_bounds.assign( 3 * (Ncam + N3D + Nobs),
    {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()});
  // Change the offset to be positive
  for (size_t k = 0; k < 3*Nobs; ++k)
    vec_bounds[offsetStart + k].first = 0;

  size_t rowPos = 0;
  // Add the cheirality conditions (R_i*X_j + T_i)_3 >= 1
  for (size_t k = 0; k < Nobs; ++k)
  {
    const size_t indexPt3D = M(2,k);
    const size_t indexCam  = M(3,k);
    const Mat3 & R = Ri[indexCam];

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = R(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = R(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = R(2,2);
    A.coeffRef(rowPos, TVAR(indexCam, 2)) = 1.0;
    A.coeffRef(rowPos, OFSVAR(k, 2)) = 1.0;
    C(rowPos) = 1.0;
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;
  } // end for (k)

  // Add conic constraint:
  for (size_t k = 0; k < Nobs; ++k)
  {
    const Vec2 pt   = M.block<2,1>(0,k);
    const double u = pt(0);
    const double v = pt(1);
    const size_t indexPt3D = M(2,k);
    const size_t indexCam  = M(3,k);
    const Mat3 & R = Ri[indexCam];

    // x-residual =>
    // (R_i*X_j + T_i)_1 / (R_i*X_j + T_i)_3 - u >= -sigma
    // (R_i*X_j + T_i)_1 - u * (R_i*X_j + T_i)_3  + sigma (R_i*X_j + T_i)_3  >= 0.0
    // R_i_3 * (sigma-u) + R_i_1 + t_i_1 + t_i_3 * (sigma-u) >= 0

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = R(0,0) + (sigma-u) * R(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = R(0,1) + (sigma-u) * R(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = R(0,2) + (sigma-u) * R(2,2);
    A.coeffRef(rowPos, TVAR(indexCam, 0)) = 1.0;
    A.coeffRef(rowPos, TVAR(indexCam, 2)) = sigma-u;
    A.coeffRef(rowPos, OFSVAR(k, 0)) = 1.0;
    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = R(0,0) - (sigma+u) * R(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = R(0,1) - (sigma+u) * R(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = R(0,2) - (sigma+u) * R(2,2);
    A.coeffRef(rowPos, TVAR(indexCam, 0)) = 1.0;
    A.coeffRef(rowPos, TVAR(indexCam, 2)) = -(sigma + u);
    A.coeffRef(rowPos, OFSVAR(k, 0)) = -1.0;
    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;

    // y-residual =>
    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = R(1,0) + (sigma-v) * R(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = R(1,1) + (sigma-v) * R(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = R(1,2) + (sigma-v) * R(2,2);
    A.coeffRef(rowPos, OFSVAR(k, 1)) = 1.0;
    A.coeffRef(rowPos, TVAR(indexCam, 1)) = 1.0;
    A.coeffRef(rowPos, TVAR(indexCam, 2)) = sigma-v;

    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = R(1,0) - (sigma+v) * R(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = R(1,1) - (sigma+v) * R(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = R(1,2) - (sigma+v) * R(2,2);
    A.coeffRef(rowPos, TVAR(indexCam, 1)) = 1.0;
    A.coeffRef(rowPos, TVAR(indexCam, 2)) = -(sigma + v);
    A.coeffRef(rowPos, OFSVAR(k, 1)) = -1.0;
    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;
  }
  // Fix the translation ambiguity. (set first cam at (0,0,0)
  //LP_EQUAL
  A.coeffRef(rowPos, TVAR(0, 0)) = 1.0; C(rowPos) = 0.0; vec_sign[rowPos] = LP_Constraints::LP_EQUAL; ++rowPos;
  A.coeffRef(rowPos, TVAR(0, 1)) = 1.0; C(rowPos) = 0.0; vec_sign[rowPos] = LP_Constraints::LP_EQUAL; ++rowPos;
  A.coeffRef(rowPos, TVAR(0, 2)) = 1.0; C(rowPos) = 0.0; vec_sign[rowPos] = LP_Constraints::LP_EQUAL; ++rowPos;

# undef TVAR
# undef XVAR
# undef OFSVAR
}

TiXi_withNoise_L1_ConstraintBuilder::TiXi_withNoise_L1_ConstraintBuilder(
  const std::vector<Mat3> & vec_Ri,
  const Mat & M)
{
  _M = M;
  _vec_Ri = vec_Ri;
}

/// Setup constraints for the translation and structure problem,
///  in the LP_Constraints object.
bool TiXi_withNoise_L1_ConstraintBuilder::Build
(
  double gamma,
  LP_Constraints_Sparse & constraint
)
{
  EncodeTiXi_withNoise(_M, _vec_Ri,
    gamma,
    constraint.constraint_mat_,
    constraint.constraint_objective_,
    constraint.vec_sign_,
    constraint.vec_cost_,
    constraint.vec_bounds_);

  constraint.bminimize_ = true;

  //-- Setup additional information about the Linear Program constraint
  // We look for nb translations and nb 3D points.
  const size_t N3D  = (size_t) _M.row(2).maxCoeff() + 1;
  const size_t Nobs = _M.cols();
  const size_t NNoiseVar = Nobs; // noise over Zdepth, (x,y) residual
  const size_t Ncam = (size_t) _M.row(3).maxCoeff() + 1;

  constraint.nbParams_ = (Ncam + N3D + NNoiseVar)*3;

  return true;
}

} // namespace lInfinityCV
} // namespace openMVG
