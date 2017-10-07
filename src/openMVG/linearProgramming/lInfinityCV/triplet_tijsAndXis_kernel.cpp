// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/linearProgramming/lInfinityCV/triplet_tijsAndXis_kernel.hpp"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/numeric/numeric.h"

// Linear programming solver(s)
#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri.hpp"

namespace openMVG {
namespace trifocal {
namespace kernel {

/// A trifocal tensor seen as 3 projective cameras
double TrifocalTensorModel::Error
(
  const TrifocalTensorModel & t,
  const Vec2 & pt1,
  const Vec2 & pt2,
  const Vec2 & pt3
)
{
  // Triangulate
  const std::vector<Mat34> poses {t.P1, t.P2, t.P3};
  const std::vector<Vec3> Xs {pt1.homogeneous(),
                              pt2.homogeneous(),
                              pt3.homogeneous()};
  Eigen::Map<const Mat3> bearing_matrix(Xs[0].data());

  Vec4 Xhomogeneous;
  TriangulateNViewAlgebraic
  (
    bearing_matrix,
    poses, // Ps are projective cameras.
    &Xhomogeneous);
  const Vec3 X = Xhomogeneous.hnormalized();

  // Return the maximum observed reprojection error
  const double pt1ReProj = (Project(t.P1, X) - pt1).squaredNorm();
  const double pt2ReProj = (Project(t.P2, X) - pt2).squaredNorm();
  const double pt3ReProj = (Project(t.P3, X) - pt3).squaredNorm();

  return std::max(pt1ReProj, std::max(pt2ReProj,pt3ReProj));
}

}  // namespace kernel
}  // namespace trifocal
}  // namespace openMVG

namespace openMVG{

using openMVG::trifocal::kernel::TrifocalTensorModel;

/// Solve the computation of the "tensor".
void translations_Triplet_Solver::Solve
(
  const Mat &pt0,
  const Mat & pt1,
  const Mat & pt2,
  const std::vector<Mat3> & vec_KR,
  std::vector<TrifocalTensorModel> *P,
  const double ThresholdUpperBound
)
{
  const int n_obs = pt0.cols();
  if (n_obs < MINIMUM_SAMPLES)
    return;

  // Build the megaMatMatrix (compact form of point coordinates & their view index)
  Mat4X megaMat(4, n_obs*3);
  {
    size_t cpt = 0;
    for (int i = 0; i  < n_obs; ++i)
    {
      megaMat.col(cpt++) << pt0.col(i)(0), pt0.col(i)(1), (double)i, 0.0;
      megaMat.col(cpt++) << pt1.col(i)(0), pt1.col(i)(1), (double)i, 1.0;
      megaMat.col(cpt++) << pt2.col(i)(0), pt2.col(i)(1), (double)i, 2.0;
    }
  }
  //-- Solve the LInfinity translation and structure from Rotation and points data.
  std::vector<double> vec_solution((3 + MINIMUM_SAMPLES)*3);

  using namespace openMVG::lInfinityCV;

  OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));

  Translation_Structure_L1_ConstraintBuilder cstBuilder(vec_KR, megaMat);
  double gamma;
  if (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
    LPsolver,
    cstBuilder,
    &vec_solution,
    ThresholdUpperBound,
    0.0, 1e-8, 2, &gamma, false))
  {
    const std::vector<Vec3> vec_tis {
      {vec_solution[0], vec_solution[1], vec_solution[2]},
      {vec_solution[3], vec_solution[4], vec_solution[5]},
      {vec_solution[6], vec_solution[7], vec_solution[8]}};

    TrifocalTensorModel PTemp;
    PTemp.P1 = HStack(vec_KR[0], vec_tis[0]);
    PTemp.P2 = HStack(vec_KR[1], vec_tis[1]);
    PTemp.P3 = HStack(vec_KR[2], vec_tis[2]);

    P->push_back(PTemp);
  }
}

// Compute the residual of reprojections
double translations_Triplet_Solver::Error
(
  const TrifocalTensorModel & Tensor,
  const Vec2 & pt0,
  const Vec2 & pt1,
  const Vec2 & pt2
)
{
  return TrifocalTensorModel::Error(Tensor, pt0, pt1, pt2);
}

} // namespace openMVG
