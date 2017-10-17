// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_essential_spherical.hpp"
#include "openMVG/numeric/nullspace.hpp"
#include "openMVG/numeric/extract_columns.hpp"

namespace openMVG {
namespace spherical_cam {

template<typename TMatX, typename TMatA>
static inline void EncodeEpipolarEquation
(
  const TMatX &x1,
  const TMatX &x2,
  TMatA *A
)
{
  for (int i = 0; i < x1.cols(); ++i)
  {
    A->row(i) <<
      x2(0, i) * x1.col(i).transpose(),
      x2(1, i) * x1.col(i).transpose(),
      x2(2, i) * x1.col(i).transpose();
  }
}

void EightPointRelativePoseSolver::Solve
(
  const Mat &x1,
  const Mat &x2,
  std::vector<Mat3> *pvec_E
)
{
  assert(3 == x1.rows());
  assert(8 <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  MatX9 epipolar_constraint(x1.cols(), 9);
  epipolar_constraint.fill(0.0);
  EncodeEpipolarEquation(x1, x2, &epipolar_constraint);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>> solver
    (epipolar_constraint.transpose() * epipolar_constraint);
  Vec9 e = solver.eigenvectors().leftCols<1>();
  Mat3 E = Map<RMat3>(e.data());

  // Find the closest essential matrix to E in frobenius norm
  // E = UD'VT
  if (x1.cols() > 8) {
    Eigen::JacobiSVD<Mat3> USV(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vec3 d = USV.singularValues();
    const double a = d[0];
    const double b = d[1];
    d << (a + b) / 2., (a + b) / 2., 0.0;
    E = USV.matrixU() * d.asDiagonal() * USV.matrixV().transpose();
  }
  pvec_E->emplace_back(E);
}


// Return the angular error between [0; PI/2]
double AngularError::Error
(
  const Mat3 &model,
  const Vec3 &x1,
  const Vec3 &x2
)
{
  const Vec3 Em1 = (model * x1).normalized();
  const double angleVal = (x2.transpose() * Em1);
  return std::abs(asin(angleVal));
}

EssentialKernel_spherical::EssentialKernel_spherical
(
  const Mat &x1,
  const Mat &x2
): x1_(x1), x2_(x2)
{
}

void EssentialKernel_spherical::Fit
(
  const std::vector<size_t> &samples,
  std::vector<Model> *models
) const
{
  const Mat x1 = ExtractColumns(x1_, samples);
  const Mat x2 = ExtractColumns(x2_, samples);

  assert(3 == x1.rows());
  assert(MINIMUM_SAMPLES <= x1.cols());
  assert(x1.rows() == x2.rows());
  assert(x1.cols() == x2.cols());

  EightPointRelativePoseSolver::Solve(x1, x2, models);
}

size_t EssentialKernel_spherical::NumSamples() const {return x1_.cols();}

/// Return the angular error (between 0 and PI/2)
double EssentialKernel_spherical::Error
(
  size_t sample,
  const Model &model
) const
{
  return AngularError::Error(model, x1_.col(sample), x2_.col(sample));
}


} // namespace spherical_cam
} // namespace openMVG
