// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/numeric/nullspace.hpp"

#include <limits>

namespace openMVG {

void TriangulateNView
(
  const Mat3X &x,
  const std::vector< Mat34 > &Ps,
  Vec4 *X
)
{
  const Mat2X::Index nviews = x.cols();
  assert(static_cast<size_t>(nviews) == Ps.size());

  Mat A = Mat::Zero(3 * nviews, 4 + nviews);
  for (int i = 0; i < nviews; ++i)
  {
    A.block<3, 4>(3 * i, 0)     = -Ps[i];
    A.block<3,1> (3 * i, 4 + i) = x.col(i);
  }
  Vec X_and_alphas(4 + nviews);
  Nullspace(A, X_and_alphas);
  *X = X_and_alphas.head(4);
}

void TriangulateNViewAlgebraic
(
  const Mat3X &x,
  const std::vector< Mat34 > &Ps,
  Vec4 *X
)
{
  const Mat2X::Index nviews = x.cols();
  assert(static_cast<size_t>(nviews) == Ps.size());

  // Solve the stacking of this equation:
  // [cross(x0,P0) X = 0]
  // [cross(x1,P1) X = 0]
  // [cross(x2,P2) X = 0]
  // ...
  Mat A(2*nviews, 4);
  for (int i = 0; i < nviews; ++i)
  {
    A.block<2, 4>(2*i, 0) <<
      x.col(i)[0] * Ps[i].row(2) - x.col(i)[2] * Ps[i].row(0),
      x.col(i)[1] * Ps[i].row(2) - x.col(i)[2] * Ps[i].row(1);
  }
  Nullspace(A, *X);
}

void Triangulation::add
(
  const Mat34& projMatrix,
  const Vec2 & p
)
{
  views.emplace_back( projMatrix, p );
}

size_t Triangulation::size() const
{
  return views.size();
}

void Triangulation::clear()
{
  views.clear();
}

double Triangulation::error(const Vec3 &X) const
{
  double squared_reproj_error = 0.0;
  for (const auto & view_it : views)
  {
    const Mat34& PMat = view_it.first;
    const Vec2 & xy = view_it.second;
    squared_reproj_error += (xy - Project(PMat, X)).norm();
  }
  return squared_reproj_error;
}

// Camera triangulation using the iterated linear method
Vec3 Triangulation::compute(int iter) const
{
  const int nviews(views.size());
  assert(nviews>=2);

  // Iterative weighted linear least squares
  Mat3 AtA;
  Vec3 Atb, X;
  Vec weights = Vec::Constant(nviews, 1.0);
  for (int it=0; it<iter; ++it)
  {
    AtA.fill(0.0);
    Atb.fill(0.0);
    for (int i=0; i<nviews; ++i)
    {
      const Mat34& PMat = views[i].first;
      const Vec2 & p = views[i].second;
      const double w = weights[i];

      Vec3 v1, v2;
      for (int j=0; j<3; ++j)
      {
        v1[j] = w * ( PMat(0,j) - p(0) * PMat(2,j) );
        v2[j] = w * ( PMat(1,j) - p(1) * PMat(2,j) );
        Atb[j] += w *
               ( v1[j] * ( p(0) * PMat(2,3) - PMat(0,3) )
               + v2[j] * ( p(1) * PMat(2,3) - PMat(1,3) ) );
      }

      for (int k=0; k<3; ++k)
      {
        for (int j=0; j<=k; ++j)
        {
          const double v = v1[j] * v1[k] + v2[j] * v2[k];
          AtA(j,k) += v;
          if (j<k) AtA(k,j) += v;
        }
      }
    }

    X = AtA.inverse() * Atb;

    // Compute reprojection error, min and max depth, and update weights
    zmin = std::numeric_limits<double>::max();
    zmax = std::numeric_limits<double>::lowest();
    err = 0.0;
    for (int i=0; i<nviews; ++i)
    {
      const Mat34& PMat = views[i].first;
      const Vec2 & p = views[i].second;
      const Vec3 xProj = PMat * X.homogeneous();
      const double z = xProj(2);
      if (z < zmin) zmin = z;
      else if (z > zmax) zmax = z;
      err += (p - xProj.hnormalized()).norm(); // residual error
      weights[i] = 1.0 / z;
    }
  }
  return X;
}


}  // namespace openMVG
