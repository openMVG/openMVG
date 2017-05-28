// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/geometry/Similarity3_Kernel.hpp"
#include "openMVG/geometry/rigid_transformation3D_srt.hpp"

namespace openMVG{
namespace geometry{
namespace kernel{

void Similarity3Solver::Solve
(
  const Mat &x,
  const Mat &y,
  std::vector<Similarity3> *sims
)
{
  assert(3 == x.rows());
  assert(x.rows() == y.rows());
  assert(x.cols() == y.cols());

  double S;
  Vec3 t;
  Mat3 R;
  if ( FindRTS(x, y, &S, &t, &R) )
  {
    // Emplace back the Similarity3
    sims->emplace_back(geometry::Pose3(R, -R.transpose()* t/S), S);
  }
}

Vec Similarity3ErrorSquaredMetric::ErrorVec
(
  const Similarity3 &S,
  const Mat3X &x1,
  const Mat3X &x2
)
{
  return (x2 - S(x1)).colwise().squaredNorm();
}

double Similarity3ErrorSquaredMetric::Error
(
  const Similarity3 &S,
  const Vec3 &x1,
  const Vec3 &x2
)
{
  return (x2 - S(x1)).squaredNorm();
}


} // namespace kernel
} // namespace geometry
} // namespace openMVG
