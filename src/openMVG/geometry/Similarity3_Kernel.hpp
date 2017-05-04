// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_SIMILARITY3_KERNEL_HPP
#define OPENMVG_GEOMETRY_SIMILARITY3_KERNEL_HPP

#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/multiview/two_view_kernel.hpp"

#include <vector>

namespace openMVG
{
namespace geometry
{
namespace kernel
{

struct Similarity3Solver
{
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 1 };

  /**
   * Computes the 3D similarity transform between two point cloud
   *
   * \param x  A 3xN matrix of column vectors.
   * \param y  A 3xN matrix of column vectors.
   * \param sim The found similarity
   *
   * The estimated 3D similarity should approximately hold the condition y = sim(x).
   */
  static void Solve
  (
    const Mat &x,
    const Mat &y,
    std::vector<Similarity3> *sims
  );
};

struct Similarity3ErrorSquaredMetric
{
  // Return the Squared error between a collection of points (stored as column)
  static Vec ErrorVec
  (
    const Similarity3 &S,
    const Mat3X &x1,
    const Mat3X &x2
  );

  // Return the Squared error between the point x2 and the transformed point S(x1)
  static double Error
  (
    const Similarity3 &S,
    const Vec3 &x1,
    const Vec3 &x2
  );
};

// Define a Kernel to solve a robust 3D similarity between point cloud
using Similarity3_Kernel = two_view::kernel::Kernel
  <
    Similarity3Solver,              // The model solver
    Similarity3ErrorSquaredMetric,  // The datum to model error metric
    Similarity3                     // The model type
  >;

} // namespace kernel
} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_SIMILARITY3_KERNEL_HPP
