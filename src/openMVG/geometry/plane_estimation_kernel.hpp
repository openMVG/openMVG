// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_PLANE_ESTIMATION_KERNEL_HPP
#define OPENMVG_GEOMETRY_PLANE_ESTIMATION_KERNEL_HPP

#include <cmath>
#include <utility>

#include "openMVG/geometry/half_space_intersection.hpp"
#include "openMVG/numeric/extract_columns.hpp"

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(
  std::pair<openMVG::geometry::halfPlane::Half_plane,
            openMVG::Vec3>)

namespace openMVG
{
namespace geometry
{
namespace plane
{

struct PlaneSolver
{
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 1 };

  /**
  * Compute a Half_plane equation from 3 points.
  *
  * \param[in] x        The input points (3x3 matrix).
  * \param[out] plane   The computed half_plane.
  */
  static void Solve
  (
    const Mat3X &x,
    halfPlane::Half_planes *plane
  )
  {
    plane->push_back(halfPlane::Half_plane_p(x.col(0), x.col(1), x.col(2)));
  }

  /**
  * Compute a Half_plane equation from 3 points and store its centroid.
  *
  * \param[in] x            The input points (3x3 matrix).
  * \param[out] local_plane The computed half_plane and centroid of the points.
  */
  static void Solve
  (
    const Mat3X &x,
    std::vector<std::pair<halfPlane::Half_plane, Vec3>> *local_plane
  )
  {
    local_plane->push_back(
      {
        halfPlane::Half_plane_p(x.col(0), x.col(1), x.col(2)),
        (x.col(0) + x.col(1) + x.col(2)) / 3.0
      }
    );
  }
};

struct AbsDistanceError
{
  /**
  * Compute the absolute distance of the point to the plane.
  *
  * \param[in] half_plane   The plane.
  * \param[in] point        The input point.
  */
  static double Error
  (
    const halfPlane::Half_plane &half_plane,
    const Vec3 &point
  )
  {
    return half_plane.absDistance(point);
  }
};

struct AbsAngularError
{
  /**
  * Compute an angular error (radians) between the direction of the point
  * and the plane normal.
  * The metric is biased by the distance:
  *  - This cost function will prefer plane with large region support.
  *
  * \param[in] half_plane   The plane and its centroid.
  * \param[in] point        The input point.
  */
  static double Error
  (
    const std::pair<halfPlane::Half_plane, Vec3> &local_plane,
    const Vec3 &point
  )
  {
    const Vec3 bearing_direction = (point - local_plane.second).normalized();
    return std::abs(acos(bearing_direction.dot(local_plane.first.normal())) - M_PI / 2.0);
  }
};

// Define a kernel to be able to perform robust estimation of Half_plane.
// Use a half_plane Solver and a point to plane distance cost function.
struct HaflPlaneKernel
{
  using Model = halfPlane::Half_plane;
  enum { MINIMUM_SAMPLES = 3 };

  explicit HaflPlaneKernel(const Mat3X &points) : points_(points) {}

  size_t NumSamples() const { return static_cast<size_t> (points_.cols()); }

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<Model> *planes
  )
  const
  {
    assert(samples.size() >= (unsigned int)MINIMUM_SAMPLES);
    const Mat3X sampled_xs = ExtractColumns(points_, samples);

    PlaneSolver::Solve(sampled_xs, planes);
  }

  double Error
  (
    uint32_t sample,
    const Model &model
  )
  const
  {
    return AbsDistanceError::Error(model, points_.col(sample));
  }

  const Mat3X &points_;
};

// Define a kernel to be able to perform robust estimation of Half_plane.
// Use a half_plane and centroid Solver and a point to plane angular cost function.
struct HaflPlaneKernelAngular
{
  using Model = std::pair<halfPlane::Half_plane, Vec3>;
  enum { MINIMUM_SAMPLES = 3 };

  explicit HaflPlaneKernelAngular(const Mat3X &points) : points_(points) {}

  size_t NumSamples() const { return static_cast<size_t> (points_.cols()); }

  void Fit
  (
    const std::vector<uint32_t> &samples,
    std::vector<Model> *planes
  )
  const
  {
    assert(samples.size() >= (unsigned int)MINIMUM_SAMPLES);
    const Mat3X sampled_xs = ExtractColumns(points_, samples);

    PlaneSolver::Solve(sampled_xs, planes);
  }

  double Error
  (
    uint32_t sample,
    const Model &model
  )
  const
  {
    return AbsAngularError::Error(model, points_.col(sample));
  }

  const Mat3X &points_;
};


} // namespace plane
} // namespace geometry
} // namespace openMVG
#endif // OPENMVG_GEOMETRY_PLANE_ESTIMATION_KERNEL_HPP
