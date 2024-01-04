// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_RESECTION_METRICS_HPP
#define OPENMVG_MULTIVIEW_RESECTION_METRICS_HPP

namespace openMVG {
namespace resection {

struct PixelReprojectionError {
  // Compute the residual of the projection distance(x, P(X))
  static inline double Error
  (
    const Mat34 & P,
    const Vec2 & x,
    const Vec3 & X
  )
  {
    return (x - (P * X.homogeneous()).hnormalized()).norm();
  }
};

struct SquaredPixelReprojectionError {
  // Compute the Square residual of the projection distance(x, P(X))
  static inline double Error
  (
    const Mat34 & P,
    const Vec2 & x,
    const Vec3 & X
  )
  {
    return (x - (P * X.homogeneous()).hnormalized()).squaredNorm();
  }
};

struct AngularReprojectionError {
  // Compute the angular residual of the projection and the bearing vector
  static inline double Error
  (
    const Mat34 & P,
    const Vec3 & bearing_vector,
    const Vec3 & X
  )
  {
    const auto new_bearing = (P * X.homogeneous()).normalized();
    return 1.0 - (bearing_vector.dot(new_bearing));
  }
};

}  // namespace resection
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_RESECTION_METRICS_HPP
