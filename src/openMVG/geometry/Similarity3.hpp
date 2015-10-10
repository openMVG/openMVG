// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_SIMILARITY3_H_
#define OPENMVG_GEOMETRY_SIMILARITY3_H_


namespace openMVG {
namespace geometry {

// Define a 3D Similarity transform encoded as a 3D pose plus a scale
struct Similarity3
{
  Pose3 _pose;
  double _scale;

  Similarity3(): _pose(Pose3()), _scale(1.0) {}
  Similarity3(const Pose3 & pose, const double scale) : _pose(pose), _scale(scale) {}

  // Operators
  Vec3 operator () (const Vec3 & point) const
  {
    return _scale * _pose(point);
  }

  Pose3 operator () (const Pose3 & pose) const
  {
    return Pose3( pose.rotation() * _pose.rotation().transpose(), this->operator()(pose.center()));
  }
};

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_SIMILARITY3_H_
