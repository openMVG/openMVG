// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_POSE3_H_
#define OPENMVG_GEOMETRY_POSE3_H_

#include "openMVG/multiview/projection.hpp"
namespace openMVG {
namespace geometry {

// Define a 3D Pose as a 3D transform : [R|C] t = -RC
class Pose3
{
  protected:
    Mat3 _rotation;
    Vec3 _center;
 
  public:
    // Constructors
    Pose3() : _rotation(Mat3::Identity()), _center(Vec3::Zero()) {}
    Pose3(const Mat3& r, const Vec3& c) : _rotation(r), _center(c) {}
    
    // Accessors
    const Mat3& rotation() const { return _rotation; }
    Mat3& rotation() { return _rotation; }
    const Vec3& center() const { return _center; }
    Vec3& center() { return _center; }
    
    // Translation vector t = -RC
    inline Vec3 translation() const { return -(_rotation * _center); }
    
    // Apply pose
    inline Vec3 operator () (const Vec3& p) const
    {
      return _rotation * (p - _center);
    }
 
    // Composition
    Pose3 operator * (const Pose3& P) const
    {
      return Pose3(_rotation * P._rotation, P._center + P._rotation.transpose() * _center );
    }
 
    // Inverse
    Pose3 inverse() const
    {
      return Pose3(_rotation.transpose(),  -(_rotation * _center));
    }
};
} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_POSE3_H_
