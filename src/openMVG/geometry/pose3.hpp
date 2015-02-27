// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_POSE3_H_
#define OPENMVG_GEOMETRY_POSE3_H_

#include "openMVG/multiview/projection.hpp"
#include <cereal/cereal.hpp>

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

    // Serialization
    template <class Archive>
    void save( Archive & ar) const
    {
      const std::vector<std::vector<double>> mat =
        {
          { _rotation(0,0), _rotation(0,1), _rotation(0,2) },
          { _rotation(1,0), _rotation(1,1), _rotation(1,2) },
          { _rotation(2,0), _rotation(2,1), _rotation(2,2) }
        };

      ar(cereal::make_nvp("rotation", mat));

      const std::vector<double> vec = { _center(0), _center(1), _center(2) };
      ar(cereal::make_nvp("center", vec ));
    }

    template <class Archive>
    void load( Archive & ar)
    {
      std::vector<std::vector<double>> mat(3, std::vector<double>(3));
      ar(cereal::make_nvp("rotation", mat));
      // copy back to the rotation
      _rotation <<
        mat[0][0], mat[0][1], mat[0][2],
        mat[1][0], mat[1][1], mat[1][2],
        mat[2][0], mat[2][1], mat[2][2];

      std::vector<double> vec(3);
      ar(cereal::make_nvp("center", vec));
      _center << vec[0], vec[1], vec[2];
    }
};
} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_POSE3_H_
