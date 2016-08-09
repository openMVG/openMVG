// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_POSE3_H_
#define OPENMVG_GEOMETRY_POSE3_H_

#include "openMVG/multiview/projection.hpp"
#include <cereal/cereal.hpp> // Serialization

namespace openMVG {
namespace geometry {

// Define a 3D Pose as a 3D transform: [R|C] t = -RC
class Pose3
{
  protected:
    Mat3 _rotation;
    Vec3 _center;

  public:
    // Constructors
    Pose3() : _rotation(Mat3::Identity()), _center(Vec3::Zero()) {}
    Pose3(const Mat3& r, const Vec3& c) : _rotation(r), _center(c) {}
    Pose3(const Mat34& Rt)
    : _rotation(Rt.block<3,3>(0,0))
    {
      Vec3 t = Rt.block<3, 1>(0,3);
      _center = -_rotation.transpose() * t;
    }

    // Accessors
    const Mat3& rotation() const { return _rotation; }
    Mat3& rotation() { return _rotation; }
    const Vec3& center() const { return _center; }
    Vec3& center() { return _center; }

    // Translation vector t = -RC
    inline Vec3 translation() const { return -(_rotation * _center); }

    // Apply pose
    inline Mat3X operator () (const Mat3X& p) const
    {
      return _rotation * (p.colwise() - _center);
    }

    // Composition
    Pose3 operator * (const Pose3& P) const
    {
      return Pose3(_rotation * P._rotation, P._center + P._rotation.transpose() * _center );
    }

    // Operator ==
    bool operator==(const Pose3& other) const
    {
      return AreMatNearEqual(_rotation, other._rotation, 1e-6) &&
              AreVecNearEqual(_center, other._center, 1e-6);
    }

    // Inverse
    Pose3 inverse() const
    {
      return Pose3(_rotation.transpose(),  -(_rotation * _center));
    }

    /// Return the depth (distance) of a point respect to the camera center
    double depth(const Vec3 &X) const
    {
      return (_rotation * (X - _center))[2];
    }

    /// Return the pose with the Scale, Rotation and translation applied.
    Pose3 transformSRt(const double S, const Mat3 & R, const Vec3 & t) const
    {
      Pose3 pose;
      pose._center = S * R * _center + t;
      pose._rotation = _rotation * R.transpose();
      return pose;
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
      _rotation.row(0) = Eigen::Map<const Vec3>(&(mat[0][0]));
      _rotation.row(1) = Eigen::Map<const Vec3>(&(mat[1][0]));
      _rotation.row(2) = Eigen::Map<const Vec3>(&(mat[2][0]));

      std::vector<double> vec(3);
      ar(cereal::make_nvp("center", vec));
      _center = Eigen::Map<const Vec3>(&vec[0]);
    }
};

/**
 * @brief Build a pose from a rotation and a translation.
 * @param[in] R The 3x3 rotation.
 * @param[in] t The 3x1 translation.
 * @return The pose as [R, -R'*t]
 */
inline Pose3 poseFromRT(const Mat3& R, const Vec3& t) 
{
  return Pose3(R, -R.transpose()*t);
}

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_POSE3_H_
