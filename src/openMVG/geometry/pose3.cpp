// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/projection.hpp"

#include <vector>

namespace openMVG{
namespace geometry{

    Pose3::Pose3
    (
      const Mat3& r,
      const Vec3& c
    )
    : rotation_( r ), center_( c ) {}

    const Mat3& Pose3::rotation() const
    {
      return rotation_;
    }

    Mat3& Pose3::rotation()
    {
      return rotation_;
    }

    const Vec3& Pose3::center() const
    {
      return center_;
    }

    Vec3& Pose3::center()
    {
      return center_;
    }

    Vec3 Pose3::translation() const
    {
      return -( rotation_ * center_ );
    }


    Mat3X Pose3::operator () ( const Mat3X& p ) const
    {
      return rotation_ * ( p.colwise() - center_ );
    }


    Pose3 Pose3::operator * ( const Pose3& P ) const
    {
      return Pose3( rotation_ * P.rotation_, P.center_ + P.rotation_.transpose() * center_ );
    }


    Pose3 Pose3::inverse() const
    {
      return Pose3( rotation_.transpose(),  -( rotation_ * center_ ) );
    }


    double Pose3::depth( const Vec3 &X ) const
    {
      return ( rotation_ * ( X - center_ ) )[2];
    }

} // namespace geometry
} // namespace openMVG


