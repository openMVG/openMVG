// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/Similarity3.hpp"

namespace openMVG{
namespace geometry{

Similarity3::Similarity3()
  : pose_( Pose3() ),
    scale_( 1.0 )
{

}

Similarity3::Similarity3( const Pose3 & pose, const double scale )
  : pose_( pose ),
    scale_( scale )
{

}

Mat3X Similarity3::operator () ( const Mat3X & point ) const
{
  return scale_ * pose_( point );
}

Pose3 Similarity3::operator () ( const Pose3 & pose ) const
{
  return Pose3( pose.rotation() * pose_.rotation().transpose(), this->operator()( pose.center() ) );
}

Similarity3 Similarity3::inverse() const
{
  return Similarity3(pose_.inverse(), 1.0 / scale_);
}

} // namespace geometry
} // namespace openMVG
