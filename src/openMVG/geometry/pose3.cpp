// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/geometry/pose3.hpp"
#include "openMVG/multiview/projection.hpp"

#include <cereal/types/vector.hpp>
#include <cereal/cereal.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

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

template <class Archive>
void Pose3::save( Archive & ar ) const
{
  const std::vector<std::vector<double>> mat =
  {
    { rotation_( 0, 0 ), rotation_( 0, 1 ), rotation_( 0, 2 ) },
    { rotation_( 1, 0 ), rotation_( 1, 1 ), rotation_( 1, 2 ) },
    { rotation_( 2, 0 ), rotation_( 2, 1 ), rotation_( 2, 2 ) }
  };

  ar( cereal::make_nvp( "rotation", mat ) );

  const std::vector<double> vec = { center_( 0 ), center_( 1 ), center_( 2 ) };
  ar( cereal::make_nvp( "center", vec ) );
}

template <class Archive>
void Pose3::load( Archive & ar )
{
  std::vector<std::vector<double>> mat( 3, std::vector<double>( 3 ) );
  ar( cereal::make_nvp( "rotation", mat ) );
  // copy back to the rotation
  rotation_.row( 0 ) = Eigen::Map<const Vec3>( &( mat[0][0] ) );
  rotation_.row( 1 ) = Eigen::Map<const Vec3>( &( mat[1][0] ) );
  rotation_.row( 2 ) = Eigen::Map<const Vec3>( &( mat[2][0] ) );

  std::vector<double> vec( 3 );
  ar( cereal::make_nvp( "center", vec ) );
  center_ = Eigen::Map<const Vec3>( &vec[0] );
}

template void Pose3::load
<cereal::JSONInputArchive>
( cereal::JSONInputArchive & ar );

template void Pose3::save
<cereal::JSONOutputArchive>
(cereal::JSONOutputArchive & ar) const;

template void Pose3::load
<cereal::XMLInputArchive>
( cereal::XMLInputArchive & ar );

template void Pose3::save
<cereal::XMLOutputArchive>
(cereal::XMLOutputArchive & ar) const;

template void Pose3::load
<cereal::BinaryInputArchive>
( cereal::BinaryInputArchive & ar );

template void Pose3::save
<cereal::BinaryOutputArchive>
(cereal::BinaryOutputArchive & ar) const;

template void Pose3::load
<cereal::PortableBinaryInputArchive>
( cereal::PortableBinaryInputArchive & ar );

template void Pose3::save
<cereal::PortableBinaryOutputArchive>
(cereal::PortableBinaryOutputArchive & ar) const;
    
} // namespace geometry
} // namespace openMVG


