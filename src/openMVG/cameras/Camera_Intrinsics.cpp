// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/stl/hash.hpp"
#include "openMVG/numeric/numeric.h"

#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG
{
namespace cameras
{

Vec2 IntrinsicBase::project(
  const geometry::Pose3 & pose,
  const Vec3 & pt3D ) const
{
  const Vec3 X = pose( pt3D ); // apply pose
  if ( this->have_disto() ) // apply disto & intrinsics
  {
    return this->cam2ima( this->add_disto( X.hnormalized() ) );
  }
  else // apply intrinsics
  {
    return this->cam2ima( X.hnormalized() );
  }
}

Vec2 IntrinsicBase::residual(
  const geometry::Pose3 & pose,
  const Vec3 & X,
  const Vec2 & x ) const
{
  const Vec2 proj = this->project( pose, X );
  return x - proj;
}

template <class Archive>
void IntrinsicBase::save( Archive & ar ) const
{
  ar( cereal::make_nvp( "width", w_ ) );
  ar( cereal::make_nvp( "height", h_ ) );
}

template <class Archive>
void IntrinsicBase::load( Archive & ar )
{
  ar( cereal::make_nvp( "width", w_ ) );
  ar( cereal::make_nvp( "height", h_ ) );
}

std::size_t IntrinsicBase::hashValue() const
{
  size_t seed = 0;
  stl::hash_combine( seed, static_cast<int>( this->getType() ) );
  stl::hash_combine( seed, w_ );
  stl::hash_combine( seed, h_ );
  const std::vector<double> params = this->getParams();
  for ( const auto & param : params )
    stl::hash_combine( seed , param );
  return seed;
}


double AngleBetweenRay(
  const geometry::Pose3 & pose1,
  const IntrinsicBase * intrinsic1,
  const geometry::Pose3 & pose2,
  const IntrinsicBase * intrinsic2,
  const Vec2 & x1, const Vec2 & x2 )
{
  // x = (u, v, 1.0)  // image coordinates
  // X = R.t() * K.inv() * x + C // Camera world point
  // getting the ray:
  // ray = X - C = R.t() * K.inv() * x
  const Vec3 ray1 = ( pose1.rotation().transpose() * intrinsic1->operator()( x1 ) ).normalized();
  const Vec3 ray2 = ( pose2.rotation().transpose() * intrinsic2->operator()( x2 ) ).normalized();
  const double mag = ray1.norm() * ray2.norm();
  const double dotAngle = ray1.dot( ray2 );
  return R2D( acos( clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 ) ) );
}


template void IntrinsicBase::load<class cereal::XMLInputArchive>(class cereal::XMLInputArchive &);
template void IntrinsicBase::save<class cereal::XMLOutputArchive>(class cereal::XMLOutputArchive &) const;
template void IntrinsicBase::load<class cereal::JSONInputArchive>(class cereal::JSONInputArchive &);
template void IntrinsicBase::save<class cereal::JSONOutputArchive>(class cereal::JSONOutputArchive &) const;
template void IntrinsicBase::load<class cereal::BinaryInputArchive>(class cereal::BinaryInputArchive &);
template void IntrinsicBase::save<class cereal::BinaryOutputArchive>(class cereal::BinaryOutputArchive &)const;
template void IntrinsicBase::load<class cereal::PortableBinaryInputArchive>(class cereal::PortableBinaryInputArchive &);
template void IntrinsicBase::save<class cereal::PortableBinaryOutputArchive>(class cereal::PortableBinaryOutputArchive &) const;

} // namespace cameras
} // namespace openMVG

#include "openMVG/cameras/cameras.hpp"

bool openMVG::cameras::IntrinsicBase::registerCameraTypes()
{
  // Force dynamic inititalization of cameras.
  // See http://uscilab.github.io/cereal/polymorphism.html#registering-from-a-source-file
  //
  // "Be careful that there are special considerations to make for placing registration in 
  // a source file, especially if you will not be explicitly referencing anything within
  // that source file."

  openMVG::cameras::Pinhole_Intrinsic p1;
  openMVG::cameras::Pinhole_Intrinsic_Brown_T2 p2;
  openMVG::cameras::Pinhole_Intrinsic_Fisheye p3;
  openMVG::cameras::Pinhole_Intrinsic_Radial_K1 p4;
  openMVG::cameras::Pinhole_Intrinsic_Radial_K3 p5;
  openMVG::cameras::Intrinsic_Spherical p6;
  static_assert(openMVG::cameras::PINHOLE_CAMERA_END==6, 
    "If you've added a camera type, you need to initialize it here.");
  return true;
}
