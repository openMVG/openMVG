// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/multiview/projection.hpp"

#include <cereal/types/vector.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG
{
namespace cameras
{

Pinhole_Intrinsic::Pinhole_Intrinsic(
  unsigned int w,
  unsigned int h,
  double focal_length_pix,
  double ppx,
  double ppy )
  : IntrinsicBase( w, h )
{
  K_ << focal_length_pix, 0., ppx, 0., focal_length_pix, ppy, 0., 0., 1.;
  Kinv_ = K_.inverse();
}

Pinhole_Intrinsic::Pinhole_Intrinsic(
  unsigned int w,
  unsigned int h,
  const Mat3& K)
  : IntrinsicBase( w, h ), K_(K)
{
  K_(0,0) = K_(1,1) = (K(0,0) + K(1,1)) / 2.0;
  Kinv_ = K_.inverse();
}

EINTRINSIC Pinhole_Intrinsic::getType() const
{
  return PINHOLE_CAMERA;
}

const Mat3& Pinhole_Intrinsic::K() const
{
  return K_;
}

const Mat3& Pinhole_Intrinsic::Kinv() const
{
  return Kinv_;
}

double Pinhole_Intrinsic::focal() const
{
  return K_( 0, 0 );
}

Vec2 Pinhole_Intrinsic::principal_point() const
{
  return Vec2( K_( 0, 2 ), K_( 1, 2 ) );
}

Vec3 Pinhole_Intrinsic::operator () ( const Vec2& p ) const
{
  return (Kinv_ * p.homogeneous()).normalized();
}

Vec2 Pinhole_Intrinsic::cam2ima( const Vec2& p ) const
{
  return focal() * p + principal_point();
}

Vec2 Pinhole_Intrinsic::ima2cam( const Vec2& p ) const
{
  return ( p -  principal_point() ) / focal();
}

bool Pinhole_Intrinsic::have_disto() const
{
  return false;
}

Vec2 Pinhole_Intrinsic::add_disto( const Vec2& p ) const
{
  return p;
}

Vec2 Pinhole_Intrinsic::remove_disto( const Vec2& p ) const
{
  return p;
}

double Pinhole_Intrinsic::imagePlane_toCameraPlaneError( double value ) const
{
  return value / focal();
}

Mat34 Pinhole_Intrinsic::get_projective_equivalent( const geometry::Pose3 & pose ) const
{
  Mat34 P;
  P_From_KRt( K(), pose.rotation(), pose.translation(), &P );
  return P;
}

std::vector<double> Pinhole_Intrinsic::getParams() const
{
  return  { K_(0, 0), K_(0, 2), K_(1, 2) };
}

bool Pinhole_Intrinsic::updateFromParams(const std::vector<double> & params)
{
  if ( params.size() == 3 )
  {
    *this = Pinhole_Intrinsic( w_, h_, params[0], params[1], params[2] );
    return true;
  }
  else
  {
    return false;
  }
}

std::vector<int> Pinhole_Intrinsic::subsetParameterization
(
  const Intrinsic_Parameter_Type & parametrization) const
{
  std::vector<int> constant_index;
  const int param = static_cast<int>(parametrization);
  if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH)
       || param & (int)Intrinsic_Parameter_Type::NONE )
  {
    constant_index.push_back(0);
  }
  if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT)
      || param & (int)Intrinsic_Parameter_Type::NONE )
  {
    constant_index.push_back(1);
    constant_index.push_back(2);
  }
  return constant_index;
}

Vec2 Pinhole_Intrinsic::get_ud_pixel( const Vec2& p ) const
{
  return p;
}

Vec2 Pinhole_Intrinsic::get_d_pixel( const Vec2& p ) const
{
  return p;
}

template <class Archive>
void Pinhole_Intrinsic::save( Archive & ar ) const
{
  IntrinsicBase::save(ar);
  ar( cereal::make_nvp( "focal_length", K_( 0, 0 ) ) );
  const std::vector<double> pp = {K_( 0, 2 ), K_( 1, 2 )};
  ar( cereal::make_nvp( "principal_point", pp ) );
}


template <class Archive>
void Pinhole_Intrinsic::load( Archive & ar )
{
  IntrinsicBase::load(ar);
  double focal_length;
  ar( cereal::make_nvp( "focal_length", focal_length ) );
  std::vector<double> pp( 2 );
  ar( cereal::make_nvp( "principal_point", pp ) );
  *this = Pinhole_Intrinsic( w_, h_, focal_length, pp[0], pp[1] );
}

IntrinsicBase * Pinhole_Intrinsic::clone( void ) const
{
  return new class_type( *this );
}

template void Pinhole_Intrinsic::load<cereal::JSONInputArchive>(class cereal::JSONInputArchive &);
template void Pinhole_Intrinsic::save<cereal::JSONOutputArchive>(class cereal::JSONOutputArchive &)const;
template void Pinhole_Intrinsic::load<cereal::XMLInputArchive>(class cereal::XMLInputArchive &);
template void Pinhole_Intrinsic::save<cereal::XMLOutputArchive>(class cereal::XMLOutputArchive &)const;
template void Pinhole_Intrinsic::load<cereal::PortableBinaryInputArchive>(class cereal::PortableBinaryInputArchive &);
template void Pinhole_Intrinsic::save<cereal::PortableBinaryOutputArchive>(class cereal::PortableBinaryOutputArchive &)const;

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic, "pinhole");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic);