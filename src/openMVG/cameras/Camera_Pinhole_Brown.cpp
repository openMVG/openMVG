// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Sida Li, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Brown.hpp"

#include <cereal/types/vector.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG
{
namespace cameras
{

Pinhole_Intrinsic_Brown_T2::Pinhole_Intrinsic_Brown_T2(
  int w,  int h,
  double focal, double ppx, double ppy,
  double k1, double k2, double k3,
  double t1, double t2 )
  : Pinhole_Intrinsic( w, h, focal, ppx, ppy )
{
  params_ = {k1, k2, k3, t1, t2};
}

EINTRINSIC Pinhole_Intrinsic_Brown_T2::getType() const
{
  return PINHOLE_CAMERA_BROWN;
}

bool Pinhole_Intrinsic_Brown_T2::have_disto() const
{
  return true;
}

Vec2 Pinhole_Intrinsic_Brown_T2::add_disto( const Vec2 & p ) const
{
  return ( p + distoFunction( params_, p ) );
}

Vec2 Pinhole_Intrinsic_Brown_T2::remove_disto( const Vec2 & p ) const
{
  const double epsilon = 1e-10; //criteria to stop the iteration
  Vec2 p_u = p;

  Vec2 d = distoFunction(params_, p_u);
  while ((p_u + d - p).lpNorm<1>() > epsilon) //manhattan distance between the two points
  {
    p_u = p - d;
    d = distoFunction(params_, p_u);
  }

  return p_u;
}

std::vector<double> Pinhole_Intrinsic_Brown_T2::getParams() const
{
  std::vector<double> params = Pinhole_Intrinsic::getParams();
  params.insert(params.end(), std::begin(params_), std::end(params_));
  return params;
}

bool Pinhole_Intrinsic_Brown_T2::updateFromParams( const std::vector<double> & params )
{
  if ( params.size() == 8 )
  {
    *this = Pinhole_Intrinsic_Brown_T2(
              w_, h_,
              params[0], params[1], params[2], // focal, ppx, ppy
              params[3], params[4], params[5], // K1, K2, K3
              params[6], params[7] );          // T1, T2
    return true;
  }
  else
  {
    return false;
  }
}

std::vector<int> Pinhole_Intrinsic_Brown_T2::subsetParameterization
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
  if ( !(param & (int)Intrinsic_Parameter_Type::ADJUST_DISTORTION)
      || param & (int)Intrinsic_Parameter_Type::NONE )
  {
    constant_index.push_back(3);
    constant_index.push_back(4);
    constant_index.push_back(5);
    constant_index.push_back(6);
    constant_index.push_back(7);
  }
  return constant_index;
}

Vec2 Pinhole_Intrinsic_Brown_T2::get_ud_pixel( const Vec2& p ) const
{
  return cam2ima( remove_disto( ima2cam( p ) ) );
}

Vec2 Pinhole_Intrinsic_Brown_T2::get_d_pixel( const Vec2& p ) const
{
  return cam2ima( add_disto( ima2cam( p ) ) );
}

template <class Archive>
void Pinhole_Intrinsic_Brown_T2::save( Archive & ar ) const
{
  ar(cereal::base_class<Pinhole_Intrinsic>(this));
  ar( cereal::make_nvp( "disto_t2", params_ ) );
}

template <class Archive>
void Pinhole_Intrinsic_Brown_T2::load( Archive & ar )
{
  ar(cereal::base_class<Pinhole_Intrinsic>(this));
  ar( cereal::make_nvp( "disto_t2", params_ ) );
}

IntrinsicBase * Pinhole_Intrinsic_Brown_T2::clone( void ) const
{
  return new class_type( *this );
}

Vec2 Pinhole_Intrinsic_Brown_T2::distoFunction( const std::vector<double> & params, const Vec2 & p )
{
  const double k1 = params[0], k2 = params[1], k3 = params[2], t1 = params[3], t2 = params[4];
  const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
  const double r4 = r2 * r2;
  const double r6 = r4 * r2;
  const double k_diff = ( k1 * r2 + k2 * r4 + k3 * r6 );
  const double t_x = t2 * ( r2 + 2 * p( 0 ) * p( 0 ) ) + 2 * t1 * p( 0 ) * p( 1 );
  const double t_y = t1 * ( r2 + 2 * p( 1 ) * p( 1 ) ) + 2 * t2 * p( 0 ) * p( 1 );
  Vec2 d( p( 0 ) * k_diff + t_x, p( 1 ) * k_diff + t_y );
  return d;
}

template void Pinhole_Intrinsic_Brown_T2::load<cereal::JSONInputArchive>(class cereal::JSONInputArchive &);
template void Pinhole_Intrinsic_Brown_T2::save<cereal::JSONOutputArchive>(class cereal::JSONOutputArchive &)const;
template void Pinhole_Intrinsic_Brown_T2::load<cereal::XMLInputArchive>(class cereal::XMLInputArchive &);
template void Pinhole_Intrinsic_Brown_T2::save<cereal::XMLOutputArchive>(class cereal::XMLOutputArchive &)const;
template void Pinhole_Intrinsic_Brown_T2::load<cereal::PortableBinaryInputArchive>(class cereal::PortableBinaryInputArchive &);
template void Pinhole_Intrinsic_Brown_T2::save<cereal::PortableBinaryOutputArchive>(class cereal::PortableBinaryOutputArchive &)const;

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME( openMVG::cameras::Pinhole_Intrinsic_Brown_T2, "pinhole_brown_t2" );
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Brown_T2)
