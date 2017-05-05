// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/numeric/numeric.h"

#include <cereal/types/vector.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG{
namespace cameras{
namespace radial_distortion{


template <class Disto_Functor>
double bisection_Radius_Solve(
  const std::vector<double> & params, // radial distortion parameters
  double r2, // targeted radius
  Disto_Functor & functor,
  double epsilon = 1e-10 // criteria to stop the bisection
)
{
  // Guess plausible upper and lower bound
  double lowerbound = r2, upbound = r2;
  while ( functor( params, lowerbound ) > r2 )
  {
    lowerbound /= 1.05;
  }
  while ( functor( params, upbound ) < r2 )
  {
    upbound *= 1.05;
  }

  // Perform a bisection until epsilon accuracy is not reached
  while ( epsilon < upbound - lowerbound )
  {
    const double mid = .5 * ( lowerbound + upbound );
    if ( functor( params, mid ) > r2 )
    {
      upbound = mid;
    }
    else
    {
      lowerbound = mid;
    }
  }
  return .5 * ( lowerbound + upbound );
}

} // namespace radial_distortion

Pinhole_Intrinsic_Radial_K1::Pinhole_Intrinsic_Radial_K1(
    int w, int h,
    double focal, double ppx, double ppy,
    double k1 )
    : Pinhole_Intrinsic( w, h, focal, ppx, ppy )
{
    params_ = {k1};
}

EINTRINSIC Pinhole_Intrinsic_Radial_K1::Pinhole_Intrinsic_Radial_K1::getType() const
{
    return PINHOLE_CAMERA_RADIAL1;
}

bool Pinhole_Intrinsic_Radial_K1::have_disto() const
{
    return true;
}

Vec2 Pinhole_Intrinsic_Radial_K1::add_disto( const Vec2 & p ) const
{

    const double k1 = params_[0];

    const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
    const double r_coeff = ( 1. + k1 * r2 );

    return ( p * r_coeff );
}

Vec2 Pinhole_Intrinsic_Radial_K1::remove_disto( const Vec2& p ) const
{
    // Compute the radius from which the point p comes from thanks to a bisection
    // Minimize disto(radius(p')^2) == actual Squared(radius(p))

    const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
    const double radius = ( r2 == 0 ) ?
                        1. :
                        ::sqrt( radial_distortion::bisection_Radius_Solve( params_, r2, distoFunctor ) / r2 );
    return radius * p;
}

std::vector<double> Pinhole_Intrinsic_Radial_K1::getParams() const
{
    std::vector<double> params = Pinhole_Intrinsic::getParams();
    params.insert(params.end(), std::begin(params_), std::end(params_));
    return params;
}

bool Pinhole_Intrinsic_Radial_K1::updateFromParams( const std::vector<double> & params )
{
    if ( params.size() == 4 )
    {
    *this = Pinhole_Intrinsic_Radial_K1(
                w_, h_,
                params[0], params[1], params[2], // focal, ppx, ppy
                params[3] ); //K1
    return true;
    }
    else
    {
    return false;
    }
}

std::vector<int> Pinhole_Intrinsic_Radial_K1::subsetParameterization
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
    }
    return constant_index;
}

Vec2 Pinhole_Intrinsic_Radial_K1::get_ud_pixel( const Vec2& p ) const
{
    return cam2ima( remove_disto( ima2cam( p ) ) );
}

Vec2 Pinhole_Intrinsic_Radial_K1::get_d_pixel( const Vec2& p ) const
{
    return cam2ima( add_disto( ima2cam( p ) ) );
}

template <class Archive>
void Pinhole_Intrinsic_Radial_K1::save( Archive & ar ) const
{
    Pinhole_Intrinsic::save( ar );
    ar( cereal::make_nvp( "disto_k1", params_ ) );
}

template <class Archive>
void Pinhole_Intrinsic_Radial_K1::load( Archive & ar )
{
    Pinhole_Intrinsic::load(ar);
    ar( cereal::make_nvp( "disto_k1", params_ ) );
}

IntrinsicBase * Pinhole_Intrinsic_Radial_K1::clone( void ) const
{
    return new class_type( *this );
}


double Pinhole_Intrinsic_Radial_K1::distoFunctor( const std::vector<double> & params, double r2 )
{
    const double & k1 = params[0];
    return r2 * Square( 1. + r2 * k1 );
}

Pinhole_Intrinsic_Radial_K3::Pinhole_Intrinsic_Radial_K3(
    int w, int h,
    double focal, double ppx, double ppy,
    double k1, double k2, double k3 )
    : Pinhole_Intrinsic( w, h, focal, ppx, ppy )
{
    params_ = {k1, k2, k3};
}

EINTRINSIC Pinhole_Intrinsic_Radial_K3::getType() const
{
    return PINHOLE_CAMERA_RADIAL3;
}

bool Pinhole_Intrinsic_Radial_K3::have_disto() const
{
    return true;
}

Vec2 Pinhole_Intrinsic_Radial_K3::add_disto( const Vec2 & p ) const
{
    const double & k1 = params_[0], & k2 = params_[1], & k3 = params_[2];

    const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
    const double r4 = r2 * r2;
    const double r6 = r4 * r2;
    const double r_coeff = ( 1. + k1 * r2 + k2 * r4 + k3 * r6 );

    return ( p * r_coeff );
}
Vec2 Pinhole_Intrinsic_Radial_K3::remove_disto( const Vec2& p ) const
{
    // Compute the radius from which the point p comes from thanks to a bisection
    // Minimize disto(radius(p')^2) == actual Squared(radius(p))

    const double r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
    const double radius = ( r2 == 0 ) ? //1. : ::sqrt(bisectionSolve(_params, r2) / r2);
                        1. :
                        ::sqrt( radial_distortion::bisection_Radius_Solve( params_, r2, distoFunctor ) / r2 );
    return radius * p;
}

std::vector<double> Pinhole_Intrinsic_Radial_K3::getParams() const
{
    std::vector<double> params = Pinhole_Intrinsic::getParams();
    params.insert( params.end(), std::begin(params_), std::end(params_));
    return params;
}

bool Pinhole_Intrinsic_Radial_K3::updateFromParams( const std::vector<double> & params )
{
    if ( params.size() == 6 )
    {
    *this = Pinhole_Intrinsic_Radial_K3(
                w_, h_,
                params[0], params[1], params[2], // focal, ppx, ppy
                params[3], params[4], params[5] ); // K1, K2, K3
    return true;
    }
    else
    {
    return false;
    }
}

std::vector<int> Pinhole_Intrinsic_Radial_K3::subsetParameterization
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
    }
    return constant_index;
}

Vec2 Pinhole_Intrinsic_Radial_K3::get_ud_pixel( const Vec2& p ) const
{
    return cam2ima( remove_disto( ima2cam( p ) ) );
}

Vec2 Pinhole_Intrinsic_Radial_K3::get_d_pixel( const Vec2& p ) const
{
    return cam2ima( add_disto( ima2cam( p ) ) );
}

template <class Archive>
void Pinhole_Intrinsic_Radial_K3::save( Archive & ar ) const
{
    Pinhole_Intrinsic::save(ar);
    ar( cereal::make_nvp( "disto_k3", params_ ) );
}

template <class Archive>
void Pinhole_Intrinsic_Radial_K3::load( Archive & ar )
{
    Pinhole_Intrinsic::load(ar);
    ar( cereal::make_nvp( "disto_k3", params_ ) );
}

IntrinsicBase * Pinhole_Intrinsic_Radial_K3::clone( void ) const
{
    return new class_type( *this );
}

double Pinhole_Intrinsic_Radial_K3::distoFunctor( const std::vector<double> & params, double r2 )
{
    const double & k1 = params[0], & k2 = params[1], & k3 = params[2];
    return r2 * Square( 1. + r2 * ( k1 + r2 * ( k2 + r2 * k3 ) ) );
}

template void Pinhole_Intrinsic_Radial_K1::load<cereal::JSONInputArchive>(class cereal::JSONInputArchive &);
template void Pinhole_Intrinsic_Radial_K1::save<cereal::JSONOutputArchive>(class cereal::JSONOutputArchive &)const;

template void Pinhole_Intrinsic_Radial_K3::load<cereal::JSONInputArchive>(class cereal::JSONInputArchive &);
template void Pinhole_Intrinsic_Radial_K3::save<cereal::JSONOutputArchive>(class cereal::JSONOutputArchive &)const;

template void Pinhole_Intrinsic_Radial_K1::load<cereal::XMLInputArchive>(class cereal::XMLInputArchive &);
template void Pinhole_Intrinsic_Radial_K1::save<cereal::XMLOutputArchive>(class cereal::XMLOutputArchive &)const;

template void Pinhole_Intrinsic_Radial_K3::load<cereal::XMLInputArchive>(class cereal::XMLInputArchive &);
template void Pinhole_Intrinsic_Radial_K3::save<cereal::XMLOutputArchive>(class cereal::XMLOutputArchive &)const;

template void Pinhole_Intrinsic_Radial_K1::load<cereal::PortableBinaryInputArchive>(class cereal::PortableBinaryInputArchive &);
template void Pinhole_Intrinsic_Radial_K1::save<cereal::PortableBinaryOutputArchive>(class cereal::PortableBinaryOutputArchive &)const;

template void Pinhole_Intrinsic_Radial_K3::load<cereal::PortableBinaryInputArchive>(class cereal::PortableBinaryInputArchive &);
template void Pinhole_Intrinsic_Radial_K3::save<cereal::PortableBinaryOutputArchive>(class cereal::PortableBinaryOutputArchive &)const;

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_Radial_K1, "pinhole_radial_k1");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Radial_K1);
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Pinhole_Intrinsic_Radial_K3, "pinhole_radial_k3");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Pinhole_Intrinsic_Radial_K3);
