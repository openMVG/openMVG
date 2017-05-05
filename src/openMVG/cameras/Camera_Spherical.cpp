// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 <Zillow Inc.> Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/Camera_Spherical.hpp"

#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG{
namespace cameras{


Intrinsic_Spherical::Intrinsic_Spherical
(
  unsigned int w,
  unsigned int h
)
: IntrinsicBase(w, h)
{
}

EINTRINSIC Intrinsic_Spherical::getType() const
{
  return CAMERA_SPHERICAL;
}

std::vector<double> Intrinsic_Spherical::getParams() const
{
  return {};
}

bool Intrinsic_Spherical::updateFromParams(const std::vector<double> &params)
{
  return true;
}

std::vector<int> Intrinsic_Spherical::subsetParameterization
(
  const Intrinsic_Parameter_Type & parametrization
) const
{
  return {};
}

Vec2 Intrinsic_Spherical::cam2ima(const Vec2 &p) const
{
  const size_t size = std::max(w(), h());
  return {
    p.x() * size + w() / 2.0,
    p.y() * size + h() / 2.0 };
}

Vec2 Intrinsic_Spherical::ima2cam(const Vec2 &p) const
{
  const size_t size = std::max(w(), h());
  return {
    (p.x() - w() / 2.0) / size,
    (p.y() - h() / 2.0) / size };
}

Vec3 Intrinsic_Spherical::operator () ( const Vec2& p ) const
{
  const Vec2 uv = ima2cam(p);

  const double
    lon = uv.x() * 2 * M_PI,
    lat = uv.y() * 2 * M_PI;

  return {
    cos(lat) * sin(lon),
    -sin(lat),
    cos(lat) * cos(lon)};
}

Vec2 Intrinsic_Spherical::project(
  const geometry::Pose3 & pose,
  const Vec3 & pt3D ) const
{
  const Vec3 X = pose( pt3D ); // apply pose
  const double lon = atan2(X.x(), X.z()); // Horizontal normalization of the  X-Z component
  const double lat = atan2(-X.y(), sqrt(X.x()*X.x() + X.z()*X.z())); // Tilt angle
  // denormalization (angle to pixel value)
  return cam2ima({lon / (2 * M_PI), lat / (2 * M_PI)});
}

bool Intrinsic_Spherical::have_disto() const 
{
    return false; 
}

Vec2 Intrinsic_Spherical::add_disto(const Vec2 &p) const
{ 
    return p; 
}

Vec2 Intrinsic_Spherical::remove_disto(const Vec2 &p) const
{ 
    return p; 
}

Vec2 Intrinsic_Spherical::get_ud_pixel(const Vec2 &p) const
{ 
    return p; 
}

Vec2 Intrinsic_Spherical::get_d_pixel(const Vec2 &p) const
{ 
    return p; 
}

double Intrinsic_Spherical::imagePlane_toCameraPlaneError(double value) const
{
    return value; 
}

Mat34 Intrinsic_Spherical::get_projective_equivalent(const geometry::Pose3 &pose) const
{
  return HStack(pose.rotation(), pose.translation());
}

template <class Archive>
void Intrinsic_Spherical::save( Archive & ar ) const
{
  ar(cereal::base_class<IntrinsicBase>(this));
}

template <class Archive>
void Intrinsic_Spherical::load( Archive & ar )
{
  ar(cereal::base_class<IntrinsicBase>(this));
}

IntrinsicBase * Intrinsic_Spherical::clone( void ) const
{
  return new class_type( *this );
}

template void Intrinsic_Spherical::load<cereal::JSONInputArchive>(class cereal::JSONInputArchive &);
template void Intrinsic_Spherical::save<cereal::JSONOutputArchive>(class cereal::JSONOutputArchive &)const;
template void Intrinsic_Spherical::load<cereal::XMLInputArchive>(class cereal::XMLInputArchive &);
template void Intrinsic_Spherical::save<cereal::XMLOutputArchive>(class cereal::XMLOutputArchive &)const;
template void Intrinsic_Spherical::load<cereal::PortableBinaryInputArchive>(class cereal::PortableBinaryInputArchive &);
template void Intrinsic_Spherical::save<cereal::PortableBinaryOutputArchive>(class cereal::PortableBinaryOutputArchive &)const;

} // namespace cameras
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::cameras::Intrinsic_Spherical, "spherical");

namespace cereal
{
  // This struct specialization will tell cereal which is the right way to serialize the ambiguity
  template <class Archive> struct specialize<Archive, openMVG::cameras::Intrinsic_Spherical, cereal::specialization::member_load_save> {};
}

CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::cameras::IntrinsicBase, openMVG::cameras::Intrinsic_Spherical);
