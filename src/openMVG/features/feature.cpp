// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/feature.hpp"

#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

namespace openMVG {
namespace features {

PointFeature::PointFeature(float x, float y): coords_(x, y) {}

float PointFeature::x() const { return coords_(0); }
float PointFeature::y() const { return coords_(1); }
const Vec2f & PointFeature::coords() const { return coords_;}

float& PointFeature::x() { return coords_(0); }
float& PointFeature::y() { return coords_(1); }
Vec2f& PointFeature::coords() { return coords_;}


//with overloaded operators:
std::ostream& operator<<(std::ostream& out, const PointFeature& obj)
{
  return out << obj.coords_(0) << " " << obj.coords_(1);
}

std::istream& operator>>(std::istream& in, PointFeature& obj)
{
  return in >> obj.coords_(0) >> obj.coords_(1);
}


SIOPointFeature::SIOPointFeature
(
  float x,
  float y,
  float scale,
  float orient
):
  PointFeature(x,y),
  scale_(scale),
  orientation_(orient)
{

}

float SIOPointFeature::scale() const { return scale_; }
float& SIOPointFeature::scale() { return scale_; }
float SIOPointFeature::orientation() const { return orientation_; }
float& SIOPointFeature::orientation() { return orientation_; }

bool SIOPointFeature::operator ==(const SIOPointFeature& b) const {
  return (scale_ == b.scale()) &&
         (orientation_ == b.orientation()) &&
         (x() == b.x()) && (y() == b.y());
};

bool SIOPointFeature::operator !=(const SIOPointFeature& b) const {
  return !((*this)==b);
};

std::ostream& operator<<(std::ostream& out, const SIOPointFeature& obj)
{
  const PointFeature *pf = static_cast<const PointFeature*>(&obj);
  return out << *pf << " " << obj.scale_ << " " << obj.orientation_;
}

std::istream& operator>>(std::istream& in, SIOPointFeature& obj)
{
  PointFeature *pf = static_cast<PointFeature*>(&obj);
  return in >> *pf >> obj.scale_ >> obj.orientation_;
}

/// Return the coterminal angle between [0;2*PI].
/// Angle value must be in Radian.
float getCoterminalAngle(float angle)
{
  const float f2PI = 2.f*M_PI;
  while (angle > f2PI) {
    angle -= f2PI;
  }
  while (angle < 0.0f) {
    angle += f2PI;
  }
  return angle;
}

AffinePointFeature::AffinePointFeature
(
  float x,
  float y,
  float a,
  float b,
  float c
) : PointFeature(x, y), a_(a), b_(b), c_(c)
{
  l1_ = (a + c - std::sqrt(a*a + c*c + 4 * b*b - 2 * a*c)) / 2.f;
  l2_ = (a + c + std::sqrt(a*a + c*c + 4 * b*b - 2 * a*c)) / 2.f;
  l1_ = 1.f / std::sqrt(l1_);
  l2_ = 1.f / std::sqrt(l2_);

  phi_ = 0.f;
  if (b == 0)
  {
    if (a > c)
      phi_ = M_PI / 2; // else 0
  }
  else
  {
    const double t = std::atan(2 * b / (a - c));
    if (a < c)
      phi_ = t / 2;
    else
      phi_ = t / 2 + ((b > 0) ? -M_PI / 2 : M_PI / 2);
  }

  if (l1_ > l2_)
  {
    std::swap(l1_, l2_);
    phi_ = getCoterminalAngle(M_PI / 2 - phi_);
  }
}

float AffinePointFeature::l1() const { return l1_; }
float AffinePointFeature::l2() const { return l2_; }
float AffinePointFeature::orientation() const { return phi_; }

bool AffinePointFeature::operator ==(const AffinePointFeature& b) const {
  return ((x() == b.x()) && (y() == b.y() &&
    (l1_ == b.l1_) && (l2_ == b.l2_) && (phi_ == b.phi_)));
};

bool AffinePointFeature::operator !=(const AffinePointFeature& rhs) const {
  return !((*this) == rhs);
};

float AffinePointFeature::a() const { return a_; }
float AffinePointFeature::b() const { return b_; }
float AffinePointFeature::c() const { return c_; }


std::ostream& operator<<(std::ostream& out, const AffinePointFeature& rhs)
{
  const PointFeature *pf = static_cast<const PointFeature*>(&rhs);
  return out << *pf << " " << rhs.l1_ << " " << rhs.l2_ << " " << rhs.phi_
    << " " << rhs.a_ << " " << rhs.b_ << " " << rhs.c_;
}

std::istream& operator>>(std::istream& in, AffinePointFeature& rhs)
{
  PointFeature *pf = static_cast<PointFeature*>(&rhs);
  return in >> *pf >> rhs.l1_ >> rhs.l2_ >> rhs.phi_ >> rhs.a_ >> rhs.b_ >> rhs.c_;
}

} // namespace features
} // namespace openMVG

