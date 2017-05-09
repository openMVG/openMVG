// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_landmark.hpp"

#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG {
namespace sfm {

Observation::Observation():id_feat(UndefinedIndexT)
{ 
}

Observation::Observation(const Vec2 & p, IndexT idFeat): x(p), id_feat(idFeat) 
{
}

// Serialization
template <class Archive>
void Observation::save( Archive & ar) const
{
  ar(cereal::make_nvp("id_feat", id_feat ));
  const std::vector<double> pp = { x(0), x(1) };
  ar(cereal::make_nvp("x", pp));
}

// Serialization
template <class Archive>
void Observation::load( Archive & ar)
{
  ar(cereal::make_nvp("id_feat", id_feat ));
  std::vector<double> p(2);
  ar(cereal::make_nvp("x", p));
  x = Eigen::Map<const Vec2>(&p[0]);
}

template void Observation::load
<cereal::JSONInputArchive>
( cereal::JSONInputArchive & ar );

template void Observation::save
<cereal::JSONOutputArchive>
(cereal::JSONOutputArchive & ar) const;

template void Observation::load
<cereal::XMLInputArchive>
( cereal::XMLInputArchive & ar );

template void Observation::save
<cereal::XMLOutputArchive>
(cereal::XMLOutputArchive & ar) const;

template void Observation::load
<cereal::BinaryInputArchive>
( cereal::BinaryInputArchive & ar );

template void Observation::save
<cereal::BinaryOutputArchive>
(cereal::BinaryOutputArchive & ar) const;

template void Observation::load
<cereal::PortableBinaryInputArchive>
(cereal::PortableBinaryInputArchive & ar );

template void Observation::save
<cereal::PortableBinaryOutputArchive>
(cereal::PortableBinaryOutputArchive & ar) const;

// Serialization
template <class Archive>
void Landmark::save( Archive & ar) const
{
  const std::vector<double> point = { X(0), X(1), X(2) };
  ar(cereal::make_nvp("X", point ));
  //ar(cereal::make_nvp("observations", obs));
}

template <class Archive>
void Landmark::load( Archive & ar)
{
  std::vector<double> point(3);
  ar(cereal::make_nvp("X", point ));
  X = Eigen::Map<const Vec3>(&point[0]);
  //ar(cereal::make_nvp("observations", obs));
}

template void Landmark::load
<cereal::JSONInputArchive>
( cereal::JSONInputArchive & ar );

template void Landmark::save
<cereal::JSONOutputArchive>
(cereal::JSONOutputArchive & ar) const;

template void Landmark::load
<cereal::XMLInputArchive>
( cereal::XMLInputArchive & ar );

template void Landmark::save
<cereal::XMLOutputArchive>
(cereal::XMLOutputArchive & ar) const;

template void Landmark::load
<cereal::BinaryInputArchive>
( cereal::BinaryInputArchive & ar );

template void Landmark::save
<cereal::BinaryOutputArchive>
(cereal::BinaryOutputArchive & ar) const;

template void Landmark::load
<cereal::PortableBinaryInputArchive>
(cereal::PortableBinaryInputArchive & ar );

template void Landmark::save
<cereal::PortableBinaryOutputArchive>
(cereal::PortableBinaryOutputArchive & ar) const;

} // namespace sfm
} // namespace openMVG
