// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_LANDMARK_IO_HPP
#define OPENMVG_SFM_SFM_LANDMARK_IO_HPP

#include "openMVG/sfm/sfm_landmark.hpp"

#include <cereal/cereal.hpp> // Serialization

#include <vector>

template <class Archive>
void openMVG::sfm::Observation::save( Archive & ar) const
{
  ar(cereal::make_nvp("id_feat", id_feat ));
  const std::vector<double> pp { x(0), x(1) };
  ar(cereal::make_nvp("x", pp));
}

template <class Archive>
void openMVG::sfm::Observation::load( Archive & ar)
{
  ar(cereal::make_nvp("id_feat", id_feat ));
  std::vector<double> p(2);
  ar(cereal::make_nvp("x", p));
  x = Eigen::Map<const Vec2>(&p[0]);
}


template <class Archive>
void openMVG::sfm::Landmark::save( Archive & ar) const
{
  const std::vector<double> point { X(0), X(1), X(2) };
  ar(cereal::make_nvp("X", point ));
  ar(cereal::make_nvp("observations", obs));
}


template <class Archive>
void openMVG::sfm::Landmark::load( Archive & ar)
{
  std::vector<double> point(3);
  ar(cereal::make_nvp("X", point ));
  X = Eigen::Map<const Vec3>(&point[0]);
  ar(cereal::make_nvp("observations", obs));
}

#endif // OPENMVG_SFM_SFM_LANDMARK_IO_HPP
