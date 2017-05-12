// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_IO_HPP
#define OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_IO_HPP

#include "nonFree/sift/SIFT_describer.hpp"

#include <cereal/types/polymorphic.hpp>

template<class Archive>
void openMVG::features::SIFT_Image_describer::Params::serialize( Archive & ar )
{
  ar(
    cereal::make_nvp("first_octave", _first_octave),
    cereal::make_nvp("num_octaves",_num_octaves),
    cereal::make_nvp("num_scales",_num_scales),
    cereal::make_nvp("edge_threshold",_edge_threshold),
    cereal::make_nvp("peak_threshold",_peak_threshold),
    cereal::make_nvp("root_sift",_root_sift));
}

template<class Archive>
void openMVG::features::SIFT_Image_describer::serialize( Archive & ar )
{
  ar(
   cereal::make_nvp("params", _params),
   cereal::make_nvp("bOrientation", _bOrientation));
}

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_Image_describer, "SIFT_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::SIFT_Image_describer)

#endif // OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_IO_HPP
