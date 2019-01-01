// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_IO_HPP
#define OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_IO_HPP

#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/akaze/AKAZE_io.hpp"

#include <cereal/types/polymorphic.hpp>

template<class Archive>
void openMVG::features::AKAZE_Image_describer::Params::serialize(Archive & ar)
{
  ar(options_, eAkazeDescriptor_);
}

template<class Archive>
void openMVG::features::AKAZE_Image_describer::serialize(Archive & ar)
{
  ar(
   cereal::make_nvp("params", params_),
   cereal::make_nvp("bOrientation", bOrientation_));
}


CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Image_describer, "AKAZE_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::AKAZE_Image_describer)

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Image_describer_SURF, "AKAZE_Image_describer_SURF");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::AKAZE_Image_describer_SURF)

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Image_describer_LIOP, "AKAZE_Image_describer_LIOP");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::AKAZE_Image_describer_LIOP)

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Image_describer_MLDB, "AKAZE_Image_describer_MLDB");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::AKAZE_Image_describer_MLDB)

#endif // OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_IO_HPP
