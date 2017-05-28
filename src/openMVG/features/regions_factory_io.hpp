// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_REGIONS_FACTORY_IO_HPP
#define OPENMVG_FEATURES_REGIONS_FACTORY_IO_HPP

#include "openMVG/features/regions_factory.hpp"
//--
// Register region type for serialization
//--
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_Regions, "SIFT_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::SIFT_Regions)
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Float_Regions, "AKAZE_Float_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::AKAZE_Float_Regions)
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Liop_Regions, "AKAZE_Liop_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::AKAZE_Liop_Regions)
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Binary_Regions, "AKAZE_Binary_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::AKAZE_Binary_Regions)

#endif // OPENMVG_FEATURES_REGIONS_FACTORY_IO_HPP