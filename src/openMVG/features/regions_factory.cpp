// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/regions_factory.hpp"

//--
// Register region type for serialization
//--
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/json.hpp>

CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_Regions, "SIFT_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::SIFT_Regions)
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Float_Regions, "AKAZE_Float_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::AKAZE_Float_Regions)
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Liop_Regions, "AKAZE_Liop_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::AKAZE_Liop_Regions)
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Binary_Regions, "AKAZE_Binary_Regions");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Regions, openMVG::features::AKAZE_Binary_Regions)


bool openMVG::features::registerRegionFormats()
{
  // Force dynamic inititalization of region formats.
  // See http://uscilab.github.io/cereal/polymorphism.html#registering-from-a-source-file
  //
  // "Be careful that there are special considerations to make for placing registration in 
  // a source file, especially if you will not be explicitly referencing anything within
  // that source file."

    openMVG::features::SIFT_Regions r1;
    openMVG::features::AKAZE_Float_Regions r2;
    openMVG::features::AKAZE_Liop_Regions r3;
    openMVG::features::AKAZE_Binary_Regions r4;
    return true;
}