// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_REGIONS_FACTORY_HPP
#define OPENMVG_FEATURES_REGIONS_FACTORY_HPP

#include "openMVG/features/binary_regions.hpp"
#include "openMVG/features/scalar_regions.hpp"

namespace openMVG {
namespace features {

/// Define the classic SIFT Keypoint
using SIFT_Regions = Scalar_Regions<SIOPointFeature, unsigned char, 128>;

/// Define the AKAZE Keypoint (with a float descriptor)
using AKAZE_Float_Regions = Scalar_Regions<SIOPointFeature, float, 64>;
/// Define the AKAZE Keypoint (with a LIOP descriptor)
using AKAZE_Liop_Regions = Scalar_Regions<SIOPointFeature, unsigned char, 144>;
/// Define the AKAZE Keypoint (with a binary descriptor saved in an uchar array)
using AKAZE_Binary_Regions = Binary_Regions<SIOPointFeature, 64>;

} // namespace features
} // namespace openMVG

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(openMVG::features::SIFT_Regions)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(openMVG::features::AKAZE_Float_Regions)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(openMVG::features::AKAZE_Liop_Regions)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(openMVG::features::AKAZE_Binary_Regions)

#endif // OPENMVG_FEATURES_REGIONS_FACTORY_HPP
