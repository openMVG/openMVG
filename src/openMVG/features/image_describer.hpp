// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_IMAGE_DESCRIBER_HPP

#include <memory>
#include <string>

#include "openMVG/features/regions.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG { namespace image { template<typename Type> class Image; } }

namespace openMVG {
namespace features {

enum EDESCRIBER_PRESET
{
  NORMAL_PRESET,
  HIGH_PRESET,
  ULTRA_PRESET
};
/// A pure virtual class for image description computation
class Image_describer
{
public:
  Image_describer() = default;
  virtual ~Image_describer() = default;

  /**
  @brief Use a preset to control the number of detected regions
  @param preset The preset configuration
  @return True if configuration succeed.
  */
  virtual bool Set_configuration_preset
  (
    EDESCRIBER_PRESET preset
  ) = 0;

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe
  (
    const image::Image<unsigned char> & image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = nullptr
  )
  {
    regions = Describe(image, mask);
    return regions != nullptr;
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  @return The detected regions and attributes
  */
  virtual std::unique_ptr<Regions> Describe(
    const image::Image<unsigned char> & image,
    const image::Image<unsigned char> * mask = nullptr) = 0;

  /// Allocate regions depending of the Image_describer
  virtual std::unique_ptr<Regions> Allocate() const = 0;

  //--
  // IO - one file for region features, one file for region descriptors
  //--

  virtual bool Load
  (
    Regions * regions,
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs
  ) const
  {
    return regions->Load(sfileNameFeats, sfileNameDescs);
  }

  virtual bool Save
  (
    const Regions * regions,
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs
  ) const
  {
    return regions->Save(sfileNameFeats, sfileNameDescs);
  };

  virtual bool LoadFeatures
  (
    Regions * regions,
    const std::string& sfileNameFeats
  ) const
  {
    return regions->LoadFeatures(sfileNameFeats);
  }
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_IMAGE_DESCRIBER_HPP
