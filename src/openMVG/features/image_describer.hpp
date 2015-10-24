
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_IMAGE_DESCRIBER_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/features/regions.hpp"
#include "openMVG/image/image_container.hpp"
#include <memory>
#include <cereal/cereal.hpp> // Serialization

#include <exception>
#include <string>

namespace openMVG {
namespace features {

enum EDESCRIBER_PRESET
{
  LOW_PRESET = 0,
  MEDIUM_PRESET,
  NORMAL_PRESET,
  HIGH_PRESET,
  ULTRA_PRESET
};

inline EDESCRIBER_PRESET describerPreset_stringToEnum(const std::string& sPreset)
{
  if(sPreset == "LOW")
    return LOW_PRESET;
  if (sPreset == "MEDIUM")
    return MEDIUM_PRESET;
  if(sPreset == "NORMAL")
    return NORMAL_PRESET;
  if (sPreset == "HIGH")
    return HIGH_PRESET;
  if (sPreset == "ULTRA")
    return ULTRA_PRESET;
  return NORMAL_PRESET;
}

inline std::string describerPreset_enumToString(const EDESCRIBER_PRESET preset)
{
  if(preset == LOW_PRESET)
    return "LOW";
  if (preset == MEDIUM_PRESET)
    return "MEDIUM";
  if(preset == NORMAL_PRESET)
    return "NORMAL";
  if (preset == HIGH_PRESET)
    return "HIGH";
  if (preset == ULTRA_PRESET)
    return "ULTRA";
  throw std::invalid_argument("Unrecognized EDESCRIBER_PRESET "+std::to_string(preset));
}


/// A pure virtual class for image description computation
class Image_describer
{
public:
  Image_describer() {}
  virtual ~Image_describer() {}

  /**
  @brief Use a preset to control the number of detected regions
  @param preset The preset configuration
  @return True if configuration succeed.
  */
  virtual bool Set_configuration_preset(EDESCRIBER_PRESET preset) = 0;

  bool Set_configuration_preset(const std::string& preset)
  {
    return Set_configuration_preset(describerPreset_stringToEnum(preset));
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  virtual bool Describe(const image::Image<unsigned char> & image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = NULL) = 0;

  /// Allocate regions depending of the Image_describer
  virtual void Allocate(std::unique_ptr<Regions> &regions) const = 0;

  //--
  // IO - one file for region features, one file for region descriptors
  //--

  virtual bool Load(Regions * regions,
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return regions->Load(sfileNameFeats, sfileNameDescs);
  }

  virtual bool Save(const Regions * regions,
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return regions->Save(sfileNameFeats, sfileNameDescs);
  };

  virtual bool LoadFeatures(Regions * regions,
    const std::string& sfileNameFeats) const
  {
    return regions->LoadFeatures(sfileNameFeats);
  }
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_IMAGE_DESCRIBER_HPP
