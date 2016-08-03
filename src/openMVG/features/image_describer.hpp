
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

#include <string>
#include <iostream>

namespace openMVG {
namespace features {

/**
 * @brief The preset to control the number of detected regions
 */
enum EDESCRIBER_PRESET
{
  LOW_PRESET = 0,
  MEDIUM_PRESET,
  NORMAL_PRESET,
  HIGH_PRESET,
  ULTRA_PRESET
};

/**
 * @brief It returns the preset from a string.
 * @param[in] sPreset the input string.
 * @return the associated describer preset.
 */
EDESCRIBER_PRESET describerPreset_stringToEnum(const std::string& sPreset);

/**
 * @brief It converts a preset to a string.
 * @param[in] preset the describer preset enum to convert.
 * @return the string associated to the describer preset.
 */
std::string describerPreset_enumToString(const EDESCRIBER_PRESET preset);

/**
 * @brief It write a describer preset into a stream by converting it to a string. 
 * @param[in] os the stream where to write the preset.
 * @param[in] p the preset to write.
 * @return the modified stream.
 */
std::ostream& operator<<(std::ostream& os, EDESCRIBER_PRESET p);

/**
 * @brief It read a describer preset from a stream. 
 * @param[in] in the stream from which the preset is read.
 * @param[out] p the preset read from the stream.
 * @return the modified stream without the read preset.
 */
std::istream& operator>>(std::istream& in, EDESCRIBER_PRESET& p);


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
    const image::Image<unsigned char> * mask = nullptr) = 0;

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
