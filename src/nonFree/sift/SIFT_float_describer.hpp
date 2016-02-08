/* 
 * File:   SIFT_Float_describer.hpp
 * Author: sgaspari
 *
 * Created on September 23, 2015, 5:52 PM
 */

#pragma once

#include "SIFT_describer.hpp"
#include <openMVG/features/descriptor.hpp>
#include <openMVG/features/image_describer.hpp>
#include <openMVG/features/regions_factory.hpp>
#include <nonFree/sift/sift.hpp>

#include <cereal/cereal.hpp>

#include <iostream>
#include <numeric>

extern "C" {
#include "nonFree/sift/vl/sift.h"
}

namespace openMVG {
namespace features {

class SIFT_float_describer : public Image_describer
{
public:
  SIFT_float_describer(const SiftParams & params = SiftParams(), bool bOrientation = true)
    :Image_describer(), _params(params), _bOrientation(bOrientation) {}

  ~SIFT_float_describer() {}

  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
    return _params.setPreset(preset);
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = NULL)
  {
    return extractSIFT<float>(image, regions, _params, _bOrientation, mask);
  }

  /*/// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
      regions.reset( new SIFT_Float_Regions );
  }*/

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(
     cereal::make_nvp("params", _params),
     cereal::make_nvp("bOrientation", _bOrientation));
  }
  
  void Allocate(std::unique_ptr<Regions>& regions) const
  {
    regions.reset(new SIFT_Float_Regions);
  }
  
private:
  SiftParams _params;
  bool _bOrientation;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_float_describer, "SIFT_float_describer");




