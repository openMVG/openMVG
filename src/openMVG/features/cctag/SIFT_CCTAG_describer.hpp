#ifdef HAVE_CCTAG

#pragma once

#include "CCTAG_describer.hpp"

#include <openMVG/features/image_describer.hpp>
#include <openMVG/features/regions_factory.hpp>

// Sift related
#include <nonFree/sift/SIFT_describer.hpp>

// CCTag related
#include <cctag/ICCTag.hpp>

#include <cereal/cereal.hpp>
#include <iostream>
#include <numeric>

//#include <cereal/cereal.hpp>

#include <iostream>
#include <numeric>

namespace openMVG {
namespace features {

class SIFT_CCTAG_Image_describer : public Image_describer
{
public:
  SIFT_CCTAG_Image_describer(const SiftParams & params = SiftParams(), bool bOrientation = true, std::size_t nRings = 3);

  ~SIFT_CCTAG_Image_describer(){}

  bool Set_configuration_preset(EDESCRIBER_PRESET preset);

  void Set_cctag_use_cuda(bool use_cuda)
  {
    _cctagDescriber.Set_use_cuda(use_cuda);
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
    const image::Image<unsigned char> * mask = NULL);

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
    regions.reset( new SIFT_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    //ar(
    // cereal::make_nvp("params", _paramsSift),
    // cereal::make_nvp("bOrientation", _bOrientation));
  }

private:
  SIFT_Image_describer _siftDescriber;
  CCTAG_Image_describer _cctagDescriber;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_CCTAG_Image_describer, "SIFT_CCTAG_Image_describer");

#endif //HAVE_CCTAG
