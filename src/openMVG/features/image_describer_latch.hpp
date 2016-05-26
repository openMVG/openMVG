// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_LATCH_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_LATCH_IMAGE_DESCRIBER_HPP

#include <iostream>
#include <numeric>

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/features/latch/LatchClassifierOpenMVG.hpp"
#include <cereal/cereal.hpp>

using namespace std;

namespace openMVG {
namespace features {

class LATCH_Image_describer : public Image_describer
{
public:
  LATCH_Image_describer(
  ):Image_describer(),
    latch(LatchClassifierOpenMVG()){
        latch.setImageSize(4000, 3000);
        }

  // Don't need to really define this for the LATCH class yet, until more descriptors come out.
  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
  {
    return true;
  }

  void setImageSize(int width, int height) {
    latch.setImageSize(width, height);
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
    const image::Image<unsigned char> * mask = nullptr ) override
  {
    std::vector<LatchClassifierKeypoint> kpts = latch.identifyFeaturePointsOpenMVG(image.GetMat());


    Allocate(regions);

        // Build alias to cached data
        LATCH_Unsigned_Int_Regions * regionsCasted = dynamic_cast<LATCH_Unsigned_Int_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          LatchClassifierKeypoint ptLatch = kpts[i];

          // Feature masking
          if (mask)
          {
            const image::Image<unsigned char> & maskIma = *mask;
            if (maskIma(ptLatch.y, ptLatch.x) == 0)
              continue;
          }
          regionsCasted->Features()[i] =
            SIOPointFeature(ptLatch.x, ptLatch.y, ptLatch.size, ptLatch.angle);
          // Compute descriptors
          for (int j = 0; j < 16; j++) {
            const unsigned int index = i * 64 + j;
            regionsCasted->Descriptors()[i][j] = static_cast<unsigned int>(latch.getDescriptorSet1()[index]);
          }
        }

    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const override
  {
    return regions.reset(new LATCH_Unsigned_Int_Regions);
  }

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(
     cereal::make_nvp("bOrientation", bOrientation_));
  }


private:
  bool bOrientation_;
  LatchClassifierOpenMVG latch;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::LATCH_Image_describer, "LATCH_Image_describer");

#endif // OPENMVG_FEATURES_LATCH_IMAGE_DESCRIBER_HPP
