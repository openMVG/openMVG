// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_DEEP_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_DEEP_IMAGE_DESCRIBER_HPP

#include <iostream>
#include <limits>
#include <numeric>
#include <cereal/cereal.hpp>

#include "openMVG/features/deep/src/DeepClassifier.hpp"
#include "openMVG/features/deep/src/DeepClassifierTHNets.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/features/tbmr/tbmr.hpp"

#include "opencv2/core.hpp"

using namespace std;

namespace openMVG {
namespace features {

enum EDEEP_DESCRIPTOR
{
  SIAM_DESC_NOTRE_DAME,
  SIAM_DESC_YOSEMITE,
  SIAM_DESC_LIBERTY,

  SIAM_2_STREAM_DESC_NOTRE_DAME,
  SIAM_2_STREAM_DESC_YOSEMITE,
  SIAM_2_STREAM_DESC_LIBERTY,

  PNNET
};

struct DEEPParams
{
  // Because of the absurdities of C++ not having an enum -> string generator
  const std::string DeepDescriptorFileName(EDEEP_DESCRIPTOR descriptor) {
    const std::string& modelDir = "/home/nomoko/Downloads/debug/openMVGOxford/src/openMVG/features/deep/networks/";
    switch (descriptor) {
      case SIAM_2_STREAM_DESC_NOTRE_DAME: {
        const std::string modelStr(modelDir + "siam2stream/siam2stream_desc_notredame.bin");
        return modelStr;
      }
      case SIAM_2_STREAM_DESC_YOSEMITE: {
        const std::string modelStr(modelDir + "siam2stream/siam2stream_desc_yosemite.bin");
        return modelStr;
      }
      case SIAM_2_STREAM_DESC_LIBERTY: {
        const std::string modelStr(modelDir + "siam2stream/siam2stream_desc_liberty.bin");
        return modelStr;
      }
      case SIAM_DESC_NOTRE_DAME: {
        const std::string modelStr(modelDir + "siam/siam_desc_notredame.bin");
        return modelStr;
      }
      case SIAM_DESC_YOSEMITE: {
        const std::string modelStr(modelDir + "siam/siam_desc_yosemite.bin");
        return modelStr;
      }
      case SIAM_DESC_LIBERTY: {
        const std::string modelStr(modelDir + "siam/siam_desc_liberty.bin");
        return modelStr;
      }
      case PNNET: {
        const std::string modelStr(modelDir + "pnnet/");
        return modelStr;
      }
      default: {
        const std::string& emptyStr("");
        return emptyStr;
      }
    }
  }

  DEEPParams(
    EDEEP_DESCRIPTOR eDeepDescriptor = SIAM_2_STREAM_DESC_YOSEMITE
  ):eDeepDescriptor_(eDeepDescriptor),
  eDeepDescriptorFileName_(DeepDescriptorFileName(eDeepDescriptor)){}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(eDeepDescriptor_);
  }

  // Parameters
  EDEEP_DESCRIPTOR eDeepDescriptor_;
  const std::string eDeepDescriptorFileName_;
};

class DEEP_Image_describer : public Image_describer
{
public:
  DEEP_Image_describer(
  const DEEPParams & params = DEEPParams()
  ):Image_describer(), params_(params) {
  deepClassifier_ = new DeepClassifierTHNets(params_.eDeepDescriptorFileName_);
  }
  ~DEEP_Image_describer() {
    delete deepClassifier_;
  }
  
  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
  {
    return true;
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = nullptr ) override
  {
    
    std::vector<cv::KeyPoint> kpts;
    float* descriptors = deepClassifier_->describeOpenMVG(image, kpts);

    Allocate(regions);
    switch (params_.eDeepDescriptor_)
    {
      case SIAM_2_STREAM_DESC_NOTRE_DAME:
      case SIAM_2_STREAM_DESC_YOSEMITE:
      case SIAM_2_STREAM_DESC_LIBERTY:
      {
        // Build alias to cached data
        DEEP_Float_512_Regions * regionsCasted = dynamic_cast<DEEP_Float_512_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          const SIOPointFeature& ptDeep = SIOPointFeature(kpts[i].pt.x, kpts[i].pt.y, kpts[i].size, kpts[i].angle);

          // Feature masking
          if (mask)
          {
          const image::Image<unsigned char> & maskIma = *mask;
          if (maskIma(ptDeep.y(), ptDeep.x()) == 0)
            continue;
          }
          regionsCasted->Features()[i] = ptDeep;
          // Store descriptors
          for (int j = 0; j < 512; j++) {
            const unsigned int index = i * 512 + j;
            const float descriptor = static_cast<float>(descriptors[index]);
            regionsCasted->Descriptors()[i][j] = descriptor;
          }
        }
        break;
      }
      case SIAM_DESC_NOTRE_DAME:
      case SIAM_DESC_YOSEMITE:
      case SIAM_DESC_LIBERTY:
      {
        // Build alias to cached data
        DEEP_Float_256_Regions * regionsCasted = dynamic_cast<DEEP_Float_256_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          const SIOPointFeature& ptDeep = SIOPointFeature(kpts[i].pt.x, kpts[i].pt.y, kpts[i].size, kpts[i].angle);

          // Feature masking
          if (mask)
          {
          const image::Image<unsigned char> & maskIma = *mask;
          if (maskIma(ptDeep.y(), ptDeep.x()) == 0)
            continue;
          }
          regionsCasted->Features()[i] = ptDeep;
          // Store descriptors
          for (int j = 0; j < 256; j++) {
            const unsigned int index = i * 256 + j;
            float descriptor = static_cast<float>(descriptors[index]);
            regionsCasted->Descriptors()[i][j] = descriptor;
          }
        }
        break;
      }
      case PNNET:
      {
        // Build alias to cached data
        DEEP_Float_128_Regions * regionsCasted = dynamic_cast<DEEP_Float_128_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          const SIOPointFeature& ptDeep = SIOPointFeature(kpts[i].pt.x, kpts[i].pt.y, kpts[i].size, kpts[i].angle);

          // Feature masking
          if (mask)
          {
          const image::Image<unsigned char> & maskIma = *mask;
          if (maskIma(ptDeep.y(), ptDeep.x()) == 0)
            continue;
          }
          regionsCasted->Features()[i] = ptDeep;
          // Store descriptors
          for (int j = 0; j < 128; j++) {
            const unsigned int index = i * 128 + j;
            float descriptor = static_cast<float>(descriptors[index]);
            regionsCasted->Descriptors()[i][j] = descriptor;
          }
        }
        break;
      }
      default:
      {
        break;
      }
    }
    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const override
  {
  switch(params_.eDeepDescriptor_)
  {
      case SIAM_2_STREAM_DESC_NOTRE_DAME:
      case SIAM_2_STREAM_DESC_YOSEMITE:
      case SIAM_2_STREAM_DESC_LIBERTY:
        return regions.reset(new DEEP_Float_512_Regions);
      case SIAM_DESC_NOTRE_DAME:
      case SIAM_DESC_YOSEMITE:
      case SIAM_DESC_LIBERTY:
        return regions.reset(new DEEP_Float_256_Regions);
      case PNNET:
        return regions.reset(new DEEP_Float_128_Regions);
      break;
  }
  }

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(
     cereal::make_nvp("params", params_),
     cereal::make_nvp("bOrientation", bOrientation_));
  }


private:
  DEEPParams params_;
  bool bOrientation_;
  DeepClassifier* deepClassifier_;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::DEEP_Image_describer, "DEEP_Image_describer");

#endif // OPENMVG_FEATURES_DEEP_IMAGE_DESCRIBER_HPP
