// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP

#include <iostream>
#include <numeric>

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/features/akaze/AKAZE.hpp"
#include "openMVG/features/akaze/msurf_descriptor.hpp"
#include "openMVG/features/akaze/mldb_descriptor.hpp"
#include "openMVG/features/liop/liop_descriptor.hpp"
#include <cereal/cereal.hpp>

using namespace std;

namespace openMVG {
namespace features {

enum EAKAZE_DESCRIPTOR
{
  AKAZE_MSURF,
  AKAZE_LIOP,
  AKAZE_MLDB
};

class AKAZE_Image_describer : public Image_describer
{
public:

  struct Params
  {
    Params(
      const features::AKAZE::Params config = features::AKAZE::Params(),
      EAKAZE_DESCRIPTOR eAkazeDescriptor = AKAZE_MSURF
    ):options_(config), eAkazeDescriptor_(eAkazeDescriptor){}

    template<class Archive>
    void serialize(Archive & ar)
    {
      ar(options_, eAkazeDescriptor_);
    }

    // Parameters
    features::AKAZE::Params options_;
    EAKAZE_DESCRIPTOR eAkazeDescriptor_;
  };

  AKAZE_Image_describer(
    const Params & params = Params(),
    bool bOrientation = true
  ):Image_describer(), params_(params), bOrientation_(bOrientation) {}


  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
  {
    switch(preset)
    {
    case NORMAL_PRESET:
      params_.options_.fThreshold = features::AKAZE::Params().fThreshold;
    break;
    case HIGH_PRESET:
      params_.options_.fThreshold = features::AKAZE::Params().fThreshold/10.;
    break;
    case ULTRA_PRESET:
     params_.options_.fThreshold = features::AKAZE::Params().fThreshold/100.;
    break;
    default:
      return false;
    }
    return true;
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe
  (
    const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = nullptr
  ) override
  {
    params_.options_.fDesc_factor =
      (params_.eAkazeDescriptor_ == AKAZE_MSURF ||
      params_.eAkazeDescriptor_ == AKAZE_LIOP) ? 10.f*sqrtf(2.f)
      : 11.f*sqrtf(2.f); // MLDB

    AKAZE akaze(image, params_.options_);
    akaze.Compute_AKAZEScaleSpace();
    std::vector<AKAZEKeypoint> kpts;
    kpts.reserve(5000);
    akaze.Feature_Detection(kpts);
    akaze.Do_Subpixel_Refinement(kpts);

    Allocate(regions);

    switch(params_.eAkazeDescriptor_)
    {
      case AKAZE_MSURF:
      {
        // Build alias to cached data
        AKAZE_Float_Regions * regionsCasted = dynamic_cast<AKAZE_Float_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          AKAZEKeypoint ptAkaze = kpts[i];

          // Feature masking
          if (mask)
          {
            const image::Image<unsigned char> & maskIma = *mask;
            if (maskIma(ptAkaze.y, ptAkaze.x) == 0)
              continue;
          }

          const TEvolution & cur_slice = akaze.getSlices()[ptAkaze.class_id];

          if (bOrientation_)
            akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);
          else
            ptAkaze.angle = 0.0f;

          regionsCasted->Features()[i] =
            SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

          ComputeMSURFDescriptor(cur_slice.Lx, cur_slice.Ly, ptAkaze.octave,
            regionsCasted->Features()[i],
            regionsCasted->Descriptors()[i]);
        }
      }
      break;
      case AKAZE_LIOP:
      {
        // Build alias to cached data
        AKAZE_Liop_Regions * regionsCasted = dynamic_cast<AKAZE_Liop_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

        // Init LIOP extractor
        LIOP::Liop_Descriptor_Extractor liop_extractor;

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          AKAZEKeypoint ptAkaze = kpts[i];

          // Feature masking
          if (mask)
          {
            const image::Image<unsigned char> & maskIma = *mask;
            if (maskIma(ptAkaze.y, ptAkaze.x) > 0)
              continue;
          }

          const TEvolution & cur_slice = akaze.getSlices()[ptAkaze.class_id];

          if (bOrientation_)
            akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);
          else
            ptAkaze.angle = 0.0f;

          regionsCasted->Features()[i] =
            SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

          // Compute LIOP descriptor (do not need rotation computation, since
          //  LIOP descriptor is rotation invariant).
          // Rescale for LIOP patch extraction
          const SIOPointFeature fp =
              SIOPointFeature(ptAkaze.x, ptAkaze.y,
              ptAkaze.size/2.0, ptAkaze.angle);

          float desc[144];
          liop_extractor.extract(image, fp, desc);
          for(int j=0; j < 144; ++j)
            regionsCasted->Descriptors()[i][j] =
              static_cast<unsigned char>(desc[j] * 255.f +.5f);
        }
      }
      break;
      case AKAZE_MLDB:
      {
        // Build alias to cached data
        AKAZE_Binary_Regions * regionsCasted = dynamic_cast<AKAZE_Binary_Regions*>(regions.get());
        regionsCasted->Features().resize(kpts.size());
        regionsCasted->Descriptors().resize(kpts.size());

      #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for
      #endif
        for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
        {
          AKAZEKeypoint ptAkaze = kpts[i];

          // Feature masking
          if (mask)
          {
            const image::Image<unsigned char> & maskIma = *mask;
            if (maskIma(ptAkaze.y, ptAkaze.x) > 0)
              continue;
          }

          const TEvolution & cur_slice = akaze.getSlices()[ptAkaze.class_id];

          if (bOrientation_)
            akaze.Compute_Main_Orientation(ptAkaze, cur_slice.Lx, cur_slice.Ly);
          else
            ptAkaze.angle = 0.0f;

          regionsCasted->Features()[i] =
            SIOPointFeature(ptAkaze.x, ptAkaze.y, ptAkaze.size, ptAkaze.angle);

          // Compute MLDB descriptor
          Descriptor<bool,486> desc;
          ComputeMLDBDescriptor(cur_slice.cur, cur_slice.Lx, cur_slice.Ly,
            ptAkaze.octave, regionsCasted->Features()[i], desc);
          // convert the bool vector to the binary unsigned char array
          unsigned char * ptr = reinterpret_cast<unsigned char*>(&regionsCasted->Descriptors()[i]);
          memset(ptr, 0, regionsCasted->DescriptorLength()*sizeof(unsigned char));
          // For each byte
          for (int j = 0; j < std::ceil(486./8.); ++j, ++ptr) {
            // set the corresponding 8bits to the good values
            for (int iBit = 0; iBit < 8 && j*8+iBit < 486; ++iBit)  {
              *ptr |= desc[j*8+iBit] << iBit;
            }
          }
        }
      }
      break;
    }
    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const override
  {
    switch(params_.eAkazeDescriptor_)
    {
      case AKAZE_MSURF:
        return regions.reset(new AKAZE_Float_Regions);
      break;
      case AKAZE_LIOP:
        return regions.reset(new AKAZE_Liop_Regions);
      case AKAZE_MLDB:
       return regions.reset(new AKAZE_Binary_Regions);
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
  Params params_;
  bool bOrientation_;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::AKAZE_Image_describer, "AKAZE_Image_describer");

#endif // OPENMVG_FEATURES_AKAZE_IMAGE_DESCRIBER_HPP
