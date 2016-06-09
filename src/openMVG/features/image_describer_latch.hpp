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

enum ELATCH_DESCRIPTOR
{
	LATCH_UNSIGNED,
	LATCH_BINARY
};

struct LATCHParams
{
  LATCHParams(
    ELATCH_DESCRIPTOR eLatchDescriptor = LATCH_UNSIGNED
  ):eLatchDescriptor_(eLatchDescriptor){}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(eLatchDescriptor_);
  }

  // Parameters
  ELATCH_DESCRIPTOR eLatchDescriptor_;
};

class LATCH_Image_describer : public Image_describer
{
public:
  LATCH_Image_describer(
	const LATCHParams & params = LATCHParams()
  ):Image_describer(), params_(params) {
        }

  // Don't need to really define this for the LATCH class yet, until more descriptors come out.
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
	std::vector<LatchClassifierKeypoint> kpts;
	unsigned int* descriptors;
	{
		LatchClassifierOpenMVG latchLocal;
		kpts = latchLocal.identifyFeaturePointsOpenMVG(image.GetMat());
		descriptors = std::move(latchLocal.getDescriptorSet1());
	}
	Allocate(regions);
	switch (params_.eLatchDescriptor_)
	{
	  case LATCH_UNSIGNED:
	  {
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
		  // Store descriptors
		  #ifdef OPENMVG_USE_OPENMP
			#pragma omp parallel for
		  #endif
		  for (int j = 0; j < 64; j++) {
			const unsigned int index = i * 64 + j;
			unsigned int descriptor = static_cast<unsigned int>(descriptors[index]);
				  regionsCasted->Descriptors()[i][j] = descriptor;
		  }
		}
	  }
	  break;
	  case LATCH_BINARY:
	  {
		// Build alias to cached data
		AKAZE_Binary_Regions * regionsCasted = dynamic_cast<AKAZE_Binary_Regions*>(regions.get());
		regionsCasted->Features().resize(kpts.size());
		regionsCasted->Descriptors().resize(kpts.size());

	  #ifdef OPENMVG_USE_OPENMP
		#pragma omp parallel for
	  #endif
		for (size_t i = 0; i < static_cast<int>(kpts.size()); ++i)
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
		  // Store descriptors
		  for (size_t j = 0; j < 16; j++) {
			const unsigned int descriptorIndex = i * 64 + j;
				  unsigned int descriptor = static_cast<unsigned int>(descriptors[descriptorIndex]);
				  for (size_t k = 0; k < 4; k++) {
					unsigned char descriptorRawByte = static_cast<unsigned char>((descriptor >> (24-k*3)) & 0xFF);
					const unsigned int descriptorRawIndex = j * 4 + k;
			  regionsCasted->Descriptors()[i][descriptorRawIndex] = descriptorRawByte;
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
	switch(params_.eLatchDescriptor_)
	{
      case LATCH_UNSIGNED:
        return regions.reset(new LATCH_Unsigned_Int_Regions);
      break;
      case LATCH_BINARY:
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
  LATCHParams params_;
  bool bOrientation_;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::LATCH_Image_describer, "LATCH_Image_describer");

#endif // OPENMVG_FEATURES_LATCH_IMAGE_DESCRIBER_HPP
