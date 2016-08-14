// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_ORB_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_ORB_IMAGE_DESCRIBER_HPP

#include <iostream>
#include <numeric>

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/features/orb/OrbClassifierOpenMVG.hpp"
#include <cereal/cereal.hpp>

using namespace std;

namespace openMVG {
namespace features {

enum EORB_DESCRIPTOR
{
	ORB
};

struct ORBParams
{
  ORBParams(
    EORB_DESCRIPTOR eOrbDescriptor = ORB
  ):eOrbDescriptor_(eOrbDescriptor){}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(eOrbDescriptor_);
  }

  // Parameters
  EORB_DESCRIPTOR eOrbDescriptor_;
};

class ORB_Image_describer : public Image_describer
{
public:
  ORB_Image_describer(
	const ORBParams & params = ORBParams()
  ):Image_describer(), params_(params), orb_() {
        }

  // Don't need to really define this for the ORB class yet, until more descriptors come out.
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
		std::vector<OrbClassifierKeypoint> kpts;
		Eigen::Matrix<unsigned char, Eigen::Dynamic, 32, Eigen::RowMajor> descriptors;
	
  	kpts = orb_.identifyFeaturePointsOpenMVG(image.GetMat());
		descriptors = orb_.returnDescriptors();
	
	Allocate(regions);
	switch (params_.eOrbDescriptor_)
	{
	  case ORB:
	  {
		// Build alias to cached data
		ORB_Binary_Regions * regionsCasted = dynamic_cast<ORB_Binary_Regions*>(regions.get());
		regionsCasted->Features().resize(kpts.size());
		regionsCasted->Descriptors().resize(kpts.size());

    #ifdef OPENMVG_USE_OPENMP
		#pragma omp parallel for
	  #endif
		for (int i = 0; i < static_cast<int>(kpts.size()); ++i)
		{
		  const OrbClassifierKeypoint ptLatch = kpts[i];

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
			for (int j = 0; j < 32; j++) {
				regionsCasted->Descriptors()[i][j] = descriptors(i,j);
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
	switch(params_.eOrbDescriptor_)
	{
      case ORB:
        return regions.reset(new ORB_Binary_Regions);
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
  ORBParams params_;
  bool bOrientation_;
  OrbClassifierOpenMVG orb_;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::ORB_Image_describer, "ORB_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::ORB_Image_describer)
#endif // OPENMVG_FEATURES_ORB_IMAGE_DESCRIBER_HPP
