// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/features/sift/hierarchical_gaussian_scale_space.hpp"
#include "openMVG/features/sift/sift_DescriptorExtractor.hpp"
#include "openMVG/features/sift/sift_KeypointExtractor.hpp"

#include <cereal/cereal.hpp>

namespace openMVG {
namespace features {

SIFT_Anatomy_Image_describer::Params::Params(
  int first_octave,
  int num_octaves,
  int num_scales,
  float edge_threshold,
  float peak_threshold,
  bool root_sift
):
  first_octave_(first_octave),
  num_octaves_(num_octaves),
  num_scales_(num_scales),
  edge_threshold_(edge_threshold),
  peak_threshold_(peak_threshold),
  root_sift_(root_sift) 
{
}

template<class Archive>
void SIFT_Anatomy_Image_describer::Params::serialize( Archive & ar )
{
  ar(
    cereal::make_nvp("first_octave", first_octave_),
    cereal::make_nvp("num_octaves",num_octaves_),
    cereal::make_nvp("num_scales",num_scales_),
    cereal::make_nvp("edge_threshold",edge_threshold_),
    cereal::make_nvp("peak_threshold",peak_threshold_),
    cereal::make_nvp("root_sift",root_sift_));
}

SIFT_Anatomy_Image_describer::SIFT_Anatomy_Image_describer
(
  const Params params
)
:Image_describer(), params_(params)
{
}

bool SIFT_Anatomy_Image_describer::Set_configuration_preset(EDESCRIBER_PRESET preset)
{
  switch(preset)
  {
  case NORMAL_PRESET:
    params_.peak_threshold_ = 0.04f;
  break;
  case HIGH_PRESET:
    params_.peak_threshold_ = 0.01f;
  break;
  case ULTRA_PRESET:
    params_.peak_threshold_ = 0.01f;
    params_.first_octave_ = -1;
  break;
  default:
    return false;
  }
  return true;
}

bool SIFT_Anatomy_Image_describer::Describe
(
  const image::Image<unsigned char>& image,
  std::unique_ptr<Regions> &regions,
  const image::Image<unsigned char> * mask
)
{
  // Convert to float in range [0;1]
  const image::Image<float> If(image.GetMat().cast<float>()/255.0f);

  // compute sift keypoints
  Allocate(regions);

  // Build alias to cached data
  SIFT_Regions * regionsCasted = dynamic_cast<SIFT_Regions*>(regions.get());
  {
    using namespace openMVG::features::sift;
    const int supplementary_images = 3;
    // => in order to ensure each gaussian slice is used in the process 3 extra images are required:
    // +1 for dog computation
    // +2 for 3d discrete extrema definition

    HierarchicalGaussianScaleSpace octave_gen(
      params_.num_octaves_,
      params_.num_scales_,
      (params_.first_octave_ == -1)
      ? GaussianScaleSpaceParams(1.6f/2.0f, 1.0f/2.0f, 0.5f, supplementary_images)
      : GaussianScaleSpaceParams(1.6f, 1.0f, 0.5f, supplementary_images));
    octave_gen.SetImage( If );

    std::vector<Keypoint> keypoints;
    keypoints.reserve(5000);
    Octave octave;
    while ( octave_gen.NextOctave( octave ) )
    {
      std::vector< Keypoint > keys;
      // Find Keypoints
      SIFT_KeypointExtractor keypointDetector(
        params_.peak_threshold_ / octave_gen.NbSlice(),
        params_.edge_threshold_);
      keypointDetector(octave, keys);
      // Find Keypoints orientation and compute their description
      Sift_DescriptorExtractor descriptorExtractor;
      descriptorExtractor(octave, keys);

      // Concatenate the found keypoints
      std::move(keys.begin(), keys.end(), std::back_inserter(keypoints));
    }
    for (const auto & k : keypoints)
    {
      // Feature masking
      if (mask)
      {
        const image::Image<unsigned char> & maskIma = *mask;
        if (maskIma(k.y, k.x) == 0)
          continue;
      }

      Descriptor<unsigned char, 128> descriptor;
      descriptor << (k.descr.cast<unsigned char>());
      {
        regionsCasted->Descriptors().emplace_back(descriptor);
        regionsCasted->Features().emplace_back(k.x, k.y, k.sigma, k.theta);
      }
    }
  }
  return true;
}

void SIFT_Anatomy_Image_describer::Allocate(std::unique_ptr<Regions> &regions) const
{
  regions.reset( new SIFT_Regions );
}

template<class Archive>
void SIFT_Anatomy_Image_describer::serialize( Archive & ar )
{
  ar(cereal::make_nvp("params", params_));
}

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_Anatomy_Image_describer, "SIFT_Anatomy_Image_describer");
CEREAL_REGISTER_POLYMORPHIC_RELATION(openMVG::features::Image_describer, openMVG::features::SIFT_Anatomy_Image_describer)
