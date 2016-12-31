// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/*

== Patent Warning and License =================================================

The SIFT method is patented

    [2] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89

 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.

The implementation is based on

    [1] "Anatomy of the SIFT Method."
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_anatomy_sift/

Changes are:
- The hierarchical scale space code can be run on it's own,
- Keypoint detection and description is split in two separate modules,
- the code can run per Octave (less memory consuming),
- some computation can be run in parallel.

*/

#ifndef OPENMVG_FEATURES_SIFT_SIFT_ANATOMY_IMAGE_DESCRIBER_HPP
#define OPENMVG_FEATURES_SIFT_SIFT_ANATOMY_IMAGE_DESCRIBER_HPP

#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/hierarchical_gaussian_scale_space.hpp"
#include "openMVG/features/sift/sift_DescriptorExtractor.hpp"
#include "openMVG/features/sift/sift_keypoint.hpp"
#include "openMVG/features/sift/sift_KeypointExtractor.hpp"

#include <cereal/cereal.hpp>

#include <iostream>
#include <numeric>


namespace openMVG {
namespace features {

class SIFT_Anatomy_Image_describer : public Image_describer
{
public:
  struct Params
  {
    Params(
      int first_octave = 0,
      int num_octaves = 6,
      int num_scales = 3,
      float edge_threshold = 10.0f,
      float peak_threshold = 0.04f,
      bool root_sift = true
    ):
      first_octave_(first_octave),
      num_octaves_(num_octaves),
      num_scales_(num_scales),
      edge_threshold_(edge_threshold),
      peak_threshold_(peak_threshold),
      root_sift_(root_sift) {}

    template<class Archive>
    void serialize( Archive & ar )
    {
      ar(
        cereal::make_nvp("first_octave", first_octave_),
        cereal::make_nvp("num_octaves",num_octaves_),
        cereal::make_nvp("num_scales",num_scales_),
        cereal::make_nvp("edge_threshold",edge_threshold_),
        cereal::make_nvp("peak_threshold",peak_threshold_),
        cereal::make_nvp("root_sift",root_sift_));
    }

    // Parameters
    int first_octave_;      // Use original image, or perform an upscale if == -1
    int num_octaves_;       // Max octaves count
    int num_scales_;        // Scales per octave
    float edge_threshold_;  // Max ratio of Hessian eigenvalues
    float peak_threshold_;  // Min contrast
    bool root_sift_;        // see [1]
  };

  SIFT_Anatomy_Image_describer
  (
    const Params params = Params()
  )
  :Image_describer(), params_(params)
  {}

  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
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
    const int w = image.Width(), h = image.Height();
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
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const override
  {
    regions.reset( new SIFT_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(cereal::make_nvp("params", params_));
  }

private:
  Params params_;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_Anatomy_Image_describer, "SIFT_Anatomy_Image_describer");

#endif // OPENMVG_FEATURES_SIFT_SIFT_ANATOMY_IMAGE_DESCRIBER_HPP