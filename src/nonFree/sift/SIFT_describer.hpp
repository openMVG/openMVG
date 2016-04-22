// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_H
#define OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_H

#include <cereal/cereal.hpp>
#include <iostream>
#include <numeric>

extern "C" {
#include "nonFree/sift/vl/sift.h"
}

namespace openMVG {
namespace features {

// Bibliography:
// [1] R. Arandjelović, A. Zisserman.
// Three things everyone should know to improve object retrieval. CVPR2012.

inline void siftDescToUChar(
  vl_sift_pix descr[128],
  Descriptor<unsigned char,128> & descriptor,
  bool brootSift = false)
{
  if (brootSift)  {
    // rootsift = sqrt( sift / sum(sift) );
    const float sum = accumulate(descr, descr+128, 0.0f);
    for (int k=0;k<128;++k)
      descriptor[k] = static_cast<unsigned char>(512.f*sqrt(descr[k]/sum));
  }
  else
    for (int k=0;k<128;++k)
    descriptor[k] = static_cast<unsigned char>(512.f*descr[k]);
}

struct SiftParams
{
  SiftParams(
    int first_octave = 0,
    int num_octaves = 6,
    int num_scales = 3,
    float edge_threshold = 10.0f,
    float peak_threshold = 0.04f,
    bool root_sift = true
  ):
    _first_octave(first_octave),
    _num_octaves(num_octaves),
    _num_scales(num_scales),
    _edge_threshold(edge_threshold),
    _peak_threshold(peak_threshold),
    _root_sift(root_sift) {}

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(
      cereal::make_nvp("first_octave", _first_octave),
      cereal::make_nvp("num_octaves",_num_octaves),
      cereal::make_nvp("num_scales",_num_scales),
      cereal::make_nvp("edge_threshold",_edge_threshold),
      cereal::make_nvp("peak_threshold",_peak_threshold),
      cereal::make_nvp("root_sift",_root_sift));
  }

  // Parameters
  int _first_octave;      // Use original image, or perform an upscale if == -1
  int _num_octaves;       // Max octaves count
  int _num_scales;        // Scales per octave
  float _edge_threshold;  // Max ratio of Hessian eigenvalues
  float _peak_threshold;  // Min contrast
  bool _root_sift;        // see [1]
};

class SIFT_Image_describer : public Image_describer
{
public:
  SIFT_Image_describer(const SiftParams & params = SiftParams(), bool bOrientation = true)
    :Image_describer(), _params(params), _bOrientation(bOrientation)
  {
    vl_constructor();
  }

  ~SIFT_Image_describer()
  {
    vl_destructor();
  }

  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
    switch(preset)
    {
    case NORMAL_PRESET:
      _params._peak_threshold = 0.04f;
    break;
    case HIGH_PRESET:
      _params._peak_threshold = 0.01f;
    break;
    case ULTRA_PRESET:
      _params._peak_threshold = 0.01f;
      _params._first_octave = -1;
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
  bool Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = NULL)
  {
    const int w = image.Width(), h = image.Height();
    //Convert to float
    const image::Image<float> If(image.GetMat().cast<float>());

    VlSiftFilt *filt = vl_sift_new(w, h,
      _params._num_octaves, _params._num_scales, _params._first_octave);
    if (_params._edge_threshold >= 0)
      vl_sift_set_edge_thresh(filt, _params._edge_threshold);
    if (_params._peak_threshold >= 0)
      vl_sift_set_peak_thresh(filt, 255*_params._peak_threshold/_params._num_scales);

    Descriptor<vl_sift_pix, 128> descr;
    Descriptor<unsigned char, 128> descriptor;

    // Process SIFT computation
    vl_sift_process_first_octave(filt, If.data());

    Allocate(regions);

    // Build alias to cached data
    SIFT_Regions * regionsCasted = dynamic_cast<SIFT_Regions*>(regions.get());
    // reserve some memory for faster keypoint saving
    regionsCasted->Features().reserve(2000);
    regionsCasted->Descriptors().reserve(2000);

    while (true) {
      vl_sift_detect(filt);

      VlSiftKeypoint const *keys  = vl_sift_get_keypoints(filt);
      const int nkeys = vl_sift_get_nkeypoints(filt);

      // Update gradient before launching parallel extraction
      vl_sift_update_gradient(filt);

      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for private(descr, descriptor)
      #endif
      for (int i = 0; i < nkeys; ++i) {

        // Feature masking
        if (mask)
        {
          const image::Image<unsigned char> & maskIma = *mask;
          if (maskIma(keys[i].y, keys[i].x) == 0)
            continue;
        }

        double angles [4] = {0.0, 0.0, 0.0, 0.0};
        int nangles = 1; // by default (1 upright feature)
        if (_bOrientation)
        { // compute from 1 to 4 orientations
          nangles = vl_sift_calc_keypoint_orientations(filt, angles, keys+i);
        }

        for (int q=0 ; q < nangles ; ++q) {
          vl_sift_calc_keypoint_descriptor(filt, &descr[0], keys+i, angles[q]);
          const SIOPointFeature fp(keys[i].x, keys[i].y,
            keys[i].sigma, static_cast<float>(angles[q]));

          siftDescToUChar(&descr[0], descriptor, _params._root_sift);
          #ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
          #endif
          {
            regionsCasted->Descriptors().push_back(descriptor);
            regionsCasted->Features().push_back(fp);
          }
        }
      }
      if (vl_sift_process_next_octave(filt))
        break; // Last octave
    }
    vl_sift_delete(filt);

    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
    regions.reset( new SIFT_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(
     cereal::make_nvp("params", _params),
     cereal::make_nvp("bOrientation", _bOrientation));
  }

private:
  SiftParams _params;
  bool _bOrientation;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFT_Image_describer, "SIFT_Image_describer");

#endif // OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_H
