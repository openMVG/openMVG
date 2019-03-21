#ifndef RICH_SIFT_DESCRIBER_HPP
#define RICH_SIFT_DESCRIBER_HPP

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/image/image_container.hpp"
#include "nonFree/sift/SIFT_describer.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>

extern "C" {
#include "nonFree/sift/vl/sift.h"
}

namespace openMVG {
namespace features {

class Rich_SIFT_Image_describer : public Image_describer
{
public:

  using Regions_type = Rich_SIFT_Regions;

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
      _first_octave(first_octave),
      _num_octaves(num_octaves),
      _num_scales(num_scales),
      _edge_threshold(edge_threshold),
      _peak_threshold(peak_threshold),
      _root_sift(root_sift) {}

    template<class Archive>
    inline void serialize( Archive & ar );

    // Parameters
    int _first_octave;      // Use original image, or perform an upscale if == -1
    int _num_octaves;       // Max octaves count
    int _num_scales;        // Scales per octave
    float _edge_threshold;  // Max ratio of Hessian eigenvalues
    float _peak_threshold;  // Min contrast
    bool _root_sift;        // see [1]
  };

  //--
  // Constructor
  //--
  Rich_SIFT_Image_describer
  (
    const Params & params = Params(),
    bool bOrientation = true
  ):Image_describer(), _params(params), _bOrientation(bOrientation)
  {
    vl_constructor();
  }

  ~Rich_SIFT_Image_describer()
  {
    vl_destructor();
  }

  bool Set_configuration_preset(EDESCRIBER_PRESET preset) override
  {
    switch (preset)
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
  @param regions The detected regions and attributes
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe
  (
    const image::Image<unsigned char> & image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = nullptr
  ) override
  {
    Regions_type *regions_type = static_cast<Regions_type*>(regions.get());
    if (regions_type == nullptr || regions_type->RegionCount() == 0)
    {
      regions = DescribeSIFT(image, mask);
    }
    else
    {
      DescribeSIFTGivenRegion(image, regions, mask);
    }
    return regions != nullptr;
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  @return regions The detected regions and attributes (the caller must delete the allocated data)
  */
  std::unique_ptr<Regions> Describe(
    const image::Image<unsigned char>& image,
    const image::Image<unsigned char>* mask = nullptr
  ) override
  {
    return DescribeSIFT(image, mask);
  }

  /**
  @brief compute given regions' attributes (description)
  @param image Image.
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  (the caller must delete the allocated data)
  */
  void DescribeSIFTGivenRegion(
      const image::Image<unsigned char> &image,
      std::unique_ptr<Regions> &regions,
      const image::Image<unsigned char>* mask = nullptr
  )
  {
    Regions_type *regions_type = static_cast<Regions_type*>(regions.get());

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

    // reserve some memory for faster keypoint saving
    auto& features = regions_type->Features();
    regions_type->Descriptors().clear();
    regions_type->Descriptors().reserve(2000);

    int curr_feature_index = 0;
    int octave_index = 0;
    while (true) {
      // Update gradient before launching parallel extraction
      vl_sift_update_gradient(filt);

      while (curr_feature_index < features.size() &&
        features[curr_feature_index].octave() == octave_index) {

        const auto& f = features[curr_feature_index];
        VlSiftKeypoint key;
        key.o = f.octave();
        key.ix = f.ix();
        key.iy = f.iy();
        key.is = f.is();
        key.x = f.x();
        key.y = f.y();
        key.s = f.s();
        key.sigma = f.scale();

        // Feature masking
        if (mask)
        {
          const image::Image<unsigned char> & maskIma = *mask;
          if (maskIma(key.y, key.x) == 0)
          {
            ++curr_feature_index;
            continue;
          }
        }

        vl_sift_calc_keypoint_descriptor(filt, &descr[0], &key, f.orientation());

        siftDescToUChar(&descr[0], descriptor, _params._root_sift);
        regions_type->Descriptors().push_back(descriptor);

        ++curr_feature_index;
      }
      if (vl_sift_process_next_octave(filt))
        break; // Last octave

      ++octave_index;
    }
    vl_sift_delete(filt);
  }

  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  @return regions The detected regions and attributes (the caller must delete the allocated data)
  */
  std::unique_ptr<Regions_type> DescribeSIFT(
      const image::Image<unsigned char>& image,
      const image::Image<unsigned char>* mask = nullptr
  )
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

    // Build alias to cached data
    auto regions = std::unique_ptr<Regions_type>(new Regions_type);

    // reserve some memory for faster keypoint saving
    regions->Features().reserve(2000);
    regions->Descriptors().reserve(2000);

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
          const Rich_SIOPointFeature fp(keys[i].x, keys[i].y,
            keys[i].sigma, static_cast<float>(angles[q]), keys[i].o,
            keys[i].ix, keys[i].iy, keys[i].is, keys[i].s);

          siftDescToUChar(&descr[0], descriptor, _params._root_sift);
          #ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
          #endif
          {
            regions->Descriptors().push_back(descriptor);
            regions->Features().push_back(fp);
          }
        }
      }
      if (vl_sift_process_next_octave(filt))
        break; // Last octave
    }
    vl_sift_delete(filt);

    return regions;
  }

  std::unique_ptr<Regions> Allocate() const override
  {
    return std::unique_ptr<Regions_type>(new Regions_type);
  }

  template<class Archive>
  inline void serialize( Archive & ar );

private:
  Params _params;
  bool _bOrientation;
};

} // namespace features
} // namespace openMVG

#endif // RICH_SIFT_DESCRIBER_HPP
