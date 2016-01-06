#pragma once

#include <openMVG/features/descriptor.hpp>
#include <openMVG/features/image_describer.hpp>
#include <openMVG/features/regions_factory.hpp>

extern "C" {
#include "nonFree/sift/vl/sift.h"
}

#include <cereal/cereal.hpp>

#include <iostream>
#include <numeric>

namespace openMVG {
namespace features {

struct SiftParams
{
  SiftParams(
    int first_octave = 0,
    int num_octaves = 6,
    int num_scales = 3,
    float edge_threshold = 10.0f,
    float peak_threshold = 0.04f,
    //
    std::size_t gridSize = 4,
    std::size_t maxTotalKeypoints = 1000,
    //
    bool root_sift = true
  ):
    _first_octave(first_octave),
    _num_octaves(num_octaves),
    _num_scales(num_scales),
    _edge_threshold(edge_threshold),
    _peak_threshold(peak_threshold),
    //
    _gridSize(gridSize),
    _maxTotalKeypoints(maxTotalKeypoints),
    //
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
      //
      cereal::make_nvp("grid_size", _gridSize),
      cereal::make_nvp("max_total_keypoints", _maxTotalKeypoints),
      //
      cereal::make_nvp("root_sift",_root_sift));
  }

  // Parameters
  int _first_octave;      // Use original image, or perform an upscale if == -1
  int _num_octaves;       // Max octaves count
  int _num_scales;        // Scales per octave
  float _edge_threshold;  // Max ratio of Hessian eigenvalues
  float _peak_threshold;  // Min contrast
  //
  std::size_t _gridSize;
  std::size_t _maxTotalKeypoints;
  //
  bool _root_sift;        // see [1]
  
  bool setPreset(EDESCRIBER_PRESET preset)
  {
    switch(preset)
    {
    case LOW_PRESET:
    {
      _maxTotalKeypoints = 1000;
      _peak_threshold = 0.04f;
      _first_octave = 2;
      break;
    }
    case MEDIUM_PRESET:
    {
      _maxTotalKeypoints = 5000;
      _peak_threshold = 0.04f;
      _first_octave = 1;
      break;
    }
    case NORMAL_PRESET:
    {
      _maxTotalKeypoints = 10000;
      _peak_threshold = 0.04f;
      break;
    }
    case HIGH_PRESET:
    {
      _maxTotalKeypoints = 20000;
      _peak_threshold = 0.01f;
      break;
    }
    case ULTRA_PRESET:
    {
      _maxTotalKeypoints = 40000;
      _peak_threshold = 0.01f;
      _first_octave = -1;
      break;
    }
    default:
      return false;
    }
    return true;
  }
  
  
};

//convertSIFT
//////////////////////////////
template < typename TOut > 
inline void convertSIFT(
  const vl_sift_pix* descr,
  Descriptor<TOut,128> & descriptor,
  bool brootSift = false
  );

template <> 
inline void convertSIFT<float>(
  const vl_sift_pix* descr,
  Descriptor<float,128> &descriptor,
  bool brootSift)
{
  if(brootSift)
  {
    const float sum = std::accumulate(descr, descr + 128, 0.0f);
    for(int k = 0; k < 128; ++k)
      descriptor[k] = std::floor(512.f*sqrt(descr[k] / sum));
  }
  else
  {
    for(int k = 0; k < 128; ++k)
      descriptor[k] = std::floor(512.f*descr[k]);
  }
}

template <> 
inline void convertSIFT<unsigned char>(
  const vl_sift_pix* descr,
  Descriptor<unsigned char,128> & descriptor,
  bool brootSift)
{
  if (brootSift)
  {
    // rootsift = sqrt( sift / sum(sift) );
    const float sum = std::accumulate(descr, descr+128, 0.0f);
    for (int k=0;k<128;++k)
      descriptor[k] = static_cast<unsigned char>(512.f*sqrt(descr[k]/sum));
  }
  else
  {
    for (int k=0;k<128;++k)
      descriptor[k] = static_cast<unsigned char>(512.f*descr[k]);
  }
}

/**
 * @brief Extract SIFT regions (in float or unsigned char).
 * 
 * @param image
 * @param regions
 * @param params
 * @param bOrientation
 * @param mask
 * @return 
 */
template < typename T >
bool extractSIFT(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const SiftParams& params,
    bool bOrientation,
    const image::Image<unsigned char> * mask)
{
  const int w = image.Width(), h = image.Height();
  //Convert to float
  const image::Image<float> imageFloat(image.GetMat().cast<float>());

  VlSiftFilt *filt = vl_sift_new(w, h, params._num_octaves, params._num_scales, params._first_octave);
  if (params._edge_threshold >= 0)
    vl_sift_set_edge_thresh(filt, params._edge_threshold);
  if (params._peak_threshold >= 0)
    vl_sift_set_peak_thresh(filt, 255*params._peak_threshold/params._num_scales);

  Descriptor<vl_sift_pix, 128> vlFeatDescriptor;
  Descriptor<T, 128> descriptor;

  // Process SIFT computation
  vl_sift_process_first_octave(filt, imageFloat.data());

  typedef Scalar_Regions<SIOPointFeature,T,128> SIFT_Region_T;
  regions.reset( new SIFT_Region_T );
  
  // Build alias to cached data
  SIFT_Region_T * regionsCasted = dynamic_cast<SIFT_Region_T*>(regions.get());
  // reserve some memory for faster keypoint saving
  regionsCasted->Features().reserve(2000);
  regionsCasted->Descriptors().reserve(2000);

  while (true) {
    vl_sift_detect(filt);

    VlSiftKeypoint const *keys  = vl_sift_get_keypoints(filt);
    const int nkeys = vl_sift_get_nkeypoints(filt);

    // Update gradient before launching parallel extraction
    vl_sift_update_gradient(filt);

    std::cout << "TEST octave_width: " << filt->octave_width << std::endl;
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for private(vlFeatDescriptor, descriptor)
    #endif
    for (int i = 0; i < nkeys; ++i) {

      // Feature masking
      if (mask)
      {
        const image::Image<unsigned char> & maskIma = *mask;
        if (maskIma(keys[i].y, keys[i].x) > 0)
          continue;
      }

      double angles [4] = {0.0, 0.0, 0.0, 0.0};
      int nangles = 1; // by default (1 upright feature)
      if (bOrientation)
      { // compute from 1 to 4 orientations
        nangles = vl_sift_calc_keypoint_orientations(filt, angles, keys+i);
      }

      for (int q=0 ; q < nangles ; ++q) {
        vl_sift_calc_keypoint_descriptor(filt, &vlFeatDescriptor[0], keys+i, angles[q]);
        const SIOPointFeature fp(keys[i].x, keys[i].y,
          keys[i].sigma, static_cast<float>(angles[q]));

        convertSIFT<T>(&vlFeatDescriptor[0], descriptor, params._root_sift);
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

  const auto& features = regionsCasted->Features();
  const auto& descriptors = regionsCasted->Descriptors();
  
  //Sorting the extracted features according to their scale
  {
    std::vector<std::size_t> indexSort(features.size());
    std::iota(indexSort.begin(), indexSort.end(), 0);
    std::sort(indexSort.begin(), indexSort.end(), [&](std::size_t a, std::size_t b){ return features[a].scale() > features[b].scale(); });
    
    std::vector<typename SIFT_Region_T::FeatureT> sortedFeatures(features.size());
    std::vector<typename SIFT_Region_T::DescriptorT> sortedDescriptors(features.size());
    for(std::size_t i: indexSort)
    {
      sortedFeatures[i] = features[indexSort[i]];
      sortedDescriptors[i] = descriptors[indexSort[i]];
    }
    regionsCasted->Features().swap(sortedFeatures);
    regionsCasted->Descriptors().swap(sortedDescriptors);
  }
    
  // Grid filtering of the keypoints to ensure a global repartition
  if(params._gridSize && params._maxTotalKeypoints)
  {
    // Only filter features if we have more features than the maxTotalKeypoints
    if(features.size() > params._maxTotalKeypoints)
    {
      std::vector<typename SIFT_Region_T::FeatureT> filtered_keypoints;
      std::vector<typename SIFT_Region_T::FeatureT> rejected_keypoints;
      filtered_keypoints.reserve(std::min(features.size(), params._maxTotalKeypoints));
      rejected_keypoints.reserve(features.size());

      const std::size_t sizeMat = params._gridSize*params._gridSize;
      std::vector<std::size_t> countFeatPerCell(sizeMat, 0);
      for (int Indice = 0; Indice < sizeMat; Indice++) {
    	  countFeatPerCell[Indice] = 0;
      }
      const std::size_t keypointsPerCell = params._maxTotalKeypoints / sizeMat;
      const double regionWidth = w / double(params._gridSize);
      const double regionHeight = h / double(params._gridSize);

      //std::cout << "Grid filtering -- keypointsPerCell: " << keypointsPerCell
      //          << ", regionWidth: " << regionWidth
      //          << ", regionHeight: " << regionHeight << std::endl;

      for(const auto& keypoint: features)
      {
        const std::size_t cellX = std::min(std::size_t(keypoint.x() / regionWidth), params._gridSize);
        const std::size_t cellY = std::min(std::size_t(keypoint.y() / regionHeight), params._gridSize);
        //std::cout << "- keypoint.pt.x: " << keypoint.pt.x << ", keypoint.pt.y: " << keypoint.pt.y << std::endl;
        //std::cout << "- cellX: " << cellX << ", cellY: " << cellY << std::endl;
        //std::cout << "- countFeatPerCell: " << countFeatPerCell << std::endl;
        //std::cout << "- gridSize: " << _params.gridSize << std::endl;

        std::size_t &count = countFeatPerCell[cellX*params._gridSize + cellY];
        ++count;
        if(count < keypointsPerCell)
          filtered_keypoints.push_back(keypoint);
        else
          rejected_keypoints.push_back(keypoint);
      }
      // If we don't have enough features (less than maxTotalKeypoints) after the grid filtering (empty regions in the grid for example).
      // We add the best other ones, without repartition constraint.
      if( filtered_keypoints.size() < params._maxTotalKeypoints )
      {
        const std::size_t remainingElements = std::min(rejected_keypoints.size(), params._maxTotalKeypoints - filtered_keypoints.size());
        std::cout << "Grid filtering -- Copy remaining points: " << remainingElements << std::endl;
        filtered_keypoints.insert(filtered_keypoints.end(), rejected_keypoints.begin(), rejected_keypoints.begin() + remainingElements);
      }

      regionsCasted->Features().swap(filtered_keypoints);
      
    }
  }
  
  return true;
}

} //namespace openMVG
} //namespace features

