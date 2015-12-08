
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/regions.hpp"
#include "openMVG/matching/metric.hpp"

#ifdef HAVE_CCTAG
#include "openMVG/features/cctag/CCTAG_describer.hpp"
#endif

#include "cereal/types/vector.hpp"

#include <string>
#include <typeinfo>

namespace openMVG {
namespace localization {

struct FeatureInImage
{
  FeatureInImage(IndexT featureIndex, IndexT point3dId)
    : _featureIndex(featureIndex)
    , _point3dId(point3dId)
  {}
  
  IndexT _featureIndex;
  IndexT _point3dId;
  
  bool operator<(const FeatureInImage& other) const
  {
    return _featureIndex < other._featureIndex;
  }
};


/// Specialization of the abstract Regions class to handle Scalar descriptors
template<typename FeatT, typename T, size_t L>
class Reconstructed_Regions
{
public:
  // Region type
  typedef FeatT FeatureT;
  // Region descriptor
  typedef features::Descriptor<T, L> DescriptorT;
  
//  template<class Archive>
//  void serialize(Archive & ar)
//  {
//    ar(_vec_feats);
//    ar(_vec_descs);
//  }
  
  void filterRegions(const std::vector<FeatureInImage>& featuresInImage)
  {
    features::Scalar_Regions<FeatT, T, L> newRegions;
    newRegions.Features().reserve(featuresInImage.size());
    newRegions.Descriptors().reserve(featuresInImage.size());
    _associated3dPoint.reserve(featuresInImage.size());
    for(std::size_t i = 0; i < featuresInImage.size(); ++i)
    {
      const FeatureInImage & feat = featuresInImage[i];
      newRegions.Features().push_back(_regions.Features()[feat._featureIndex]);
      newRegions.Descriptors().push_back(_regions.Descriptors()[feat._featureIndex]);
      _mapFullToLocal[feat._featureIndex] = i; //todo@Simone check if map is already initialized
      _associated3dPoint.push_back(feat._point3dId);
    }
    _regions.swap(newRegions);
  }
  
  
#ifdef HAVE_CCTAG
  
  void filterCCTagRegions(const std::vector<FeatureInImage>& featuresInImage)
  {
    features::Scalar_Regions<FeatT, T, L> newRegions;
    newRegions.Features().reserve(featuresInImage.size());
    newRegions.Descriptors().reserve(featuresInImage.size());
    _associated3dPoint.reserve(featuresInImage.size());
    for(std::size_t i = 0; i < featuresInImage.size(); ++i)
    {
      const FeatureInImage & feat = featuresInImage[i];  
      // if the descriptor is a CCTag Descriptor, then add the region
      IndexT cctagId = features::getCCTagId(_regions.Descriptors()[feat._featureIndex]);
      if( cctagId != UndefinedIndexT)
      {
        newRegions.Features().push_back(_regions.Features()[feat._featureIndex]);
        newRegions.Descriptors().push_back(_regions.Descriptors()[feat._featureIndex]);
        _mapFullToLocal[feat._featureIndex] = i; //@TODO check if map is already initialized
        _associated3dPoint.push_back(feat._point3dId); 
      }
    }
    _regions.swap(newRegions);
  }
  
  void updateLandmarksVisibility(std::vector<bool> & presentIds) const
  {
    assert(presentIds.size()==128);
    for (const auto & desc : _regions.Descriptors())
    {
      IndexT cctagId = features::getCCTagId(desc);
      if( cctagId != UndefinedIndexT) // todo put an assert instead ?
        presentIds[cctagId] = true;
    }
  }
#endif    

public:
  //--
  //-- internal data
  features::Scalar_Regions<FeatT, T, L> _regions;
  std::vector<IndexT> _associated3dPoint;
  std::map<IndexT, IndexT> _mapFullToLocal; 
};

} // namespace features
} // namespace openMVG
