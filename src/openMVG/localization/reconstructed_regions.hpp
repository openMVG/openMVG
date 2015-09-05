
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

//  template<class Archive>
//  void serialize(Archive & ar)
//  {
//    ar(_vec_feats);
//    ar(_vec_descs);
//  }

  void filterRegions(const std::vector<FeatureInImage>& featuresInImage)
  {
    features::Scalar_Regions<FeatT, T, L> newRegions;
    newRegions.Features().resize(featuresInImage.size());
    newRegions.Descriptors().resize(featuresInImage.size());
    _associated3dPoint.resize(featuresInImage.size());
    for(std::size_t i = 0; i < featuresInImage.size(); ++i)
    {
      const FeatureInImage & feat = featuresInImage[i];
      newRegions.Features()[i] = _regions.Features()[feat._featureIndex];
      newRegions.Descriptors()[i] = _regions.Descriptors()[feat._featureIndex];
      _mapFullToLocal[feat._featureIndex] = i;
      _associated3dPoint[i] = feat._point3dId;
    }
    _regions.swap(newRegions);
  }

public:
  //--
  //-- internal data
  features::Scalar_Regions<FeatT, T, L> _regions;
  std::vector<IndexT> _associated3dPoint;
  std::map<IndexT, IndexT> _mapFullToLocal; 
};

} // namespace features
} // namespace openMVG
