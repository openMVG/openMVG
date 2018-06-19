// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_BINARY_REGIONS_HPP
#define OPENMVG_FEATURES_BINARY_REGIONS_HPP

#include <typeinfo>

#include "openMVG/features/regions.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/matching/metric.hpp"

namespace openMVG {
namespace features {

/// Binary_Regions represented as uchar based array
/// It's better to set L as a valid %8 value
///  for memory alignement and performance issue
template<typename FeatT, size_t L>
class Binary_Regions : public Regions
{
public:

  //-- Type alias
  //--

  /// Region
  using FeatureT = FeatT;
  /// Description of a region
  using DescriptorT = Descriptor<unsigned char, L>;

  /// Container for multiple regions
  using FeatsT = std::vector<FeatureT>;
  /// Container for multiple region descriptions
  using DescsT = std::vector<DescriptorT, Eigen::aligned_allocator<DescriptorT>>;

  //-- Class functions
  //--

  bool IsScalar() const override {return false;}
  bool IsBinary() const override {return true;}
  std::string Type_id() const override {return typeid(typename DescriptorT::bin_type).name();}
  size_t DescriptorLength() const override {return static_cast<size_t>(L);}

  /// Read from files the regions and their corresponding descriptors.
  bool Load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) override
  {
    return loadFeatsFromFile(sfileNameFeats, vec_feats_)
          & loadDescsFromBinFile(sfileNameDescs, vec_descs_);
  }

  /// Export in two separate files the regions and their corresponding descriptors.
  bool Save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const override
  {
    return saveFeatsToFile(sfileNameFeats, vec_feats_)
          & saveDescsToBinFile(sfileNameDescs, vec_descs_);
  }

  bool LoadFeatures(const std::string& sfileNameFeats) override
  {
    return loadFeatsFromFile(sfileNameFeats, vec_feats_);
  }

  PointFeatures GetRegionsPositions() const override
  {
    return {vec_feats_.cbegin(), vec_feats_.cend()};
  }

  Vec2 GetRegionPosition(size_t i) const override
  {
    return Vec2f(vec_feats_[i].coords()).cast<double>();
  }

  /// Return the number of defined regions
  size_t RegionCount() const override {return vec_feats_.size();}

  /// Mutable and non-mutable FeatureT getters.
  inline FeatsT & Features() { return vec_feats_; }
  inline const FeatsT & Features() const { return vec_feats_; }

  /// Mutable and non-mutable DescriptorT getters.
  inline DescsT & Descriptors() { return vec_descs_; }
  inline const DescsT & Descriptors() const { return vec_descs_; }

  const void * DescriptorRawData() const override { return &vec_descs_[0];}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(vec_feats_, vec_descs_);
  }

  Regions * EmptyClone() const override
  {
    return new Binary_Regions;
  }

  // Return the squared Hamming distance between two descriptors
  double SquaredDescriptorDistance(size_t i, const Regions * regions, size_t j) const override
  {
    assert(i < vec_descs_.size());
    assert(regions);
    assert(j < regions->RegionCount());

    const Binary_Regions<FeatT, L> * regionsT = dynamic_cast<const Binary_Regions<FeatT, L> *>(regions);
    matching::Hamming<unsigned char> metric;
    const typename matching::Hamming<unsigned char>::ResultType descDist =
      metric(vec_descs_[i].data(), regionsT->vec_descs_[j].data(), DescriptorT::static_size);
    return descDist * descDist;
  }

  /// Add the Inth region to another Region container
  void CopyRegion(size_t i, Regions * region_container) const override
  {
    assert(i < vec_feats_.size() && i < vec_descs_.size());
    static_cast<Binary_Regions<FeatT, L> *>(region_container)->vec_feats_.push_back(vec_feats_[i]);
    static_cast<Binary_Regions<FeatT, L> *>(region_container)->vec_descs_.push_back(vec_descs_[i]);
  }

private:
  //--
  //-- internal data
  FeatsT vec_feats_; // region features
  DescsT vec_descs_; // region descriptions
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_BINARY_REGIONS_HPP
