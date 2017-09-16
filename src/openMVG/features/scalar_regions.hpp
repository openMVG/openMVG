// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_SCALAR_REGIONS_HPP
#define OPENMVG_FEATURES_SCALAR_REGIONS_HPP

#include <typeinfo>

#include "openMVG/features/regions.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/matching/metric.hpp"

namespace openMVG {
namespace features {

/// Specialization of the abstract Regions class to handle Scalar descriptors
template<typename FeatT, typename T, size_t L>
class Scalar_Regions : public Regions
{
public:

  //-- Type alias
  //--

  /// Region type
  using FeatureT = FeatT;
  /// Region descriptor
  using DescriptorT = Descriptor<T, L>;

  /// Container for multiple regions
  using FeatsT = std::vector<FeatureT>;
  /// Container for multiple regions description
  using DescsT = std::vector<DescriptorT, Eigen::aligned_allocator<DescriptorT>>;

  //-- Class functions
  //--

  bool IsScalar() const override {return true;}
  bool IsBinary() const override {return false;}
  std::string Type_id() const override {return typeid(T).name();}
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
    return PointFeatures(vec_feats_.begin(), vec_feats_.end());
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
    return new Scalar_Regions();
  }

  // Return the L2 distance between two descriptors
  double SquaredDescriptorDistance(size_t i, const Regions * regions, size_t j) const override
  {
    assert(i < vec_descs_.size());
    assert(regions);
    assert(j < regions->RegionCount());

    const Scalar_Regions<FeatT, T, L> * regionsT = dynamic_cast<const Scalar_Regions<FeatT, T, L> *>(regions);
    matching::L2<T> metric;
    return metric(vec_descs_[i].data(), regionsT->vec_descs_[j].data(), DescriptorT::static_size);
  }

  /// Add the Inth region to another Region container
  void CopyRegion(size_t i, Regions * region_container) const override
  {
    assert(i < vec_feats_.size() && i < vec_descs_.size());
    static_cast<Scalar_Regions<FeatT, T, L> *>(region_container)->vec_feats_.push_back(vec_feats_[i]);
    static_cast<Scalar_Regions<FeatT, T, L> *>(region_container)->vec_descs_.push_back(vec_descs_[i]);
  }

private:
  //--
  //-- internal data
  FeatsT vec_feats_; // region features
  DescsT vec_descs_; // region descriptions
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_SCALAR_REGIONS_HPP
