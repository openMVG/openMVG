
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_REGIONS_HPP
#define OPENMVG_FEATURES_REGIONS_HPP

#include "openMVG/numeric/numeric.h"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/matching/metric.hpp"
#include "cereal/types/vector.hpp"
#include <string>
#include <typeinfo>

namespace openMVG {
namespace features {

/// Describe an image a set of regions (position, ...) + attributes
/// Each region is described by a set of attributes (descriptor)
class Regions
{
public:

  virtual ~Regions() {}

  //--
  // IO - one file for region features, one file for region descriptors
  //--

  virtual bool Load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) = 0;

  virtual bool Save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const = 0;

  virtual bool LoadFeatures(
    const std::string& sfileNameFeats) = 0;

  //--
  //- Basic description of a descriptor [Type, Length]
  //--
  virtual bool IsScalar() const = 0;
  virtual bool IsBinary() const = 0;

  /// basis element used for description
  virtual std::string Type_id() const = 0;
  virtual size_t DescriptorLength() const = 0;
  //--

  //-- Assume that a region can always be represented at least by a 2D positions
  virtual PointFeatures GetRegionsPositions() const = 0;
  virtual Vec2 GetRegionPosition(size_t i) const = 0;

  /// Return the number of defined regions
  virtual size_t RegionCount() const = 0;

  /// Return a pointer to the first value of the descriptor array
  // Used to avoid complex template imbrication
  virtual const void * DescriptorRawData() const = 0;

  /// Return the squared distance between two descriptors
  // A default metric is used according the descriptor type:
  // - Scalar: L2,
  // - Binary: Hamming
  virtual double SquaredDescriptorDistance(size_t i, const Regions *, size_t j) const = 0;

  /// Add the Inth region to another Region container
  virtual void CopyRegion(size_t i, Regions *) const = 0;

  virtual Regions * EmptyClone() const = 0;

};

/// Specialization of the abstract Regions class to handle Scalar descriptors
template<typename FeatT, typename T, size_t L>
class Scalar_Regions : public Regions
{
public:

  //-- Typedef
  //--

  /// Region type
  typedef FeatT FeatureT;
  /// Region descriptor
  typedef Descriptor<T, L> DescriptorT;

  /// Container for multiple regions
  typedef std::vector<FeatureT> FeatsT;
  /// Container for multiple regions description
  typedef std::vector<DescriptorT > DescsT;

  //-- Class functions
  //--

  bool IsScalar() const {return true;}
  bool IsBinary() const {return false;}
  std::string Type_id() const {return typeid(T).name();}
  size_t DescriptorLength() const {return static_cast<size_t>(L);}

  /// Read from files the regions and their corresponding descriptors.
  bool Load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs)
  {
    return loadFeatsFromFile(sfileNameFeats, _vec_feats)
          & loadDescsFromBinFile(sfileNameDescs, _vec_descs);
  }

  /// Export in two separate files the regions and their corresponding descriptors.
  bool Save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return saveFeatsToFile(sfileNameFeats, _vec_feats)
          & saveDescsToBinFile(sfileNameDescs, _vec_descs);
  }

  bool LoadFeatures(const std::string& sfileNameFeats)
  {
    return loadFeatsFromFile(sfileNameFeats, _vec_feats);
  }

  PointFeatures GetRegionsPositions() const
  {
    return PointFeatures(_vec_feats.begin(), _vec_feats.end());
  }

  Vec2 GetRegionPosition(size_t i) const
  {
    return Vec2f(_vec_feats[i].coords()).cast<double>();
  }

  /// Return the number of defined regions
  size_t RegionCount() const {return _vec_feats.size();}

  /// Mutable and non-mutable FeatureT getters.
  inline std::vector<FeatureT> & Features() { return _vec_feats; }
  inline const std::vector<FeatureT> & Features() const { return _vec_feats; }

  /// Mutable and non-mutable DescriptorT getters.
  inline std::vector<DescriptorT> & Descriptors() { return _vec_descs; }
  inline const std::vector<DescriptorT> & Descriptors() const { return _vec_descs; }

  const void * DescriptorRawData() const { return &_vec_descs[0];}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(_vec_feats);
    ar(_vec_descs);
  }

  Regions * EmptyClone() const
  {
    return new Scalar_Regions();
  }

  // Return the L2 distance between two descriptors
  double SquaredDescriptorDistance(size_t i, const Regions * regions, size_t j) const
  {
    assert(i < _vec_descs.size());
    assert(regions);
    assert(j < regions->RegionCount());

    const Scalar_Regions<FeatT, T, L> * regionsT = dynamic_cast<const Scalar_Regions<FeatT, T, L> *>(regions);
    static matching::L2_Vectorized<T> metric;
    return metric(_vec_descs[i].getData(), regionsT->_vec_descs[j].getData(), DescriptorT::static_size);
  }

  /// Add the Inth region to another Region container
  void CopyRegion(size_t i, Regions * region_container) const
  {
    assert(i < _vec_feats.size() && i < _vec_descs.size());
    static_cast<Scalar_Regions<FeatT, T, L> *>(region_container)->_vec_feats.push_back(_vec_feats[i]);
    static_cast<Scalar_Regions<FeatT, T, L> *>(region_container)->_vec_descs.push_back(_vec_descs[i]);
  }

private:
  //--
  //-- internal data
  std::vector<FeatureT> _vec_feats;    // region features
  std::vector<DescriptorT> _vec_descs; // region descriptions
};

/// Binary_Regions represented as uchar based array
/// It's better to set L as a valid %8 value
///  for memory alignement and performance issue
template<typename FeatT, size_t L>
class Binary_Regions : public Regions
{
public:

  //-- Typedef
  //--

  /// Region
  typedef FeatT FeatureT;
  /// Description of a region
  typedef Descriptor<unsigned char, L> DescriptorT;

  /// Container for multiple regions
  typedef std::vector<FeatureT> FeatsT;
  /// Container for multiple region descriptions
  typedef std::vector<DescriptorT > DescsT;

  //-- Class functions
  //--

  bool IsScalar() const {return false;}
  bool IsBinary() const {return true;}
  std::string Type_id() const {return typeid(typename DescriptorT::bin_type).name();}
  size_t DescriptorLength() const {return static_cast<size_t>(L);}

  /// Read from files the regions and their corresponding descriptors.
  bool Load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs)
  {
    return loadFeatsFromFile(sfileNameFeats, _vec_feats)
          & loadDescsFromBinFile(sfileNameDescs, _vec_descs);
  }

  /// Export in two separate files the regions and their corresponding descriptors.
  bool Save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return saveFeatsToFile(sfileNameFeats, _vec_feats)
          & saveDescsToBinFile(sfileNameDescs, _vec_descs);
  }

  bool LoadFeatures(const std::string& sfileNameFeats)
  {
    return loadFeatsFromFile(sfileNameFeats, _vec_feats);
  }

  PointFeatures GetRegionsPositions() const
  {
    return PointFeatures(_vec_feats.begin(), _vec_feats.end());
  }

  Vec2 GetRegionPosition(size_t i) const
  {
    return Vec2f(_vec_feats[i].coords()).cast<double>();
  }

  /// Return the number of defined regions
  size_t RegionCount() const {return _vec_feats.size();}

  /// Mutable and non-mutable FeatureT getters.
  inline std::vector<FeatureT> & Features() { return _vec_feats; }
  inline const std::vector<FeatureT> & Features() const { return _vec_feats; }

  /// Mutable and non-mutable DescriptorT getters.
  inline std::vector<DescriptorT> & Descriptors() { return _vec_descs; }
  inline const std::vector<DescriptorT> & Descriptors() const { return _vec_descs; }

  const void * DescriptorRawData() const { return &_vec_descs[0];}

  template<class Archive>
  void serialize(Archive & ar)
  {
    ar(_vec_feats);
    ar(_vec_descs);
  }

  Regions * EmptyClone() const
  {
    return new Binary_Regions();
  }

  // Return the squared Hamming distance between two descriptors
  double SquaredDescriptorDistance(size_t i, const Regions * regions, size_t j) const
  {
    assert(i < _vec_descs.size());
    assert(regions);
    assert(j < regions->RegionCount());

    const Binary_Regions<FeatT, L> * regionsT = dynamic_cast<const Binary_Regions<FeatT, L> *>(regions);
    static matching::Hamming<unsigned char> metric;
    const typename matching::Hamming<unsigned char>::ResultType descDist =
      metric(_vec_descs[i].getData(), regionsT->_vec_descs[j].getData(), DescriptorT::static_size);
    return descDist * descDist;
  }

  /// Add the Inth region to another Region container
  void CopyRegion(size_t i, Regions * region_container) const
  {
    assert(i < _vec_feats.size() && i < _vec_descs.size());
    static_cast<Binary_Regions<FeatT, L> *>(region_container)->_vec_feats.push_back(_vec_feats[i]);
    static_cast<Binary_Regions<FeatT, L> *>(region_container)->_vec_descs.push_back(_vec_descs[i]);
  }

private:
  //--
  //-- internal data
  std::vector<FeatureT> _vec_feats; // region features
  std::vector<DescriptorT> _vec_descs; // region descriptions
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_REGIONS_HPP
