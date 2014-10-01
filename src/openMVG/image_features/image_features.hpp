// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/features/features.hpp"

#include <bitset>
#include <string>
#include <typeinfo>

namespace openMVG{
namespace image_features{

/// An image is described by a set of regions (position...)
/// and each regions is described by a set of attributes (descriptor)
class Regions
{
public:

  //--
  // IO
  //--

  virtual bool load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) = 0;

  virtual bool save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const = 0;

  //--
  //- Basic description of a descriptor [Type, Length]
  //--
  virtual bool isScalar() const = 0;
  virtual bool isBinary() const = 0;

  virtual std::string type_id() const = 0; // basis element used for description
  virtual size_t descriptorLength() const = 0;
  //--

  //-- Assume that a region can alway be represented at least as a position
  virtual std::vector<PointFeature> getRegionsPositions() const = 0;

};

template<typename T, size_t L>
class Regions_ScalarDescriptorGeneric : public Regions
{
public:

  //-- Typedef
  //--

  // Region
  typedef SIOPointFeature FeatureT;
   // Description of a region
  typedef Descriptor<T, L> DescriptorT;

  // Container for multiple regions
  typedef std::vector<FeatureT> FeatsT;
  // Container for multiple regions description
  typedef std::vector<DescriptorT > DescsT;

  //-- Class functions
  //--

  bool isScalar() const {return true;}
  bool isBinary() const {return false;}
  std::string type_id() const {return typeid(T).name();}
  size_t descriptorLength() const {return static_cast<size_t>(L);}

  /// Read from files the regions and their corresponding descriptors.
  bool load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs)
  {
    return loadFeatsFromFile(sfileNameFeats, _vec_feats)
          & loadDescsFromFile(sfileNameDescs, _vec_descs);
  }

  /// Export in two separate files the regions and their corresponding descriptors.
  bool save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return saveFeatsToFile(sfileNameFeats, _vec_feats)
          & saveDescsToFile(sfileNameDescs, _vec_descs);
  }

  std::vector<PointFeature> getRegionsPositions() const
  {
    std::vector<PointFeature> pts(_vec_feats.size());
    for(size_t i = 0; i < _vec_feats.size(); ++i)
      pts[i]= PointFeature(_vec_feats[i].x(), _vec_feats[i].y());
    return pts;
  }

  /// Mutable and non-mutable FeatureT getters.
  inline std::vector<FeatureT> & features() { return _vec_feats; }
  inline const std::vector<FeatureT> & features() const { return _vec_feats; }

  /// Mutable and non-mutable DescriptorT getters.
  inline std::vector<DescriptorT> & descriptors() { return _vec_descs; }
  inline const std::vector<DescriptorT> & descriptors() const { return _vec_descs; }

private:
  //--
  //-- internal data
  std::vector<FeatureT> _vec_feats; // regions
  std::vector<DescriptorT> _vec_descs; // regions description
};

template<size_t L>
class Regions_BinaryDescriptorGeneric : public Regions
{
public:

  //-- Typedef
  //--

  // Region
  typedef SIOPointFeature FeatureT;
   // Description of a region
  typedef Descriptor<std::bitset<L>, 1> DescriptorT;

  // Container for multiple regions
  typedef std::vector<FeatureT> FeatsT;
  // Container for multiple regions description
  typedef std::vector<DescriptorT > DescsT;

  //-- Class functions
  //--

  bool isScalar() const {return false;}
  bool isBinary() const {return true;}
  std::string type_id() const {return typeid(typename DescriptorT::bin_type).name();}
  size_t descriptorLength() const {return static_cast<size_t>(L);}

  /// Read from files the regions and their corresponding descriptors.
  bool load(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs)
  {
    return loadFeatsFromFile(sfileNameFeats, _vec_feats)
          & loadDescsFromFile(sfileNameDescs, _vec_descs);
  }

  /// Export in two separate files the regions and their corresponding descriptors.
  bool save(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return saveFeatsToFile(sfileNameFeats, _vec_feats)
          & saveDescsToFile(sfileNameDescs, _vec_descs);
  }

  std::vector<PointFeature> getRegionsPositions() const
  {
    std::vector<PointFeature> pts(_vec_feats.size());
    for(size_t i = 0; i < _vec_feats.size(); ++i)
      pts[i] = PointFeature(_vec_feats[i].x(), _vec_feats[i].y());
    return pts;
  }

  /// Mutable and non-mutable FeatureT getters.
  inline std::vector<FeatureT> & features() { return _vec_feats; }
  inline const std::vector<FeatureT> & features() const { return _vec_feats; }

  /// Mutable and non-mutable DescriptorT getters.
  inline std::vector<DescriptorT> & descriptors() { return _vec_descs; }
  inline const std::vector<DescriptorT> & descriptors() const { return _vec_descs; }

private:
  //--
  //-- internal data
  std::vector<FeatureT> _vec_feats; // regions
  std::vector<DescriptorT> _vec_descs; // regions description
};

typedef Regions_ScalarDescriptorGeneric<unsigned char,128> SIFT_Regions;
typedef Regions_ScalarDescriptorGeneric<float,20> DISSOCIATED_DIPOLES_Regions;

} // image_features
} // namespace openMVG

