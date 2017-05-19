// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_REGIONS_HPP
#define OPENMVG_FEATURES_REGIONS_HPP

#include <string>
#include <openMVG/features/feature.hpp>
#include <openMVG/features/feature_container.hpp>
#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace features {

/// Describe an image a set of regions (position, ...) + attributes
/// Each region is described by a set of attributes (descriptor)
class Regions
{
public:

  virtual ~Regions() = default;

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
  // - Scalar: SquaredL2,
  // - Binary: SquaredHamming
  virtual double SquaredDescriptorDistance(
    size_t i,
    const Regions *,
    size_t j) const = 0;

  /// Add the Inth region to another Region container
  virtual void CopyRegion(size_t i, Regions *) const = 0;

  virtual Regions * EmptyClone() const = 0;

};

std::unique_ptr<features::Regions> Init_region_type_from_file
(
  const std::string & sImage_describer_file
);

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_REGIONS_HPP
