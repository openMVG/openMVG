
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_KEYPOINTSET_HPP
#define OPENMVG_FEATURES_KEYPOINTSET_HPP

#include "openMVG/features/feature.hpp"
#include "openMVG/features/descriptor.hpp"
#include <string>

namespace openMVG {
namespace features {

/// Association storage of associated feature and descriptor for a given image.
/// Load, save, R/W accessor operation.
///
/// typedef vector<SIOPointFeature> featsT;
/// typedef vector<Descriptor<uchar,128> > descsT;
/// KeypointSet< featsT, descsT > kpSet;
template<typename FeaturesT, typename DescriptorsT>
class KeypointSet {
public:
  // Alias to stored Feature and Descriptor type
  typedef typename FeaturesT::value_type FeatureT;
  typedef typename DescriptorsT::value_type DescriptorT;

  /// Read from files the feats and their corresponding descriptors.
  bool loadFromFile(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs)
  {
    return loadFeatsFromFile(sfileNameFeats, feats_)
          & loadDescsFromFile(sfileNameDescs, descs_);
  }

  /// Export in two separate files the feats and their corresponding descriptors.
  bool saveToFile(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return saveFeatsToFile(sfileNameFeats, feats_)
          & saveDescsToFile(sfileNameDescs, descs_);
  }

  /// Read from files the feats and their corresponding descriptors
  ///  descriptor in binary to save place
  bool loadFromBinFile(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs)
  {
    return loadFeatsFromFile(sfileNameFeats, feats_)
          & loadDescsFromBinFile(sfileNameDescs, descs_);
  }

  /// Export in two separate files the feats and their corresponding descriptors
  ///  descriptor in binary to save place
  bool saveToBinFile(
    const std::string& sfileNameFeats,
    const std::string& sfileNameDescs) const
  {
    return saveFeatsToFile(sfileNameFeats, feats_)
          & saveDescsToBinFile(sfileNameDescs, descs_);
  }

  /// Mutable and non-mutable FeatureT getters.
  inline FeaturesT & features() { return feats_; }
  inline const FeaturesT & features() const { return feats_; }

  /// Mutable and non-mutable DescriptorT getters.
  inline DescriptorsT & descriptors() { return descs_; }
  inline const DescriptorsT & descriptors() const { return descs_; }

private:
  FeaturesT feats_;
  DescriptorsT descs_;
};

} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_KEYPOINTSET_HPP
