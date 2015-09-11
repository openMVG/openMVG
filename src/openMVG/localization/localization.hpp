
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "reconstructed_regions.hpp"

#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/stl/stlMap.hpp>


namespace openMVG {
namespace localization {

typedef Reconstructed_Regions<features::SIOPointFeature, unsigned char, 128> Reconstructed_RegionsT;

class VoctreeLocalizer
{
public:
  Hash_Map<IndexT, Reconstructed_RegionsT > regions_per_view;

  // loadSfmData(const std::string & sfmDataPath)

  /**
   * @brief Load all the Descriptors who have contributed to the reconstruction.
   */
  bool loadReconstructionDescriptors(
    const sfm::SfM_Data & sfm_data,
    const std::string & feat_directory);

  /**
   * @brief Load the vocabulary tree.
   */
  bool loadVoctree(const std::string & feat_directory) {}

  
};

}
}
