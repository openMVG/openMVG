// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_REGIONS_PROVIDER_CACHE_HPP
#define OPENMVG_SFM_SFM_REGIONS_PROVIDER_CACHE_HPP

#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

#include <mutex>
#include <map>
#include <string>

namespace openMVG {
namespace sfm {

/// Regions provider Cache
/// Store only a given count of regions in memory
struct Regions_Provider_Cache : public Regions_Provider
{
public:

  Regions_Provider_Cache (const unsigned int max_cache_size);

  std::shared_ptr<features::Regions> get(const IndexT x) const override;

  // Initialize the regions_provider_cache
  bool load
  (
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Regions>& region_type,
    C_Progress * 
  ) override;

private:

  mutable std::mutex mutex_; // To deal with multithread concurrent access

  std::string feat_directory_; // The regions file directory
  std::map<openMVG::IndexT, std::string> map_id_string_; // association of the view id & its basename
  const unsigned int max_cache_size_;

private:

  /// @brief Prune smart_ptr that are only referenced into the cache (not longer used externally)
  /// @return the number of removed elements
  std::size_t prune(bool bPruneOnlyOne = false) const;

  std::size_t size() const;

}; // Regions_Provider_Cache

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_REGIONS_PROVIDER_CACHE_HPP
