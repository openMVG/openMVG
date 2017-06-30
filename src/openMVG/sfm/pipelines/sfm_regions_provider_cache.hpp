// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_REGIONS_PROVIDER_CACHE_HPP
#define OPENMVG_SFM_SFM_REGIONS_PROVIDER_CACHE_HPP

#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <string>
#include <thread>

namespace openMVG {
namespace sfm {

/// Regions provider Cache
/// Store only a given count of regions in memory
struct Regions_Provider_Cache : public Regions_Provider
{
public:

  Regions_Provider_Cache (const unsigned int max_cache_size)
    : Regions_Provider(), max_cache_size_(max_cache_size)
  {
  }

  std::shared_ptr<features::Regions> get(const IndexT x) const override
  {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = cache_.find(x);
    std::shared_ptr<features::Regions> ret;

    if (it == end(cache_))
    {
      // Load the ressource link to this ID
      const std::string id =
        stlplus::create_filespec(feat_directory_, map_id_string_.at(x));
      const std::string featFile = id + ".feat";
      const std::string descFile = id + ".desc";
      ret.reset(region_type_->EmptyClone());
      if (ret->Load(featFile, descFile))
      {
        cache_[x] = ret;
      }
      else
      {
        // Invalid ressource
      }
    }
    else
    {
      ret = it->second;
    }
    // If the cache is too large:
    //  - try to prune elements that are no longer used
    if (size() > max_cache_size_)
    {
      prune(true); // Prune only one cached object
    }
    return ret;
  }

  // Initialize the regions_provider_cache
  bool load
  (
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Regions>& region_type,
    C_Progress * 
  ) override
  {
    std::cout << "Initialization of the Regions_Provider_Cache. #Elements in the cache: "<< max_cache_size_ << std::endl;

    feat_directory_ = feat_directory;
    region_type_.reset(region_type->EmptyClone());

    // Build an association table from view id to feature & descriptor files
    for (const auto & iterViews : sfm_data.GetViews())
    {
      const openMVG::IndexT id = iterViews.second->id_view;
      assert( id == iterViews.first);
      map_id_string_[id] = stlplus::basename_part(iterViews.second->s_Img_path);
    }

    return true;
  }

private:

  mutable std::mutex mutex_; // To deal with multithread concurrent access

  std::string feat_directory_; // The regions file directory
  std::map<openMVG::IndexT, std::string> map_id_string_; // association of the view id & its basename
  const unsigned int max_cache_size_;

private:

  /// @brief Prune smart_ptr that are only referenced into the cache (not longer used externally)
  /// @return the number of removed elements
  std::size_t prune(bool bPruneOnlyOne = false) const
  {
    std::size_t count = 0;
    for (auto it = begin(cache_); it != end(cache_);)
    {
      if (it->second.use_count() == 1)
      {
        cache_.erase(it++);
        ++count;
        if (bPruneOnlyOne)
        {
          break;
        }
      }
      else
      {
        ++it;
      }
    }
    return count;
  }

  std::size_t size() const { return cache_.size(); }

}; // Regions_Provider_Cache

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_REGIONS_PROVIDER_CACHE_HPP
