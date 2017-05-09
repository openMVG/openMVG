// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG {
namespace sfm {

Regions_Provider_Cache::Regions_Provider_Cache (const unsigned int max_cache_size)
  : Regions_Provider(), max_cache_size_(max_cache_size)
{
}

std::shared_ptr<features::Regions> Regions_Provider_Cache::get(const IndexT x) const
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
bool Regions_Provider_Cache::load
(
  const SfM_Data & sfm_data,
  const std::string & feat_directory,
  std::unique_ptr<features::Regions>& region_type,
  C_Progress * 
)
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

std::size_t Regions_Provider_Cache::prune(bool bPruneOnlyOne) const
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

std::size_t Regions_Provider_Cache::size() const { return cache_.size(); }


} // namespace sfm
} // namespace openMVG
