// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_REGIONS_PROVIDER_HPP
#define OPENMVG_SFM_SFM_REGIONS_PROVIDER_HPP

#include <atomic>
#include <memory>
#include <string>

#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace sfm {

/// Abstract Regions provider
/// Allow to load and return the regions related to a view
struct Regions_Provider
{
public:

  virtual ~Regions_Provider() = default;

  std::string Type_id()
  {
    if (region_type_)
      return region_type_->Type_id();
    else
    {
      return std::string("initialized regions type");
    }
  }

  bool IsScalar()
  {
    if (region_type_)
      return region_type_->IsScalar();
    else
    {
      std::cerr << "Invalid region type" << std::endl;
      return false;
    }
  }

  bool IsBinary()
  {
    if (region_type_)
      return region_type_->IsBinary();
    else
    {
      std::cerr << "Invalid region type" << std::endl;
      return false;
    }
  }

  const openMVG::features::Regions* getRegionsType() const
  {
    if (region_type_)
      return &(*region_type_);
    else
      return nullptr;
  }

  virtual std::shared_ptr<features::Regions> get(const IndexT x) const
  {
    auto it = cache_.find(x);
    std::shared_ptr<features::Regions> ret;

    if (it != end(cache_))
    {
      ret = it->second;
    }
    // else Invalid ressource
    return ret;
  }

  // Load Regions related to a provided SfM_Data View container
  virtual bool load(
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Regions>& region_type,
    C_Progress *  my_progress_bar = nullptr)
  {
    if (!my_progress_bar)
      my_progress_bar = &C_Progress::dummy();
    region_type_.reset(region_type->EmptyClone());

    my_progress_bar->restart(sfm_data.GetViews().size(), "\n- Regions Loading -\n");
    // Read for each view the corresponding regions and store them
    std::atomic<bool> bContinue(true);
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel
#endif
    for (Views::const_iterator iter = sfm_data.GetViews().begin();
      iter != sfm_data.GetViews().end() && bContinue; ++iter)
    {
        if (my_progress_bar->hasBeenCanceled())
        {
          bContinue = false;
          continue;
        }
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif
      {
        const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second->s_Img_path);
        const std::string basename = stlplus::basename_part(sImageName);
        const std::string featFile = stlplus::create_filespec(feat_directory, basename, ".feat");
        const std::string descFile = stlplus::create_filespec(feat_directory, basename, ".desc");

        std::unique_ptr<features::Regions> regions_ptr(region_type->EmptyClone());
        if (!regions_ptr->Load(featFile, descFile))
        {
          std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
          bContinue = false;
        }
        //else
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          cache_[iter->second->id_view] = std::move(regions_ptr);
        }
        ++(*my_progress_bar);
      }
    }
    return bContinue;
  }

protected:
  /// Regions per ViewId of the considered SfM_Data container
  mutable Hash_Map<IndexT, std::shared_ptr<features::Regions> > cache_;
  std::unique_ptr<openMVG::features::Regions> region_type_;
}; // Regions_Provider

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_REGIONS_PROVIDER_HPP
