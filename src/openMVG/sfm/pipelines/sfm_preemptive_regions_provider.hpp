// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_PREEMPTIVE_REGIONS_PROVIDER_HPP
#define OPENMVG_SFM_SFM_PREEMPTIVE_REGIONS_PROVIDER_HPP

#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"

namespace openMVG {
namespace sfm {

/// Abstract Regions provider
/// Allow to load and return a subset of regions related to a view
struct Preemptive_Regions_Provider : public Regions_Provider
{
public:

  /// Set the number of regions to keep per image (largest scale kept first)
  explicit Preemptive_Regions_Provider(int kept_regions_count = 100):
    Regions_Provider(), kept_regions_count_(kept_regions_count)
  {};

  // Load Regions related to a provided SfM_Data View container
  bool load(
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Regions>& region_type,
    system::ProgressInterface * my_progress_bar = nullptr) override
  {
    if (!my_progress_bar)
      my_progress_bar = &system::ProgressInterface::dummy();
    region_type_.reset(region_type->EmptyClone());

    my_progress_bar->Restart(sfm_data.GetViews().size(), "- Regions ---- Loading -");
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
          OPENMVG_LOG_ERROR << "Invalid regions files for the view: " << sImageName;
          bContinue = false;
        }
        //else
#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          // Sort regions by feature scale & keep the desired count
          regions_ptr->SortAndSelectByRegionScale(kept_regions_count_);
          cache_[iter->second->id_view] = std::move(regions_ptr);
        }
        ++(*my_progress_bar);
      }
    }
    return bContinue;
  }

protected:
  int kept_regions_count_ = 100;
}; // Preemptive_Regions_Provider

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_PREEMPTIVE_REGIONS_PROVIDER_HPP
