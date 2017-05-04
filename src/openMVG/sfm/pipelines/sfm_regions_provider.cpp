// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG {
namespace sfm {

std::string Regions_Provider::Type_id()
{
  if (region_type_)
    return region_type_->Type_id();
  else
  {
    return std::string("initialized regions type");
  }
}

bool Regions_Provider::IsScalar()
{
  if (region_type_)
    return region_type_->IsScalar();
  else
  {
    std::cerr << "Invalid region type" << std::endl;
    return false;
  }
}

bool Regions_Provider::IsBinary()
{
  if (region_type_)
    return region_type_->IsBinary();
  else
  {
    std::cerr << "Invalid region type" << std::endl;
    return false;
  }
}

const openMVG::features::Regions* Regions_Provider::getRegionsType() const
{
  if (region_type_)
    return &(*region_type_);
  else
    return nullptr;
}

std::shared_ptr<features::Regions> Regions_Provider::get(const IndexT x) const
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
bool Regions_Provider::load(
  const SfM_Data & sfm_data,
  const std::string & feat_directory,
  std::unique_ptr<features::Regions>& region_type,
  C_Progress *  my_progress_bar)
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


} // namespace sfm
} // namespace openMVG
