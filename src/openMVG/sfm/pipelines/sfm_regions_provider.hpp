
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_REGIONS_PROVIDER_HPP
#define OPENMVG_SFM_REGIONS_PROVIDER_HPP

#include <openMVG/types.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/features/regions.hpp>
#include <memory>

namespace openMVG{

/// Abstract Regions provider
/// Allow to load and return the regions related to a view
struct Regions_Provider
{
  /// Regions per ViewId of the considered SfM_Data container
  Hash_Map<IndexT, std::unique_ptr<features::Regions> > regions_per_view;

  // Load Regions related to a provided SfM_Data View container
  virtual bool load(
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Image_describer>& image_describer)
  {
    // Read for each view the corresponding regions and store them
    for (Views::const_iterator iter = sfm_data.getViews().begin();
      iter != sfm_data.getViews().end(); ++iter)
    {
      const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
      const std::string basename = stlplus::basename_part(sImageName);
      const std::string featFile = stlplus::create_filespec(feat_directory, basename, ".feat");
      const std::string descFile = stlplus::create_filespec(feat_directory, basename, ".desc");

      image_describer->Allocate(regions_per_view[iter->second.get()->id_view]);
      if (!image_describer->Load(
        regions_per_view[iter->second.get()->id_view].get(),
        featFile,
        descFile))
      {
        std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
        return false;
      }
    }
    return true;
  }

}; // Regions_Provider

} // namespace openMVG

#endif // OPENMVG_SFM_REGIONS_PROVIDER_HPP
