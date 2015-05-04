
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_FEATURES_PROVIDER_HPP
#define OPENMVG_SFM_FEATURES_PROVIDER_HPP

#include <openMVG/types.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/features/features.hpp>
#include <memory>

namespace openMVG{

/// Abstract PointFeature provider (read some feature and store them as PointFeature).
/// Allow to load and return the features related to a view
struct Features_Provider
{
  /// PointFeature array per ViewId of the considered SfM_Data container
  Hash_Map<IndexT, PointFeatures> feats_per_view;

  // Load features related to a provided SfM_Data View container
  virtual bool load(
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Image_describer>& image_describer)
  {
    // Read for each view the corresponding features and store them as PointFeatures
    std::unique_ptr<features::Regions> regions;
    image_describer->Allocate(regions);
    for (Views::const_iterator iter = sfm_data.getViews().begin();
      iter != sfm_data.getViews().end(); ++iter)
    {
      const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
      const std::string basename = stlplus::basename_part(sImageName);
      const std::string featFile = stlplus::create_filespec(feat_directory, basename, ".feat");

      if (!image_describer->LoadFeatures(regions.get(), featFile))
      {
        std::cerr << "Invalid feature files for the view: " << sImageName << std::endl;
        return false;
      }
      // save loaded Features as PointFeature
      feats_per_view[iter->second.get()->id_view] = regions->GetRegionsPositions();
    }
    return true;
  }

  virtual bool load(
    const SfM_Data & sfm_data,
    const std::string & feat_directory,
    std::unique_ptr<features::Regions>& region_type)
  {
    // Read for each view the corresponding features and store them as PointFeatures
    std::unique_ptr<features::Regions> regions(region_type->EmptyClone());
    for (Views::const_iterator iter = sfm_data.getViews().begin();
      iter != sfm_data.getViews().end(); ++iter)
    {
      const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
      const std::string basename = stlplus::basename_part(sImageName);
      const std::string featFile = stlplus::create_filespec(feat_directory, basename, ".feat");

      if (!regions->LoadFeatures(featFile))
      {
        std::cerr << "Invalid feature files for the view: " << sImageName << std::endl;
        return false;
      }
      // save loaded Features as PointFeature
      feats_per_view[iter->second.get()->id_view] = regions->GetRegionsPositions();
    }
    return true;
  }


  /// Return the PointFeatures belonging to the View, if the view does not exist
  ///  it return an empty PointFeature array.
  const PointFeatures & getFeatures(const IndexT & id_view) const
  {
    // Have an empty feature set in order to deal with non existing view_id
    static const PointFeatures emptyFeats = PointFeatures();

    Hash_Map<IndexT, PointFeatures>::const_iterator it = feats_per_view.find(id_view);
    if (it != feats_per_view.end())
      return it->second;
    else
      return emptyFeats;
  }
}; // Features_Provider

} // namespace openMVG

#endif // OPENMVG_SFM_FEATURES_PROVIDER_HPP
