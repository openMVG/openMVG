// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_VIEW_METADATA_HPP
#define OPENMVG_SFM_VIEW_METADATA_HPP

#include "openMVG/sfm/sfm_view.hpp"

namespace openMVG {
namespace sfm {

/// A view with all the input image metadata
struct View_Metadata : public View
{
  // Map for metadata
  std::map<std::string, std::string> metadata;

  // Constructor (use unique index for the view_id)
  View_Metadata(
    const std::string & sImgPath = "",
    IndexT view_id = UndefinedIndexT,
    IndexT intrinsic_id = UndefinedIndexT,
    IndexT pose_id = UndefinedIndexT,
    IndexT width = UndefinedIndexT, IndexT height = UndefinedIndexT,
    const std::map<std::string, std::string>& metadata = std::map<std::string, std::string>())
    : View(sImgPath, view_id, intrinsic_id, pose_id, width, height), 
    metadata(metadata)
    {
    }

  // Serialization
  template <class Archive>
  void serialize(Archive & ar)
  {
    View::serialize(ar);
    ar(cereal::make_nvp("metadata", metadata));
  }
};

} // namespace sfm
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::sfm::View_Metadata, "view_metadata");

#endif // OPENMVG_SFM_VIEW_METADATA_HPP
