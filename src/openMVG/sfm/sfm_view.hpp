// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_VIEW_HPP
#define OPENMVG_SFM_VIEW_HPP

#include "openMVG/types.hpp"
#include "openMVG/numeric/numeric.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>

namespace openMVG {
namespace sfm {

/// A view define an image by a string and unique indexes for the view, the camera intrinsic & the pose
struct View
{
  // image path on disk
  std::string s_Img_path;

  // Id of the view
  IndexT id_view;

  // Index of intrinsics and the pose
  IndexT id_intrinsic, id_pose;

  // image size
  IndexT ui_width, ui_height;

  // Constructor (use unique index for the view_id)
  View(
    const std::string & sImgPath = "",
    IndexT view_id = UndefinedIndexT,
    IndexT intrinsic_id = UndefinedIndexT,
    IndexT pose_id = UndefinedIndexT,
    IndexT width = UndefinedIndexT, IndexT height = UndefinedIndexT)
    :s_Img_path(sImgPath), id_view(view_id), id_intrinsic(intrinsic_id),
    id_pose(pose_id), ui_width(width), ui_height(height)
    {}

  virtual ~View() = default;

  /**
  * Serialization out
  * @param ar Archive
  */
  template <class Archive>
  void save( Archive & ar ) const
  {
    ar(cereal::make_nvp("local_path", stlplus::folder_part(s_Img_path)),
       cereal::make_nvp("filename", stlplus::filename_part(s_Img_path)),
       cereal::make_nvp("width", ui_width),
       cereal::make_nvp("height", ui_height),
       cereal::make_nvp("id_view", id_view),
       cereal::make_nvp("id_intrinsic", id_intrinsic),
       cereal::make_nvp("id_pose", id_pose));
  }

  /**
  * @brief Serialization in
  * @param ar Archive
  */
  template <class Archive>
  void load( Archive & ar )
  {
    //Define a view with two string (base_path & basename)
    std::string local_path = s_Img_path;
    std::string filename = s_Img_path;

    ar(cereal::make_nvp("local_path", local_path),
       cereal::make_nvp("filename", filename),
       cereal::make_nvp("width", ui_width),
       cereal::make_nvp("height", ui_height),
       cereal::make_nvp("id_view", id_view),
       cereal::make_nvp("id_intrinsic", id_intrinsic),
       cereal::make_nvp("id_pose", id_pose));

    s_Img_path = stlplus::create_filespec(local_path, filename);
  }
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_VIEW_HPP
