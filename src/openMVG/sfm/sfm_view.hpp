// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_VIEW_HPP
#define OPENMVG_SFM_VIEW_HPP

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG {

/// A view define an image by a string and unique indexes for the view, the camera intrinsic & the pose
struct View
{
  std::string s_Img_path;

  // Id of the view
  IndexT id_view;

  // Index of intrinsics and the pose
  IndexT id_intrinsic, id_pose;

  // Serialization
  template <class Archive>
  void serialize( Archive & ar )
  {
    //Define a view with two string (base_path & basename)
    std::string local_path = stlplus::folder_append_separator(
      stlplus::folder_part(s_Img_path));
    std::string filename = stlplus::filename_part(s_Img_path);

    ar(cereal::make_nvp("local_path", local_path),
       cereal::make_nvp("filename", filename),
       cereal::make_nvp("id_view", id_view),
       cereal::make_nvp("id_intrinsic", id_intrinsic),
       cereal::make_nvp("id_pose", id_pose));

    s_Img_path = stlplus::create_filespec(local_path, filename);
  }
};

} // namespace openMVG

#endif // OPENMVG_SFM_VIEW_HPP
