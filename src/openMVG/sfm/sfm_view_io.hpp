// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_VIEW_IO_HPP
#define OPENMVG_SFM_SFM_VIEW_IO_HPP

#include "openMVG/sfm/sfm_view.hpp"

#include <cereal/types/polymorphic.hpp>

template <class Archive>
void openMVG::sfm::View::save( Archive & ar ) const
{
  ar(cereal::make_nvp("local_path", stlplus::folder_part(s_Img_path)),
     cereal::make_nvp("filename", stlplus::filename_part(s_Img_path)),
     cereal::make_nvp("width", ui_width),
     cereal::make_nvp("height", ui_height),
     cereal::make_nvp("id_view", id_view),
     cereal::make_nvp("id_intrinsic", id_intrinsic),
     cereal::make_nvp("id_pose", id_pose));
}


template <class Archive>
void openMVG::sfm::View::load( Archive & ar )
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


#endif // OPENMVG_SFM_SFM_VIEW_IO_HPP