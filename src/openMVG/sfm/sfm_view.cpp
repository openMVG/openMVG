// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_view.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

namespace openMVG {
namespace sfm {

View::View(
  const std::string & sImgPath,
  IndexT view_id,
  IndexT intrinsic_id,
  IndexT pose_id,
  IndexT width, 
  IndexT height)
  : s_Img_path(sImgPath)
  , id_view(view_id)
  , id_intrinsic(intrinsic_id)
  , id_pose(pose_id)
  , ui_width(width)
  , ui_height(height)
  {
  }

template <class Archive>
void View::save( Archive & ar ) const
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
void View::load( Archive & ar )
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

template void View::load
<cereal::JSONInputArchive>
( cereal::JSONInputArchive & ar );

template void View::save
<cereal::JSONOutputArchive>
(cereal::JSONOutputArchive & ar) const;

template void View::load
<cereal::XMLInputArchive>
( cereal::XMLInputArchive & ar );

template void View::save
<cereal::XMLOutputArchive>
(cereal::XMLOutputArchive & ar) const;

template void View::load
<cereal::BinaryInputArchive>
( cereal::BinaryInputArchive & ar );

template void View::save
<cereal::BinaryOutputArchive>
(cereal::BinaryOutputArchive & ar) const;

template void View::load
<cereal::PortableBinaryInputArchive>
( cereal::PortableBinaryInputArchive & ar );

template void View::save
<cereal::PortableBinaryOutputArchive>
(cereal::PortableBinaryOutputArchive & ar) const;


} // namespace sfm
} // namespace openMVG
