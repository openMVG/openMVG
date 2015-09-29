
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SFM_DATA_IO_CEREAL_HPP
#define OPENMVG_SFM_DATA_IO_CEREAL_HPP

#include "openMVG/sfm/sfm_data_io.hpp"

#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

#include <cereal/types/map.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include <iomanip>
#include <fstream>

namespace openMVG {
namespace sfm {

template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  SfM_Data & data,
  const std::string & filename,
  ESfM_Data flags_part)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  // List which part of the file must be considered
  const bool b_views = (flags_part & VIEWS) == VIEWS;
  const bool b_intrinsics = (flags_part & INTRINSICS) == INTRINSICS;
  const bool b_extrinsics = (flags_part & EXTRINSICS) == EXTRINSICS;
  const bool b_structure = (flags_part & STRUCTURE) == STRUCTURE;
  const bool b_control_point = (flags_part & CONTROL_POINTS) == CONTROL_POINTS;

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("sfm_data_version", version));
    archive(cereal::make_nvp("root_path", data.s_root_path));

    if (b_views)
      archive(cereal::make_nvp("views", data.views));
    else
      if (bBinary) { // Binary file require read all the member
        Views views;
        archive(cereal::make_nvp("views", views));
      }

    if (b_intrinsics)
      archive(cereal::make_nvp("intrinsics", data.intrinsics));
    else
      if (bBinary) { // Binary file require read all the member
        Intrinsics intrinsics;
        archive(cereal::make_nvp("intrinsics", intrinsics));
      }

    if (b_extrinsics)
      archive(cereal::make_nvp("extrinsics", data.poses));
    else
      if (bBinary) { // Binary file require read all the member
        Poses poses;
        archive(cereal::make_nvp("extrinsics", poses));
      }

    if (b_structure)
      archive(cereal::make_nvp("structure", data.structure));
    else
      if (bBinary) { // Binary file require read all the member
        Landmarks structure;
        archive(cereal::make_nvp("structure", structure));
      }

    if (version != "0.1") // fast check to assert we are at least using version 0.2
    {
      if (b_control_point)
        archive(cereal::make_nvp("control_points", data.control_points));
      else
        if (bBinary) { // Binary file require read all the member
          Landmarks control_points;
          archive(cereal::make_nvp("control_points", control_points));
        }
    }
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const SfM_Data & data,
  const std::string & filename,
  ESfM_Data flags_part)
{
  // List which part of the file must be considered
  const bool b_views = (flags_part & VIEWS) == VIEWS;
  const bool b_intrinsics = (flags_part & INTRINSICS) == INTRINSICS;
  const bool b_extrinsics = (flags_part & EXTRINSICS) == EXTRINSICS;
  const bool b_structure = (flags_part & STRUCTURE) == STRUCTURE;
  const bool b_control_point = (flags_part & CONTROL_POINTS) == CONTROL_POINTS;

  //Create the stream and check it is ok
  std::ofstream stream(filename.c_str(), std::ios::binary | std::ios::out);
  if (!stream.is_open())
    return false;

  // Data serialization
  {
    archiveType archive(stream);
    // since OpenMVG 0.9, the sfm_data version 0.2 is introduced
    //  - it adds control_points storage
    const std::string version = "0.2";
    archive(cereal::make_nvp("sfm_data_version", version));
    archive(cereal::make_nvp("root_path", data.s_root_path));

    if (b_views)
      archive(cereal::make_nvp("views", data.views));
    else
      archive(cereal::make_nvp("views", Views()));

    if (b_intrinsics)
      archive(cereal::make_nvp("intrinsics", data.intrinsics));
    else
      archive(cereal::make_nvp("intrinsics", Intrinsics()));

    if (b_extrinsics)
      archive(cereal::make_nvp("extrinsics", data.poses));
    else
      archive(cereal::make_nvp("extrinsics", Poses()));

    // Structure -> See for export in another file
    if (b_structure)
      archive(cereal::make_nvp("structure", data.structure));
    else
      archive(cereal::make_nvp("structure", Landmarks()));

    if (version != "0.1") // fast check to assert we are at least using version 0.2
    {
      if (b_control_point)
        archive(cereal::make_nvp("control_points", data.control_points));
      else
        archive(cereal::make_nvp("control_points", Landmarks()));
    }
  }
  return true;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_IO_CEREAL_HPP
