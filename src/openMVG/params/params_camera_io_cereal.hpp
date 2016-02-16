// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_PARAMS_CAMERA_IO_CEREAL_HPP
#define OPENMVG_PARAMS_CAMERA_IO_CEREAL_HPP

#include "openMVG/params/params_camera_io.hpp"

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
namespace params {

template <
// JSONInputArchive/ ...
typename archiveType
>
bool Load_Cereal(
  paramsCamera & data,
  const std::string & filename)
{
  const bool bBinary = stlplus::extension_part(filename) == "bin";

  //Create the stream and check it is ok
  std::ifstream stream(filename.c_str(), std::ios::binary | std::ios::in);
  if (!stream.is_open())
    return false;

  // Data serialization
  try
  {
    archiveType archive(stream);

    std::string version;
    archive(cereal::make_nvp("camera_data_version", version));
    archive(cereal::make_nvp("focal_px", data.focal_px));
    archive(cereal::make_nvp("kMatrix", data.kMatrix));
    archive(cereal::make_nvp("camera_type", data.camera_type));
    archive(cereal::make_nvp("shared_intrinsics", data.shared_intrinsics));
  }
  catch (const cereal::Exception & e)
  {
    std::cerr << e.what() << std::endl;
    return false;
  }
  stream.close();
  return true;
}

template <
// JSONOutputArchive/ ...
typename archiveType
>
bool Save_Cereal(
  const paramsCamera & data,
  const std::string & filename)
{
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
    archive(cereal::make_nvp("camera_data_version", version));
    archive(cereal::make_nvp("focal_px", data.focal_px));
    archive(cereal::make_nvp("kMatrix", data.kMatrix));
    archive(cereal::make_nvp("camera_type", data.camera_type));
    archive(cereal::make_nvp("shared_intrinsics", data.shared_intrinsics));
  }
  stream.close();
  return true;
}

} // namespace params
} // namespace openMVG

#endif // OPENMVG_PARAMS_CAMERA_IO_CEREAL_HPP
