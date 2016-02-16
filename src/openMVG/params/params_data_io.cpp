
// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/params/params_data_io.hpp"
#include "openMVG/params/params_data_io_cereal.hpp"

#include "openMVG/stl/stlMap.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"


namespace openMVG {
namespace params {

bool Load(paramsData & params_data, const std::string & filename)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    bStatus = Load_Cereal<cereal::JSONInputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown params_data input format: " << ext << std::endl;
    return false;
  }

  return bStatus;
}

bool Save(const paramsData & params_data, const std::string & filename)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown sfm_data export format: " << ext << std::endl;
  }
  return false;
}

} // namespace params
} // namespace openMVG

