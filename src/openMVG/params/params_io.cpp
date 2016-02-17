
// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2015 Klemen Istenic (CIRS - Universitat de Girona).

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/params/params_io.hpp>
#include <openMVG/params/params_io_cereal.hpp>
#include "openMVG/stl/stlMap.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"


namespace openMVG {
namespace params {


// --------------------
// Camera data
// --------------------
bool Load(paramsCamera & params_data, const std::string & filename)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json"){
    bStatus = Load_Cereal<cereal::JSONInputArchive>(params_data, filename);
  	params_data.valid = true;
  }
  else
  {
    std::cerr << "Unknown params_camera input format: " << ext << std::endl;
    return false;
  }

  return bStatus;
}

bool Save(const paramsCamera & params_data, const std::string & filename)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown params_camera export format: " << ext << std::endl;
  }
  return false;
}


// --------------------
// Detection data
// --------------------
bool Load(paramsDetection & params_data, const std::string & filename)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json"){
    bStatus = Load_Cereal<cereal::JSONInputArchive>(params_data, filename);
  	params_data.valid = true;
  }
  else
  {
    std::cerr << "Unknown params_detection input format: " << ext << std::endl;
    return false;
  }

  return bStatus;
}

bool Save(const paramsDetection & params_data, const std::string & filename)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown params_detection export format: " << ext << std::endl;
  }
  return false;
}


// --------------------
// Matching data
// --------------------
bool Load(paramsMatching & params_data, const std::string & filename)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json"){
    bStatus = Load_Cereal<cereal::JSONInputArchive>(params_data, filename);
  	params_data.valid = true;
  }
  else
  {
    std::cerr << "Unknown params_matching input format: " << ext << std::endl;
    return false;
  }

  return bStatus;
}

bool Save(const paramsMatching & params_data, const std::string & filename)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown params_matching export format: " << ext << std::endl;
  }
  return false;
}


// --------------------
// Incremental SfM data
// --------------------

bool Load(paramsIncrementalSfM & params_data, const std::string & filename)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json"){
    bStatus = Load_Cereal<cereal::JSONInputArchive>(params_data, filename);
  	params_data.valid = true;
  }
  else
  {
    std::cerr << "Unknown params_incrementalSfM input format: " << ext << std::endl;
    return false;
  }

  return bStatus;
}

bool Save(const paramsIncrementalSfM & params_data, const std::string & filename)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown params_incrementalSfM export format: " << ext << std::endl;
  }
  return false;
}


// --------------------
// Global SfM data
// --------------------

bool Load(paramsGlobalSfM & params_data, const std::string & filename)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json"){
    bStatus = Load_Cereal<cereal::JSONInputArchive>(params_data, filename);
  	params_data.valid = true;
  }
  else
  {
    std::cerr << "Unknown params_globalSfM input format: " << ext << std::endl;
    return false;
  }

  return bStatus;
}

bool Save(const paramsGlobalSfM & params_data, const std::string & filename)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(params_data, filename);
  else
  {
    std::cerr << "Unknown params_globalSfM export format: " << ext << std::endl;
  }
  return false;
}

} // namespace params
} // namespace openMVG

