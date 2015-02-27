
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SFM_DATA_IO_HPP
#define OPENMVG_SFM_DATA_IO_HPP

#include "openMVG/stl/stlMap.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace openMVG
{

enum ESfM_Data
{
  VIEWS       = 0x01,
  EXTRINSICS  = 0x02,
  INTRINSICS  = 0x04,
  STRUCTURE   = 0x08,
  ALL = VIEWS | EXTRINSICS | INTRINSICS | STRUCTURE
};

} // namespace openMVG

#include "openMVG/sfm/sfm_data_io_cereal.hpp"

namespace openMVG {

///Check that each pose have a valid intrinsic and pose id in the existing View ids
bool ValidIds(const SfM_Data & sfm_data)
{
  std::set<IndexT> set_id_intrinsics;
  transform(sfm_data.getIntrinsics().begin(), sfm_data.getIntrinsics().end(),
    std::inserter(set_id_intrinsics, set_id_intrinsics.begin()), std::RetrieveKey());

  std::set<IndexT> set_id_extrinsics; //unique so can use a set
  transform(sfm_data.getPoses().begin(), sfm_data.getPoses().end(),
    std::inserter(set_id_extrinsics, set_id_extrinsics.begin()), std::RetrieveKey());

  // Collect existing id_intrinsic && id_extrinsic from views
  std::set<IndexT> reallyDefined_id_intrinsics;
  std::set<IndexT> reallyDefined_id_extrinsics;
  for (Views::const_iterator iter = sfm_data.getViews().begin();
    iter != sfm_data.getViews().end();
    ++iter)
  {
    // If a pose is defined, at least the intrinsic must be valid,
    // In order to generate a valid camera.
    const IndexT id_pose = iter->second.id_pose;
    const IndexT id_intrinsic = iter->second.id_intrinsic;

    if (set_id_extrinsics.count(id_pose))
      reallyDefined_id_extrinsics.insert(id_pose); //at least it exists
    if (set_id_extrinsics.count(id_intrinsic))
      reallyDefined_id_intrinsics.insert(id_intrinsic); //at least it exists
  }
  // Check if all defined intrinsic & extrinsic are at least connected to a view
  return
    (
    set_id_intrinsics.size() == reallyDefined_id_intrinsics.size() &&
    set_id_extrinsics.size() == reallyDefined_id_extrinsics.size()
    );
}

bool Load(SfM_Data & sfm_data, const std::string & filename, ESfM_Data flags_part)
{
  bool bStatus = false;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    bStatus = Load_Cereal<cereal::JSONInputArchive>(sfm_data, filename, flags_part);
  else if (ext == "bin")
    bStatus = Load_Cereal<cereal::PortableBinaryInputArchive>(sfm_data, filename, flags_part);
  else if (ext == "xml")
    bStatus = Load_Cereal<cereal::XMLInputArchive>(sfm_data, filename, flags_part);
  else return false;

  // Assert that loaded intrinsics | extrinsics are linked to valid view
  if (
    (flags_part & VIEWS) == VIEWS && (
    (flags_part & INTRINSICS) == INTRINSICS ||
    (flags_part & EXTRINSICS) == EXTRINSICS))
  {
    return ValidIds(sfm_data);
  }
  return bStatus;
}

bool Save(const SfM_Data & sfm_data, const std::string & filename, ESfM_Data flags_part)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "json")
    return Save_Cereal<cereal::JSONOutputArchive>(sfm_data, filename, flags_part);
  else if (ext == "bin")
    return Save_Cereal<cereal::PortableBinaryOutputArchive>(sfm_data, filename, flags_part);
  else if (ext == "xml")
    return Save_Cereal<cereal::XMLOutputArchive>(sfm_data, filename, flags_part);
  return false;
}

} // namespace openMVG

#endif // OPENMVG_SFM_DATA_IO_HPP
