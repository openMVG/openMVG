
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/sfm_data_io.hpp"

#include "openMVG/stl/stlMap.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/sfm/sfm_data_io_cereal.hpp"
#include "openMVG/sfm/sfm_data_io_ply.hpp"
#include "openMVG/sfm/sfm_data_io_baf.hpp"

namespace openMVG {
namespace sfm {

///Check that each pose have a valid intrinsic and pose id in the existing View ids
bool ValidIds(const SfM_Data & sfm_data, ESfM_Data flags_part)
{
  const bool bCheck_Intrinsic = (flags_part & INTRINSICS) == INTRINSICS;
  const bool bCheck_Extrinsic = (flags_part & EXTRINSICS) == EXTRINSICS;

  std::set<IndexT> set_id_intrinsics;
  transform(sfm_data.GetIntrinsics().begin(), sfm_data.GetIntrinsics().end(),
    std::inserter(set_id_intrinsics, set_id_intrinsics.begin()), stl::RetrieveKey());

  std::set<IndexT> set_id_extrinsics; //unique so can use a set
  transform(sfm_data.GetPoses().begin(), sfm_data.GetPoses().end(),
    std::inserter(set_id_extrinsics, set_id_extrinsics.begin()), stl::RetrieveKey());

  // Collect existing id_intrinsic && id_extrinsic from views
  std::set<IndexT> reallyDefined_id_intrinsics;
  std::set<IndexT> reallyDefined_id_extrinsics;
  for (Views::const_iterator iter = sfm_data.GetViews().begin();
    iter != sfm_data.GetViews().end();
    ++iter)
  {
    // If a pose is defined, at least the intrinsic must be valid,
    // In order to generate a valid camera.
    const IndexT id_pose = iter->second.get()->id_pose;
    const IndexT id_intrinsic = iter->second.get()->id_intrinsic;

    if (set_id_extrinsics.count(id_pose))
      reallyDefined_id_extrinsics.insert(id_pose); //at least it exists

    if (set_id_intrinsics.count(id_intrinsic))
      reallyDefined_id_intrinsics.insert(id_intrinsic); //at least it exists
  }
  // Check if defined intrinsic & extrinsic are at least connected to views
  bool bRet = true;
  if (bCheck_Intrinsic)
    bRet &= set_id_intrinsics.size() == reallyDefined_id_intrinsics.size();

  if (bCheck_Extrinsic)
    bRet &= set_id_extrinsics.size() == reallyDefined_id_extrinsics.size();

  if (bRet == false)
    std::cout << "There is orphan intrinsics data or poses (do not depend on any view)" << std::endl;

  return bRet;
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
  else
  {
    std::cerr << "Unknown sfm_data input format: " << ext << std::endl;
    return false;
  }

  // Assert that loaded intrinsics | extrinsics are linked to valid view
  if ( bStatus &&
    (flags_part & VIEWS) == VIEWS && (
    (flags_part & INTRINSICS) == INTRINSICS ||
    (flags_part & EXTRINSICS) == EXTRINSICS))
  {
    return ValidIds(sfm_data, flags_part);
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
  else if (ext == "ply")
    return Save_PLY(sfm_data, filename, flags_part);
  else if (ext == "baf") // Bundle Adjustment file
    return Save_BAF(sfm_data, filename, flags_part);
  else
  {
    std::cerr << "Unknown sfm_data export format: " << ext << std::endl;
  }
  return false;
}

} // namespace sfm
} // namespace openMVG


