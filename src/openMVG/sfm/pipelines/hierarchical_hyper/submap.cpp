// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/sfm_data_BA_fixed_points.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_utilities.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

namespace openMVG {
namespace sfm {

bool ExportSubmapData(const HsfmSubmaps & submaps, const IndexT submap_id, const std::string & file_path)
{
  const HsfmSubmap & smap = submaps.at(submap_id);
  // save sfm data
  Save(smap.sfm_data,
      stlplus::create_filespec(stlplus::folder_part(file_path), stlplus::basename_part(file_path), "json"),
      ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS));
  // save sfm data bin file
  Save(smap.sfm_data,
      stlplus::create_filespec(stlplus::folder_part(file_path), stlplus::basename_part(file_path), "bin"),
      ESfM_Data(ALL));

  std::string filename = stlplus::create_filespec(stlplus::folder_part(file_path), stlplus::basename_part(file_path), "ply");

  const IndexT parent_submap_id = submaps.at(submap_id).parent_id;

  // find separator tracks
  const std::set<IndexT> & separator_tracks_indices = submaps.at(parent_submap_id).separator;

  std::set<IndexT> separator_points_indices;
  std::set<IndexT> separator_camera_indices = {};// TODO < fill this ?

  // find 3d points and camera positions
  std::vector<Vec3> points, cam_pos;
  unsigned i(0);
  std::map<IndexT, unsigned> map_track_to_point_id;
  for (const auto & landmark : smap.sfm_data.GetLandmarks())
  {
    points.push_back(landmark.second.X);
    if (separator_tracks_indices.find(landmark.first) != separator_tracks_indices.end())
    {
      separator_points_indices.insert(i);
      map_track_to_point_id[landmark.first] = i;
    }
    i++;
  }

  // find coordinates of separator points also reconstructed in sibling submap
  std::set<IndexT> common_separator_points_indices;
  const IndexT sibling_submap_id = getSiblingSubmapId(submaps, submap_id);
  for (const auto & landmark : submaps.at(sibling_submap_id).sfm_data.GetLandmarks())
  {
    const IndexT & track_id = landmark.first;
    if (separator_tracks_indices.count(track_id) != 0)
      common_separator_points_indices.insert(map_track_to_point_id[track_id]);
  }

  for (const auto & view : smap.sfm_data.GetViews())
  {
    if (smap.sfm_data.IsPoseAndIntrinsicDefined(view.second.get()))
      cam_pos.push_back(smap.sfm_data.GetPoseOrDie(view.second.get()).center());
  }

  std::set<IndexT> non_common_separator_points_indices;
  std::set_difference(separator_points_indices.begin(), separator_points_indices.end(),
      common_separator_points_indices.begin(), common_separator_points_indices.end(),
      std::inserter(non_common_separator_points_indices, non_common_separator_points_indices.begin()));

  if (!plyHelper::exportToPly_MultiHighlight(
      points, cam_pos, filename, {non_common_separator_points_indices, common_separator_points_indices}, separator_camera_indices))
  {
    return false;
  }

  return true;
}

}
}
