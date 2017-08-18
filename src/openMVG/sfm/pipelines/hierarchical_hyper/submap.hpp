// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/sfm_data.hpp"

#include <set>

namespace openMVG {
namespace sfm {

// struct used to store the submaps created by clustering.
// A HsfmSubmap can be a parent, if it has been clustered into two
// children submaps.
// It stores a sfm_data, and the tracks ids corresponding to the landmarks contained in the submap.
// If it is a parent it also stores the ids of the separator tracks (which are tracks contained in
// both children submaps), and the submap id of its children.
struct HsfmSubmap
{
  SfM_Data sfm_data;
  std::set<IndexT> track_ids;
  bool is_parent = false;
  std::set<IndexT> separator;// set of track ids in the separator (if parent)
  std::pair<IndexT, IndexT> children_submaps = {0,0}; // pair of submaps index for children (if any)
  IndexT parent_id = 0;

  /** Serialization
  * @param ar Archive
  */
  template <class Archive>
  inline void serialize ( Archive & ar );
};

using HsfmSubmaps = std::map<openMVG::IndexT, openMVG::sfm::HsfmSubmap>;

// export submap data as sfm_data.json and ply (with separator tracks highlighted)
// in output directory
bool ExportSubmapData(const HsfmSubmaps & submaps, const IndexT submap_id, const std::string & file_path);

} // namespace sfm
} // namespace openMVG
