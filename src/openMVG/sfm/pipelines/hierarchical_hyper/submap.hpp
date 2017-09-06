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

/**
 * @brief The HsfmSubmap struct, a data structure for clustered data.
 * A HsfmSubmap contains a sfm_data representing the subscene,
 * a set of track_ids corresponding to the landmarks contained in the submap.
 * It contains a flag to signify if it is a parent (i.e. if it has itself been clustered
 * into two children submaps).
 * for a parent submap :
 * - It stores the ids of the separator tracks (the tracks contained in both
 *   children submaps).
 * - It contains the submap indices of its children
 * It also stores the index of its own parent.
 */
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
