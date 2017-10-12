// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015, 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
#define OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
#include <string>

namespace openMVG {
namespace sfm {

struct SfM_Data;

///  @brief Split match_file into match_files
///  When handling a scene which may contains many disconnected graphs, 
///  for the moment OpenMVG consider only the largest connected component for SfM.
///  Especially, when the GPS data of the images known, the disconnected graphs in a scene should be 
///  merged since they all share the same coordinate system.
///  Steps handling the above scenario.
///  Step 1 : Use this function to compute the various connected component and export 
///  the matches of the each connected component in a different match file.
///  Step 2 : Run the SfM pipeline of your choice on the produced matches file.
///  Step 3 : Merge all the sfm_data into a single sfm_data and 
///  triangulate the initial match file, when the GPS data of the images known.
///         
///
///  @param[in]  sfm_data   The sfm_data.
///  @param[in]  match_file The match file which contains all nodes of the whole graph.
///  @param[in]  match_component_filename The output match_component_file contains a list of match files.
///  @param[in]  is_biedge  Set true for global SFM, false for sequential SFM.
///  @param[out] min_nodes  The minimum nodes of graph in the output match_file.Note: the value should be larger than 3.
///
///  @return True if the area can be computed
///
bool SplitMatchFileIntoMatchFiles(const SfM_Data & sfm_data, const std::string & match_file,
  const std::string & match_component_filename, bool is_biedge, int min_nodes);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
