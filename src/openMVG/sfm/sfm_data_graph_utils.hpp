// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015, 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
#define OPENMVG_SFM_DATA_GRAPH_UTILS_HPP

#include "openMVG/types.hpp"
#include "openMVG/matching/indMatch.hpp"
#include <string>

namespace openMVG {
namespace sfm {
///  @brief Split match pairs into connected match pairs
///  @param[in]  pairs   The pair sets of the images.
///  @param[in]  is_biedge  Set true for global SFM, false for sequential SFM.
///  @param[in]  min_nodes  The minimum nodes of graph in the output match_file. Note: the value should be larger than 3.
///  @param[out] subgraps_ids The pairs id of subgraphs split from the whole graph.
///
///  @return True if some components have been kept
bool PairsToConnectedComponents
(
  const Pair_Set & pairs,
  bool is_biedge,
  int min_nodes,
  std::map<IndexT, std::set<IndexT>> & subgraphs_ids
);

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
///  @param[in]  pairs   The pair sets of the images.
///  @param[in]  matches The pairwise matches of the images corresponding to pairs.
///  @param[in]  is_biedge  Set true for global SFM, false for sequential SFM.
///  @param[in]  min_nodes  The minimum nodes of graph in the output match_file. Note: the value should be larger than 3.
///  @param[out] subgraphs_matches The matches of subgraphs split from the whole graph.
///
///  @return True if subgraphs_matches is not empty
///
bool SplitMatchesIntoSubgraphMatches
(
  const Pair_Set & pairs,
  const matching::PairWiseMatches & matches,
  bool is_biedge,
  int min_nodes,
  std::vector<matching::PairWiseMatches> & subgraphs_matches
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_GRAPH_UTILS_HPP
