// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015, 2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_graph_utils.hpp"

#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_builder.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/types.hpp"

#include <fstream>
#include <iosfwd>

namespace openMVG {
namespace sfm {

/// Filter the toFilter iterable sequence into a new sequence without change the source sequence
/// (keep only the element that share a common index with the provided Ids index list).
/// See the reference: sfm_filters.hpp
static void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  const matching::PairWiseMatches & map_matches_in,
  matching::PairWiseMatches & map_matches_out)
{
  matching::PairWiseMatches map_matches_E_infered;
  for (const auto & iter_pairwise_matches : map_matches_in)
  {
    if (set_remainingIds.count(iter_pairwise_matches.first.first) &&
      set_remainingIds.count(iter_pairwise_matches.first.second))
    {
      map_matches_E_infered.insert(iter_pairwise_matches);
    }
  }
  map_matches_out.swap(map_matches_E_infered);
}

bool SplitMatchFileIntoMatchFiles(const SfM_Data & sfm_data, const std::string & match_file,
  const std::string & match_component_filename, bool is_biedge, int min_nodes)
{
  if (!stlplus::file_exists(match_file) || match_component_filename.empty())
  {
    return false;
  }

  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if (!(matches_provider->load(sfm_data, match_file)))
  {
    return false;
  }
  using Graph = graph::indexedGraph::GraphT;
  const Pair_Set pairs = matches_provider->getPairs();
  graph::indexedGraph putativeGraph(pairs);

  // For global SFM, firstly remove the not bi-edge element 
  if (is_biedge)
  {
    using EdgeMapAlias = Graph::EdgeMap<bool>;
    EdgeMapAlias cutMap(putativeGraph.g);
    if (lemon::biEdgeConnectedCutEdges(putativeGraph.g, cutMap) > 0)
    {
      // Some edges must be removed because they don't follow the biEdge condition.
      using EdgeIterator = Graph::EdgeIt;
      EdgeIterator itEdge(putativeGraph.g);
      for (EdgeMapAlias::MapIt it(cutMap); it != lemon::INVALID; ++it, ++itEdge)
      {
        if (*it)
        {
          putativeGraph.g.erase(itEdge);  // remove the not bi-edge element
        }
      }
    }
  }

  // Compute all subgraphs in the putative graph
  std::vector<std::set<IndexT>> vec_subgraph;
  const int connectedComponentCount = lemon::countConnectedComponents(putativeGraph.g);
  if (connectedComponentCount >= 1)
  {
    const std::map<IndexT, std::set<Graph::Node> > map_subgraphs = graph::exportGraphToMapSubgraphs<Graph, IndexT>(putativeGraph.g);
    for (const auto & iter_map_subgraphs : map_subgraphs)
    {
      if (iter_map_subgraphs.second.size() > min_nodes)
      {
        std::set<IndexT> subgraphNodes;
        const std::set<lemon::ListGraph::Node> & ccSet = iter_map_subgraphs.second;
        for (const auto & iter2 : ccSet)
        {
          const IndexT Id = (*putativeGraph.node_map_id)[iter2];
          subgraphNodes.insert(Id);
        }
        vec_subgraph.emplace_back(subgraphNodes);
      }
    }
  }

  // Save all of the subgraphs into match file
  std::set<std::string> set_filenames;

  const std::string &file_basename = stlplus::basename_part(match_file);
  const std::string &output_folder = stlplus::folder_part(match_component_filename);
  const std::string &match_file_extension = stlplus::extension_part(match_file);
  int index = 0;
  for (const auto & subgraph : vec_subgraph)
  {
    std::stringstream strstream_subgraph_filename;
    strstream_subgraph_filename << file_basename << "_" << index << "_" << subgraph.size() << "." << match_file_extension;
    const std::string &subgraph_filename = strstream_subgraph_filename.str();
    const std::string &subgraph_match_file_components = stlplus::create_filespec(output_folder, subgraph_filename);

    matching::PairWiseMatches subgraph_map_matches;
    KeepOnlyReferencedElement(subgraph, matches_provider->pairWise_matches_, subgraph_map_matches);

    const bool success_flag = matching::Save(subgraph_map_matches, subgraph_match_file_components);
    if (success_flag)
    {
      set_filenames.insert(subgraph_filename);
    }
    else
    {
      return false;
    }
    ++index;
  }

  // Save the match file name of subgraph into a match component file
  std::ofstream stream(match_component_filename.c_str());
  if (!stream.is_open())
  {
    return false;
  }
  for (const auto & sFileName : set_filenames)
  {
    stream << sFileName << std::endl;
  }
  stream.close();
  return true;
}

} // namespace sfm
} // namespace openMVG
