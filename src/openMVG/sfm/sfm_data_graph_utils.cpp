// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015,2017 WhuAegeanSea, Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_graph_utils.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_builder.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/types.hpp"
#include "openMVG/sfm/sfm_filters.hpp"

#include <fstream>
#include <iosfwd>

namespace openMVG {
namespace sfm {

void KeepOnlyReferencedElement(
	const std::set<IndexT> & set_remainingIds,
	const matching::PairWiseMatches& map_matches_in,
	matching::PairWiseMatches& map_matches_out)
{
	matching::PairWiseMatches map_matches_E_infered;
	for (matching::PairWiseMatches::const_iterator iter = map_matches_in.begin();
		iter != map_matches_in.end(); ++iter)
	{
		if (set_remainingIds.count(iter->first.first) &&
			set_remainingIds.count(iter->first.second))
		{
			map_matches_E_infered.insert(*iter);
		}
	}
	map_matches_out.swap(map_matches_E_infered);
}

bool SplitMatchFileIntoMatchFiles(SfM_Data & sfm_data, const std::string & match_file,
	const std::string & match_file_components, bool bBiEdge, int nMinNode)
{
	if (!stlplus::file_exists(match_file))
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

	if (bBiEdge)
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


	std::vector<std::set<IndexT>> vec_subgraph;
	const int connectedComponentCount = lemon::countConnectedComponents(putativeGraph.g);
	if (connectedComponentCount >= 1)
	{
		const std::map<IndexT, std::set<Graph::Node> > map_subgraphs = graph::exportGraphToMapSubgraphs<Graph, IndexT>(putativeGraph.g);
		for (typename std::map<IndexT, std::set<Graph::Node> >::const_iterator iter = map_subgraphs.begin();
			iter != map_subgraphs.end(); ++iter)
		{
			if (iter->second.size() > nMinNode)
			{
				std::set<IndexT> subgraphNodes;
				const std::set<lemon::ListGraph::Node> & ccSet = iter->second;
				for (const auto & iter2 : ccSet)
				{
					const IndexT Id = (*putativeGraph.node_map_id)[iter2];
					subgraphNodes.insert(Id);
				}
				vec_subgraph.emplace_back(subgraphNodes);
			}
		}
	}

	std::set<std::string> setFileName;

	const std::string &sBaseName = stlplus::basename_part(match_file);
	const std::string &sOutputFolder = stlplus::folder_part(match_file_components);
	const std::string &sExt = stlplus::extension_part(match_file);
	int nIndex = 0;
	for (auto subgraph : vec_subgraph)
	{
		std::stringstream strstream_FileName;
		strstream_FileName << sBaseName << "_" << nIndex << "_" << subgraph.size() << "." << sExt;
		const std::string &sOutputFileName = strstream_FileName.str();
		const std::string &sMatch_file_components = stlplus::create_filespec(sOutputFolder, sOutputFileName);

		matching::PairWiseMatches map_matches_out;
		KeepOnlyReferencedElement(subgraph, matches_provider->pairWise_matches_, map_matches_out);
		matching::Save(map_matches_out, sMatch_file_components);

		setFileName.insert(sOutputFileName);
		++nIndex;
	}

	std::ofstream stream(match_file_components.c_str());
	if (!stream.is_open())
	{
		return false;
	}
	for (auto sFileName : setFileName)
	{
		stream << sFileName << std::endl;
	}
	stream.close();
	return true;
}

} // namespace sfm
} // namespace openMVG
