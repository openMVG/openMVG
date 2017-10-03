// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * This file contains interfaces to the SCOTCH hypergraph partitioning library
 * that we use for hypersfm.
 */

#pragma once

#include "openMVG/types.hpp"

/**
 * @note WARNING : this function USED TO create a memory leak (the incriminated part is SCOTCH_meshGraph() in
 * the scotch library)
 * you have to use our version of the scotch library which corrects this memory leak (in third_party/scotch_6.0.4).
 * @brief Use the scotch library to partition a hypergraph made of view_id nodes and tracks set hyperedges
 * into two sets of nodes (view ids in that case)
 * @param hyper_graph The hypergraph to be partitioned, defined as a map of <view_ids,tracks_ids>.
 * @param view_id_partitions The returned partition, in the shape of a pair of view_ids sets.
 * @retval true If partitioning succeeded
 * @retval false If partitioning failed
 */
bool ScotchPartitionHyperGraph(
    const std::map<std::set<openMVG::IndexT>, std::set<size_t>> & hyper_graph,
    std::pair<std::set<openMVG::IndexT>, std::set<openMVG::IndexT>> & view_id_partitions);
