// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cereal/types/set.hpp>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io_cereal.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/submap.hpp"
#include "openMVG/tracks/tracks.hpp"

/*
 * Inspired by paper :
 * "HyperSFM", Kai Ni and Frank Dellaert,
 * International Conference on 3D Imaging, Modeling, Processing, Visualization and Transmission (3DIMPVT), 2012
 * link : http://kaini.org/assets/Ni12_3dimpvt.pdf HyperSFM (2012)
 *
 * Robert Maier's Master thesis was also of great help in order
 * to help implement hyperSFM :
 * "Out-of-Core Bundle Adjustment for 3D Workpiece Reconstruction"
 * in particular , paragraphs 2.2.4 and 4.5
 * link : https://vision.in.tum.de/_media/spezial/bib/maier2013thesis.pdf
 */

namespace openMVG {
namespace sfm {


/**
 * @brief The HyperCluster class, can cluster a sfm_data into multiple submaps.
 */
class HyperCluster
{

public:
  HyperCluster(const sfm::SfM_Data & sfm_data, const tracks::STLMAPTracks & map_tracks, const int threshold_submap_tracksize);

  /**
   * @brief cluster the scene into a binary tree by recursively partitioning submaps
   * @return [bool] success flag
   */
  bool recursivePartitioning();

  /**
   * @brief exportTreeGraph visualization of the clustering (uses graphviz)
   * @param filename where the graph is to be saved
   * @return bool success flag
   */
  bool exportTreeGraph(const std::string & filename) const;

  /**
   * @brief accessor for submaps
   * @return all submaps
   */
  HsfmSubmaps getSubMaps() const {return submaps_;}

private:

  /**
   * @brief Partition a single submap into two children submaps, using the hypergraph
   * partitioning method.
   * @param [in] submap_id
   * @param [out] the returned partitioned pair,
   * @return flag that expresses if the partition is usable or not
   * @note a partitioned pair is considered a failure if the partition contains less than the
   * MIN_NUMBER_VIEWS_PER_SUBMAPS_
   */
  bool PartitionSubmap(const IndexT submap_id, std::vector<sfm::HsfmSubmap> & partitioned_pair);

  /**
   * @brief creates the hypergraph corresponding to a submap (not necessarily the root submap)
   * @param the id of the submap
   * @return the hypergraph as a map of hyperedges(set of view ids) to the corresponding track ids (set of track ids)
   */
  std::map<std::set<IndexT>, std::set<size_t>> createSubHyperGraph(IndexT submap_id);

  bool submapsAreLeftToBePartitioned(const std::set<openMVG::IndexT> & non_partitionable_smap_ids = {});

  /**
   * @brief the input (root) sfm_data
   */
  sfm::SfM_Data root_sfm_data_;

  /**
   * @brief container for the tracks in the scene
   */
  tracks::STLMAPTracks map_tracks_;

  /**
   * @brief container for the submaps
   */
  HsfmSubmaps submaps_;

  /**
   * @brief the maximum number of tracks we want in a submap
   */
  int threshold_submap_tracksize_;

  /**
   * @brief the minimum amount of views we want in a submap
   * @note this should be AT LEAST 2 if you want to be able to do
   * any reconstruction with the submaps later on
   */
  int MIN_NUMBER_VIEWS_PER_SUBMAPS_;
};

} // namespace sfm
} // namespace openMVG
