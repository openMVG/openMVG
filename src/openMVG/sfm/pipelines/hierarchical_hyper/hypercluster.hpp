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

// Class for clustering sfm data using hypergraph partitioning
// (HyperSfM method)
class HyperCluster
{

public:
  HyperCluster(const sfm::SfM_Data & sfm_data, const tracks::STLMAPTracks & map_tracks, const int threshold_submap_tracksize);

  // create clustering tree by recursively partitioning submaps into smaller ones
  bool recursivePartitioning();

  // visualization using graphviz
  bool exportTreeGraph(const std::string & filename) const;

  // accessors
  HsfmSubmaps getSubMaps() const {return submaps_;}

private:

  // Partition a single submap into two children submaps, using a hypergraph
  // partitioning method.
  bool PartitionSubmap(const IndexT submap_id, std::vector<sfm::HsfmSubmap> & partitioned_pair);

  std::map<std::set<IndexT>, std::set<size_t>> createSubHyperGraph(IndexT submap_id);

  // input sfm data
  sfm::SfM_Data root_sfm_data_;

  // Tracks container
  tracks::STLMAPTracks map_tracks_;

  // container for submaps
  HsfmSubmaps submaps_;

  int threshold_submap_tracksize_;
};

} // namespace sfm
} // namespace openMVG
