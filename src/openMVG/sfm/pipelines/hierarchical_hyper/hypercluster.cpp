// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include "ceres/rotation.h"
#include "openMVG/sfm/pipelines/hierarchical_hyper/hypercluster.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/hypergraph_partitioning_scotch.hpp"

namespace openMVG {
namespace sfm {

SfM_Data create_sub_sfm_data(const SfM_Data & parent_sfm_data, const std::set<IndexT> & view_ids)
{
  SfM_Data out_sfm_data = parent_sfm_data;
  for (sfm::Views::iterator view_it = out_sfm_data.views.begin();
      view_it != out_sfm_data.views.end();)
  {
    if (view_ids.find(view_it->second->id_view) == view_ids.end())
    {
      view_it = out_sfm_data.views.erase(view_it);
      continue;
    }
    ++view_it;
  }

  // remove unused intrinsics
  for (sfm::Intrinsics::iterator intrinsic_it = out_sfm_data.intrinsics.begin();
      intrinsic_it != out_sfm_data.intrinsics.end();)
  {
    const IndexT & id_intrinsic = intrinsic_it->first;

    bool found_intrinsic = false;

    for (const auto & view : out_sfm_data.views)
    {
      if (view.second->id_intrinsic == id_intrinsic)
      {
        found_intrinsic = true;
        break;
      }
    }
    if (!found_intrinsic)
      intrinsic_it = out_sfm_data.intrinsics.erase(intrinsic_it);
    else
      ++intrinsic_it;
  }
  return out_sfm_data;
}

HyperCluster::HyperCluster(const sfm::SfM_Data & sfm_data, const tracks::STLMAPTracks & map_tracks, const int threshold_submap_tracksize)
  : root_sfm_data_(sfm_data), map_tracks_(map_tracks), threshold_submap_tracksize_(threshold_submap_tracksize), MIN_NUMBER_VIEWS_PER_SUBMAPS_(2)
{
  root_sfm_data_.structure.clear();
  root_sfm_data_.poses.clear();
}

bool HyperCluster::submapsAreLeftToBePartitioned(const std::set<openMVG::IndexT> & non_partitionable_smap_ids)
{
  // a partitionable submap should have the following requirements :
  // - not a parent
  // - not flagged at "non partitionable" (see input param)
  // - be above the clustering threshold

  return
    std::any_of(submaps_.begin(), submaps_.end(),
       [=](std::pair<IndexT, HsfmSubmap> a)
       {return ((!a.second.is_parent)
                && non_partitionable_smap_ids.count(a.first) == 0
                && (a.second.track_ids.size() > threshold_submap_tracksize_));});
}

bool HyperCluster::recursivePartitioning()
{
  // initialize stuff
  submaps_.clear();

  for (const auto & track : map_tracks_)
  {
    submaps_[0].track_ids.insert(track.first);
  }
  submaps_[0].sfm_data = root_sfm_data_;

  std::set<openMVG::IndexT> non_partitionable_submap_ids;

  while (submapsAreLeftToBePartitioned(non_partitionable_submap_ids))
  {
    std::cout << "some submaps are still too big ! >> keep clustering !" << std::endl;

    HsfmSubmaps new_submaps;

    std::cout << " submaps size : " << std::endl;
    for (const auto & smap : submaps_)
      std::cout << smap.first << " : " << smap.second.track_ids.size()
        << " tracks, " << smap.second.sfm_data.GetViews().size() << " views, is parent : " << smap.second.is_parent << std::endl;

    // we will append the child submaps to the end of the submaps map.
    int next_index = submaps_.rbegin()->first + 1;

    // partition the submaps which are too big
    for (auto & smap : submaps_)
    {
      HsfmSubmap & submap = smap.second;
      const openMVG::IndexT & submap_id = smap.first;

      if (!submap.is_parent && submap.track_ids.size() > threshold_submap_tracksize_)
      {
        // partition the current submap
        std::vector<HsfmSubmap> new_submap_pair(2);
        if(!PartitionSubmap(submap_id, new_submap_pair))
        {
          std::cout << "stop partitioning for submap nb " << submap_id << std::endl;
          non_partitionable_submap_ids.insert(submap_id);
          continue;
        }

        // submap becomes parent
        submap.is_parent = true;

        // determine id for children submaps
        submap.children_submaps = {next_index, next_index + 1};

        // add the two newly created submaps to the new submaps list
        new_submaps[next_index] = new_submap_pair[0];
        new_submaps[next_index + 1] = new_submap_pair[1];
        next_index += 2;
      }
    }
    std::cout << "number of submaps to be added : " << new_submaps.size() << std::endl;

    // add new submaps to existing submaps
    submaps_.insert(new_submaps.begin(), new_submaps.end());
  }
  return true;
}

bool HyperCluster::exportTreeGraph(const std::string & filename) const
{
  // Export the graph as a DOT (graph description language) file
  std::ofstream file(filename.c_str());
  file << "graph hypercluster{" << std::endl;
  file << "ranksep=.1;" << std::endl;
  file << "size = \"7.5,7.5\";" << std::endl;
  file << "node [shape=circle]" << std::endl;

  for (const auto & smap : submaps_)
  {
    if (smap.second.children_submaps.first != 0)
    {
      file << smap.first << " -- " << smap.second.children_submaps.first << "\n";
      file << smap.first << " -- " << smap.second.children_submaps.second << "\n";
    }
    else
    {
      // link to view and track number!
      file << smap.first << " -- v_" << smap.second.sfm_data.GetViews().size()
        << "_t_" << smap.second.track_ids.size() << "\n";
    }
  }
  file << "}" << std::endl;

  file.close();
  std::cout << "dot file written : convert to svg!" << std::endl;

  // convert to svg with dot
  const std::string dotcommand = "dot -Tsvg -O -Goverlap=scale -Gsplines=false -v " + filename;
  if (std::system(nullptr) != 0){
    std::system(dotcommand.c_str());
  }
  std::cout << "done!" << std::endl;
  return true;
}

bool HyperCluster::PartitionSubmap(const IndexT submap_id, std::vector<sfm::HsfmSubmap> & partitioned_pair)
{
  partitioned_pair.clear();

  sfm::SfM_Data & parent_sfm_data = submaps_.at(submap_id).sfm_data;

  // compute hyper graph corresponding to the tracks of the submap
  std::map<std::set<IndexT>, std::set<size_t>> sub_hyper_edges_and_tracks = createSubHyperGraph(submap_id);

  // use external library to partition the submap
  std::pair<std::set<IndexT>, std::set<IndexT>> view_id_partitions;

  if (!ScotchPartitionHyperGraph(sub_hyper_edges_and_tracks, view_id_partitions))
  {
    std::cerr << "SCOTCH could not partition hyper graph for submap " << submap_id << std::endl;
    return false;
  }

  if (view_id_partitions.first.size() < MIN_NUMBER_VIEWS_PER_SUBMAPS_)
  {
    std::cerr << "returned partition is too small ! (only " << view_id_partitions.first.size() << " views." << std::endl;
    return false;
  }

  // create the children submaps and fill the view ids
  // with the partitioned views
  HsfmSubmap first_submap, second_submap;
  first_submap.sfm_data = create_sub_sfm_data(parent_sfm_data, view_id_partitions.first);
  second_submap.sfm_data = create_sub_sfm_data(parent_sfm_data, view_id_partitions.second);
  first_submap.parent_id = submap_id;
  second_submap.parent_id = submap_id;
  partitioned_pair.push_back(first_submap);
  partitioned_pair.push_back(second_submap);

  // deduce separator (set of hyperedges that are contained in both submaps.)
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (const auto & hyper_edge : sub_hyper_edges_and_tracks)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      // find if hyper_edge is a subset of one of the children submaps
      // (if yes, it can not be a separator by definition since it is fully contained
      // in a single submap).
      bool is_separator = true;
      for (auto & smap : partitioned_pair)
      {
        std::set<IndexT> view_ids;
        for (const auto & view : smap.sfm_data.GetViews())
          view_ids.insert(view.first);

        if (std::includes(view_ids.begin(), view_ids.end(),
              hyper_edge.first.begin(), hyper_edge.first.end()))
        {
          is_separator = false;
#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
          smap.track_ids.insert(hyper_edge.second.begin(), hyper_edge.second.end());
          break;
        }
      }

#ifdef OPENMVG_USE_OPENMP
  #pragma omp critical
#endif
      if (is_separator)
      {
        // if it reaches this point -> it means this edge is part of the separator
        submaps_.at(submap_id).separator.insert(hyper_edge.second.begin(), hyper_edge.second.end());

        // also add the tracks ids to both children submaps (in the separator == in both submaps)
        for (auto & smap : partitioned_pair)
        {
          smap.track_ids.insert(hyper_edge.second.begin(), hyper_edge.second.end());
        }
      }
    }// omp section
  }

  return true;
}

std::map<std::set<IndexT>, std::set<size_t>> HyperCluster::createSubHyperGraph(IndexT submap_id)
{
  std::map<std::set<IndexT>, std::set<size_t>> sub_hyper_edges_and_tracks;
  if (map_tracks_.size() == 0)
    std::cerr << "no tracks ! cannot compute sub hypergraph." <<std::endl;

  // list of the view ids actually available in the current submap
  std::set<IndexT> submap_viewids;
  for (const auto & view : submaps_.at(submap_id).sfm_data.GetViews())
    submap_viewids.insert(view.first);

  // count tracks by view_ids sequences, sort result.
  // ex : the number of tracks that cover views 1,2 and 4 will
  // be the weight of the hyper edge {1,2,4}
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (std::set<IndexT>::const_iterator track_id = submaps_.at(submap_id).track_ids.begin();
      track_id != submaps_.at(submap_id).track_ids.end(); ++track_id)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      const auto & track = map_tracks_.at(*track_id);

      std::set<IndexT> view_ids;

      // list the views related to current track
      for (const auto & track_part : track)
      {
        // we don't want to add view ids that are not in the submap
        // (basically this condition is for dealing with separator tracks) TODO < do we want to deal with it this way ??
        if (submap_viewids.find(track_part.first) != submap_viewids.end())
          view_ids.insert(track_part.first);
      }

#ifdef OPENMVG_USE_OPENMP
    #pragma omp critical
#endif
      sub_hyper_edges_and_tracks[view_ids].insert(*track_id);
    }// omp section
  }

  return sub_hyper_edges_and_tracks;
}

}
}
