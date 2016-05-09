
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Implementation of [1] an efficient algorithm to compute track from pairwise
//  correspondences.
//
//  [1] Pierre Moulon and Pascal Monasse,
//    "Unordered feature tracking made fast and easy" CVMP 2012.
//
// It tracks the position of features along the series of image from pairwise
//  correspondences.
//
// From map< [imageI,ImageJ], [indexed matches array] > it builds tracks.
//
// Usage :
//  PairWiseMatches map_Matches;
//  PairedIndMatchImport(sMatchFile, map_Matches); // Load series of pairwise matches
//  //---------------------------------------
//  // Compute tracks from matches
//  //---------------------------------------
//  TracksBuilder tracksBuilder;
//  tracks::STLMAPTracks map_tracks;
//  tracksBuilder.Build(map_Matches); // Build: Efficient fusion of correspondences
//  tracksBuilder.Filter();           // Filter: Remove tracks that have conflict
//  tracksBuilder.ExportToSTL(map_tracks); // Build tracks with STL compliant type
//

#ifndef OPENMVG_TRACKS_H_
#define OPENMVG_TRACKS_H_

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/tracks/union_find.hpp"
#include "openMVG/tracks/flat_pair_map.hpp"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace openMVG  {

using namespace openMVG::matching;

namespace tracks  {

// Data structure to store a track: collection of {ImageId,FeatureId}
//  The corresponding image points with their imageId and FeatureId.
typedef std::map<uint32_t,uint32_t> submapTrack;
// A track is a collection of {trackId, submapTrack}
typedef std::map< size_t, submapTrack > STLMAPTracks;

struct TracksBuilder
{
  typedef std::pair<uint32_t, uint32_t>indexedFeaturePair;

  flat_pair_map<indexedFeaturePair, uint32_t> map_node_to_index;
  UnionFind uf_tree;

  /// Build tracks for a given series of pairWise matches
  void Build( const PairWiseMatches &  map_pair_wise_matches)
  {
    // 1. We need to know how much single set we will have.
    //   i.e each set is made of a tuple : (imageIndex, featureIndex)
    typedef std::set<indexedFeaturePair> SetIndexedPair;
    SetIndexedPair allFeatures;
    // For each couple of images list the used features
    for ( const auto & iter : map_pair_wise_matches )
    {
      const auto & I = iter.first.first;
      const auto & J = iter.first.second;
      const std::vector<IndMatch> & vec_FilteredMatches = iter.second;

      // Retrieve all shared features and add them to a set
      for( const auto & cur_filtered_match : vec_FilteredMatches )
      {
        allFeatures.emplace(I,cur_filtered_match.i_);
        allFeatures.emplace(J,cur_filtered_match.j_);
      }
    }

    // 2. Build the 'flat' representation where a tuple (the node)
    //  is attached to an unique index.
    map_node_to_index.reserve(allFeatures.size());
    unsigned int cpt = 0;
    for (const auto & feat : allFeatures)
    {
      map_node_to_index.emplace_back(feat, cpt);
      ++cpt;
    }
    // Sort the flat_pair_map
    map_node_to_index.sort();
    // Clean some memory
    allFeatures.clear();

    // 3. Add the node and the pairwise correpondences in the UF tree.
    uf_tree.InitSets(map_node_to_index.size());

    // 4. Union of the matched features corresponding UF tree sets
    for ( const auto & iter : map_pair_wise_matches )
    {
      const auto & I = iter.first.first;
      const auto & J = iter.first.second;
      const std::vector<IndMatch> & vec_FilteredMatches = iter.second;
      for (const IndMatch & match : vec_FilteredMatches)
      {
        const indexedFeaturePair pairI(I, match.i_);
        const indexedFeaturePair pairJ(J, match.j_);
        // Link feature correspondences to the corresponding containing sets.
        uf_tree.Union(map_node_to_index[pairI], map_node_to_index[pairJ]);
      }
    }
  }

  /// Remove bad tracks (too short or track with ids collision)
  bool Filter(size_t nLengthSupTo = 2)
  {
    // Remove bad tracks:
    // - track that are too short,
    // - track with id conflicts:
    //    i.e. tracks that have many times the same image index

    // From the UF tree, create tracks of the image indexes.
    //  If an image index appears two time the track must disappear
    //  If a track is too short it has to be removed.
    std::map<unsigned int, std::set<unsigned int> > tracks;

    std::set<unsigned int> problematic_track_id;
    // Build tracks from the UF tree, track problematic ids.
    for (unsigned int k = 0; k < map_node_to_index.size(); ++k)
    {
      const unsigned int & track_id = uf_tree.m_cc_parent[k];
      if (problematic_track_id.count(track_id) != 0)
        continue; // Track already marked

      const auto & feat = map_node_to_index[k];

      if (tracks[track_id].count(feat.first.first))
      {
        problematic_track_id.insert(track_id);
      }
      else
      {
        tracks[track_id].insert(feat.first.first);
      }
    }

    // - track that are too short,
    for (const auto & val : tracks)
    {
      if (val.second.size() < nLengthSupTo)
      {
        problematic_track_id.insert(val.first);
      }
    }

    for (unsigned int & root_index : uf_tree.m_cc_parent)
    {
      if (problematic_track_id.count(root_index) > 0)
      {
        // reset selected root
        uf_tree.m_cc_size[root_index] = 1;
        root_index = std::numeric_limits<unsigned int>::max();
      }
    }
    return false;
  }

  /// Return the number of connected set in the UnionFind structure (tree forest)
  size_t NbTracks() const
  {
    std::set<unsigned int> parent_id(uf_tree.m_cc_parent.begin(), uf_tree.m_cc_parent.end());
    // Erase the "special marker" that depicted rejected tracks
    parent_id.erase(std::numeric_limits<unsigned int>::max());
    return parent_id.size();
  }

  /// Export tracks as a map (each entry is a sequence of imageId and featureIndex):
  ///  {TrackIndex => {(imageIndex, featureIndex), ... ,(imageIndex, featureIndex)}
  void ExportToSTL(STLMAPTracks & map_tracks)
  {
    map_tracks.clear();
    for (unsigned int k = 0; k < map_node_to_index.size(); ++k)
    {
      const auto & feat = map_node_to_index[k];
      const unsigned int track_id = uf_tree.m_cc_parent[k];
      if
      (
        // ensure never add rejected elements (track marked as invalid)
        track_id != std::numeric_limits<unsigned int>::max()
        // ensure never add 1-length track element (it's not a track)
        && uf_tree.m_cc_size[track_id] > 1
      )
      {
        map_tracks[track_id].insert(feat.first);
      }
    }
  }
};

struct TracksUtilsMap
{
  /**
   * @brief Find common tracks between images.
   *
   * @param[in] set_imageIndex: set of images we are looking for common tracks
   * @param[in] map_tracksIn: all tracks of the scene
   * @param[out] map_tracksOut: output with only the common tracks
   */
  static bool GetTracksInImages
  (
    const std::set<size_t> & set_imageIndex,
    const STLMAPTracks & map_tracksIn,
    STLMAPTracks & map_tracksOut
  )
  {
    map_tracksOut.clear();

    // Go along the tracks
    for ( const auto & iterT : map_tracksIn )
    {
      // Look if the track contains the provided view index & save the point ids
      submapTrack map_temp;
      bool bTest = true;
      for (auto iterIndex = set_imageIndex.begin();
        iterIndex != set_imageIndex.end() && bTest; ++iterIndex)
      {
        auto iterSearch = iterT.second.find(*iterIndex);
        if (iterSearch != iterT.second.end())
          map_temp[iterSearch->first] = iterSearch->second;
        else
          bTest = false;
      }

      if (!map_temp.empty() && map_temp.size() == set_imageIndex.size())
        map_tracksOut[iterT.first] = std::move(map_temp);
    }
    return !map_tracksOut.empty();
  }

  /// Return the tracksId as a set (sorted increasing)
  static void GetTracksIdVector
  (
    const STLMAPTracks & map_tracks,
    std::set<size_t> * set_tracksIds
  )
  {
    set_tracksIds->clear();
    for ( const auto & iterT : map_tracks )
    {
      set_tracksIds->insert(iterT.first);
    }
  }

  /// Get feature index PerView and TrackId
  static bool GetFeatIndexPerViewAndTrackId
  (
    const STLMAPTracks & map_tracks,
    const std::set<size_t> & set_trackId,
    size_t nImageIndex,
    std::vector<size_t> * pvec_featIndex
  )
  {
    for (const size_t & trackId: set_trackId)
    {
      STLMAPTracks::const_iterator iterT = map_tracks.find(trackId);
      if (iterT != map_tracks.end())
      {
        //try to find imageIndex
        const submapTrack & map_ref = iterT->second;
        submapTrack::const_iterator iterSearch = map_ref.find(nImageIndex);
        if (iterSearch != map_ref.end())
        {
          pvec_featIndex->emplace_back(iterSearch->second);
        }
      }
    }
    return !pvec_featIndex->empty();
  }

  struct FunctorMapFirstEqual : public std::unary_function <STLMAPTracks , bool>
  {
    size_t id;
    FunctorMapFirstEqual(size_t val):id(val){};
    bool operator()(const std::pair<size_t, submapTrack > & val) {
      return (id == val.first);
    }
  };

  /**
   * @brief Convert a trackId to a vector of indexed Matches.
   *
   * @param[in]  map_tracks: set of tracks with only 2 elements
   *             (image A and image B) in each submapTrack.
   * @param[in]  vec_filterIndex: the track indexes to retrieve.
   *             Only track indexes contained in this filter vector are kept.
   * @param[out] pvec_index: list of matches
   *             (feature index in image A, feature index in image B).
   *
   * @warning The input tracks must be composed of only two images index.
   * @warning Image index are considered sorted (increasing order).
   */
  static void TracksToIndexedMatches
  (
    const STLMAPTracks & map_tracks,
    const std::vector<IndexT> & vec_filterIndex,
    std::vector<IndMatch> * pvec_index
  )
  {

    std::vector<IndMatch> & vec_indexref = *pvec_index;
    vec_indexref.clear();
    for ( const auto & filter_index : vec_filterIndex )
    {
      // Retrieve the track information from the current index i.
      auto itF =
        find_if(map_tracks.begin(), map_tracks.end(), FunctorMapFirstEqual( filter_index ) ) ;
      // The current track.
      const submapTrack & map_ref = itF->second;

      // We have 2 elements for a track.
      assert(map_ref.size() == 2);
      const IndexT indexI = (map_ref.begin())->second;
      const IndexT indexJ = (++map_ref.begin())->second;

      vec_indexref.emplace_back(indexI, indexJ);
    }
  }

  /// Return the occurrence of tracks length.
  static void TracksLength
  (
    const STLMAPTracks & map_tracks,
    std::map<size_t, size_t> & map_Occurence_TrackLength
  )
  {
    for ( const auto & iterT : map_tracks )
    {
      const size_t trLength = iterT.second.size();
      if (map_Occurence_TrackLength.end() ==
        map_Occurence_TrackLength.find(trLength))
      {
        map_Occurence_TrackLength[trLength] = 1;
      }
      else
      {
        map_Occurence_TrackLength[trLength] += 1;
      }
    }
  }

  /// Return a set containing the image Id considered in the tracks container.
  static void ImageIdInTracks
  (
    const STLMAPTracks & map_tracks,
    std::set<size_t> & set_imagesId
  )
  {
    for ( const auto & iterT : map_tracks )
    {
      const submapTrack & map_ref = iterT.second;
      for ( const auto & iter : map_ref )
      {
        set_imagesId.insert(iter.first);
      }
    }
  }
};

} // namespace tracks
} // namespace openMVG

#endif // OPENMVG_TRACKS_H_
