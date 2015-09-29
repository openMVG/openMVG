
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
//  tracksBuilder.Filter();           // Filter: Remove track that have conflict
//  tracksBuilder.ExportToSTL(map_tracks); // Build tracks with STL compliant type
//

#ifndef OPENMVG_TRACKS_H_
#define OPENMVG_TRACKS_H_

#include "lemon/list_graph.h"
#include "lemon/unionfind.h"
using namespace lemon;

#include "openMVG/matching/indMatch.hpp"

#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>
#include <set>
#include <map>

namespace openMVG  {

using namespace openMVG::matching;

/// Lightweight copy of the flat_map of BOOST library
/// Use a vector to speed up insertion (preallocated array)
template<typename T1, typename T2>
class flat_pair_map
{
  typedef std::pair<T1, T2> P;
public:
  typedef typename std::vector< P >::iterator iterator;

  typename std::vector< P >::iterator find(const T1 & val)  {
    return std::lower_bound(m_vec.begin(), m_vec.end(), val, superiorToFirst);
  }

  T2 & operator[](const T1 & val) {
    return std::lower_bound(m_vec.begin(), m_vec.end(), val, superiorToFirst)->second;
  }

  void sort()  {std::sort(m_vec.begin(), m_vec.end(), sortPairAscend);}
  void push_back(const P & val)  { m_vec.push_back(val);  }
  void clear()  { m_vec.clear(); }
  void reserve(size_t count)  { m_vec.reserve(count); }
private:
  std::vector< P > m_vec;

  static bool sortPairAscend(const P &a, const P &b) {return a.first<b.first;}
  static bool superiorToFirst(const P &a, const T1 &b) {return a.first<b;}
};

namespace tracks  {

// Data structure to store a track: collection of {ImageId,FeatureId}
//  The corresponding image points with their imageId and FeatureId.
typedef std::map<size_t,size_t> submapTrack;
// A track is a collection of {trackId, submapTrack}
typedef std::map< size_t, submapTrack > STLMAPTracks;

struct TracksBuilder
{
  typedef std::pair<size_t, size_t> indexedFeaturePair;
  typedef ListDigraph::NodeMap<size_t> IndexMap;
  typedef lemon::UnionFindEnum< IndexMap > UnionFindObject;

  typedef flat_pair_map< lemon::ListDigraph::Node, indexedFeaturePair> MapNodeToIndex;
  typedef flat_pair_map< indexedFeaturePair, lemon::ListDigraph::Node > MapIndexToNode;

  lemon::ListDigraph _graph; //Graph container to create the node
  MapNodeToIndex _map_nodeToIndex; //Node to index map
  std::auto_ptr<IndexMap> _index;
  std::auto_ptr<UnionFindObject> _tracksUF;

  const UnionFindObject & getUnionFindEnum() const {return *_tracksUF; }
  const MapNodeToIndex & getReverseMap() const {return _map_nodeToIndex;}

  /// Build tracks for a given series of pairWise matches
  bool Build( const PairWiseMatches &  map_pair_wise_matches)
  {
    typedef std::set<indexedFeaturePair> SetIndexedPair;
    // Set of all features of all images: (imageIndex, featureIndex)
    SetIndexedPair allFeatures;
    // For each couple of images
    for (PairWiseMatches::const_iterator iter = map_pair_wise_matches.begin();
      iter != map_pair_wise_matches.end();
      ++iter)
    {
      const size_t & I = iter->first.first;
      const size_t & J = iter->first.second;
      // Features correspondences between I and J image.
      const std::vector<IndMatch> & vec_FilteredMatches = iter->second;

      // Retrieve all features
      for( size_t k = 0; k < vec_FilteredMatches.size(); ++k)
      {
        allFeatures.insert(std::make_pair(I,vec_FilteredMatches[k]._i));
        allFeatures.insert(std::make_pair(J,vec_FilteredMatches[k]._j));
      }
    }

    // Build the node indirection for each referenced feature
    MapIndexToNode map_indexToNode;
    map_indexToNode.reserve(allFeatures.size());
    _map_nodeToIndex.reserve(allFeatures.size());
    for (SetIndexedPair::const_iterator iter = allFeatures.begin();
      iter != allFeatures.end();
      ++iter)
    {
      lemon::ListDigraph::Node node = _graph.addNode();
      map_indexToNode.push_back( std::make_pair(*iter, node));
      _map_nodeToIndex.push_back( std::make_pair(node,*iter));
    }

    // Sort the flat_pair_map
    map_indexToNode.sort();
    _map_nodeToIndex.sort();

    // Add the element of myset to the UnionFind insert method.
    _index = std::auto_ptr<IndexMap>( new IndexMap(_graph) );
    _tracksUF = std::auto_ptr<UnionFindObject>( new UnionFindObject(*_index));
    for (ListDigraph::NodeIt it(_graph); it != INVALID; ++it) {
      _tracksUF->insert(it);
    }

    // Make the union according the pair matches
    for (PairWiseMatches::const_iterator iter = map_pair_wise_matches.begin();
      iter != map_pair_wise_matches.end();
      ++iter)
    {
      const size_t & I = iter->first.first;
      const size_t & J = iter->first.second;
      const std::vector<IndMatch> & vec_FilteredMatches = iter->second;
      // We have correspondences between I and J image index.

      for( size_t k = 0; k < vec_FilteredMatches.size(); ++k)
      {
        indexedFeaturePair pairI(I,vec_FilteredMatches[k]._i);
        indexedFeaturePair pairJ(J,vec_FilteredMatches[k]._j);
        _tracksUF->join( map_indexToNode[pairI], map_indexToNode[pairJ] );
      }
    }
    return false;
  }

  /// Remove bad tracks (too short or track with ids collision)
  bool Filter(size_t nLengthSupTo = 2, bool bMultithread = true)
  {
    // Remove bad tracks:
    // - track that are too short,
    // - track with id conflicts (many times the same image index)

    std::set<int> set_classToErase;
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel if(bMultithread)
#endif
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*_tracksUF); cit != INVALID; ++cit) {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif
      {
        size_t cpt = 0;
        std::set<size_t> myset;
        for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*_tracksUF, cit); iit != INVALID; ++iit) {
          myset.insert(_map_nodeToIndex[ iit ].first);
          ++cpt;
        }
        if (myset.size() != cpt || myset.size() < nLengthSupTo)
        {
#ifdef OPENMVG_USE_OPENMP
          #pragma omp critical
#endif
          set_classToErase.insert(cit.operator int());
        }
      }
    }
    std::for_each (set_classToErase.begin(), set_classToErase.end(),
      std::bind1st( std::mem_fun( &UnionFindObject::eraseClass ), _tracksUF.get() ));
    return false;
  }

  /// Remove the pair that have too few correspondences.
  bool FilterPairWiseMinimumMatches(size_t minMatchesOccurences, bool bMultithread = true)
  {
    std::vector<size_t> vec_tracksToRemove;
    typedef std::map< size_t, std::set<size_t> > TrackIdPerImageT;
    TrackIdPerImageT map_tracksIdPerImages;

    //-- Count the number of track per image Id
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*_tracksUF); cit != INVALID; ++cit) {
      const size_t trackId = cit.operator int();
      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*_tracksUF, cit); iit != INVALID; ++iit) {
        const MapNodeToIndex::iterator iterTrackValue = _map_nodeToIndex.find(iit);
        const indexedFeaturePair & currentPair = iterTrackValue->second;
        map_tracksIdPerImages[currentPair.first].insert(trackId);
      }
    }

    //-- Compute corresponding track per image pair
#ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel if(bMultithread)
#endif
    for (TrackIdPerImageT::const_iterator iter = map_tracksIdPerImages.begin();
      iter != map_tracksIdPerImages.end();
      ++iter)
    {
#ifdef OPENMVG_USE_OPENMP
    #pragma omp single nowait
#endif
      {
        const std::set<size_t> & setA = iter->second;
        std::vector<size_t> inter;
        for (TrackIdPerImageT::const_iterator iter2 = iter;
          iter2 != map_tracksIdPerImages.end();  ++iter2)
        {
          // compute intersection of track ids
          const std::set<size_t> & setB = iter2->second;
          inter.clear();
          std::set_intersection(setA.begin(), setA.end(), setB.begin(), setB.end(), std::back_inserter(inter));
          if (inter.size() < minMatchesOccurences)
          {
#ifdef OPENMVG_USE_OPENMP
            #pragma omp critical
#endif
            {
              std::copy(inter.begin(), inter.end(), std::back_inserter(vec_tracksToRemove));
            }
          }
        }
      }
    }
    std::sort(vec_tracksToRemove.begin(), vec_tracksToRemove.end());
    std::vector<size_t>::iterator it = std::unique(vec_tracksToRemove.begin(), vec_tracksToRemove.end());
    vec_tracksToRemove.resize( std::distance(vec_tracksToRemove.begin(), it) );
    std::for_each(vec_tracksToRemove.begin(), vec_tracksToRemove.end(),
      std::bind1st(std::mem_fun(&UnionFindObject::eraseClass), _tracksUF.get()));
    return false;
  }

  bool ExportToStream(std::ostream & os)
  {
    size_t cpt = 0;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*_tracksUF); cit != INVALID; ++cit) {
      os << "Class: " << cpt++ << std::endl;
      size_t cptTrackLength = 0;
      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*_tracksUF, cit); iit != INVALID; ++iit) {
        ++cptTrackLength;
      }
      os << "\t" << "track length: " << cptTrackLength << std::endl;

      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*_tracksUF, cit); iit != INVALID; ++iit) {
        os << _map_nodeToIndex[ iit ].first << "  " << _map_nodeToIndex[ iit ].second << std::endl;
      }
    }
    return os.good();
  }

  /// Return the number of connected set in the UnionFind structure (tree forest)
  size_t NbTracks() const
  {
    size_t cpt = 0;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*_tracksUF); cit != INVALID; ++cit)
      ++cpt;
    return cpt;
  }

  /// Export tracks as a map (each entry is a sequence of imageId and featureIndex):
  ///  {TrackIndex => {(imageIndex, featureIndex), ... ,(imageIndex, featureIndex)}
  void ExportToSTL(STLMAPTracks & map_tracks)
  {
    map_tracks.clear();

    size_t cptClass = 0;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*_tracksUF); cit != INVALID; ++cit, ++cptClass) {
      std::pair<STLMAPTracks::iterator, bool> ret =
        map_tracks.insert(std::pair<size_t, submapTrack >(cptClass, submapTrack()));
      STLMAPTracks::iterator iterN = ret.first;

      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*_tracksUF, cit); iit != INVALID; ++iit) {
        const MapNodeToIndex::iterator iterTrackValue = _map_nodeToIndex.find(iit);
        const indexedFeaturePair & currentPair = iterTrackValue->second;

        iterN->second[currentPair.first] = currentPair.second;
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
   * @param[in] map_tracksIn: all tracks of the world
   * @param[out] map_tracksOut: output with only the common tracks
   */
  static bool GetTracksInImages(
    const std::set<size_t> & set_imageIndex,
    const STLMAPTracks & map_tracksIn,
    STLMAPTracks & map_tracksOut)
  {
    map_tracksOut.clear();

    // Go along the tracks
    for (STLMAPTracks::const_iterator iterT = map_tracksIn.begin();
      iterT != map_tracksIn.end(); ++iterT)  {

      // If the track contain one of the provided index save the point of the track
      submapTrack map_temp;
      for (std::set<size_t>::const_iterator iterIndex = set_imageIndex.begin();
        iterIndex != set_imageIndex.end(); ++iterIndex)
      {
        submapTrack::const_iterator iterSearch = iterT->second.find(*iterIndex);
        if (iterSearch != iterT->second.end())
          map_temp[iterSearch->first] = iterSearch->second;
      }

      if (!map_temp.empty() && map_temp.size() == set_imageIndex.size())
        map_tracksOut[iterT->first] = map_temp;
    }
    return !map_tracksOut.empty();
  }

  /// Return the tracksId as a set (sorted increasing)
  static void GetTracksIdVector(
    const STLMAPTracks & map_tracks,
    std::set<size_t> * set_tracksIds)
  {
    set_tracksIds->clear();
    for (STLMAPTracks::const_iterator iterT = map_tracks.begin();
      iterT != map_tracks.end(); ++iterT)
    {
      set_tracksIds->insert(iterT->first);
    }
  }

  /// Get feature index PerView and TrackId
  static bool GetFeatIndexPerViewAndTrackId(
    const STLMAPTracks & map_tracks,
    const std::set<size_t> & set_trackId,
    size_t nImageIndex,
    std::vector<size_t> * pvec_featIndex)
  {
    for (STLMAPTracks::const_iterator iterT = map_tracks.begin();
      iterT != map_tracks.end(); ++iterT)
    {
      const size_t trackId = iterT->first;
      if (set_trackId.find(trackId) != set_trackId.end())
      {
        //try to find imageIndex
        const submapTrack & map_ref = iterT->second;
        submapTrack::const_iterator iterSearch = map_ref.find(nImageIndex);
        if (iterSearch != map_ref.end())
        {
          pvec_featIndex->push_back(iterSearch->second);
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
      return ( id == val.first);
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
  static void TracksToIndexedMatches(const STLMAPTracks & map_tracks,
    const std::vector<IndexT> & vec_filterIndex,
    std::vector<IndMatch> * pvec_index)
  {

    std::vector<IndMatch> & vec_indexref = *pvec_index;
    vec_indexref.clear();
    for (size_t i = 0; i < vec_filterIndex.size(); ++i)
    {
      // Retrieve the track information from the current index i.
      STLMAPTracks::const_iterator itF =
        find_if(map_tracks.begin(), map_tracks.end(), FunctorMapFirstEqual(vec_filterIndex[i]));
      // The current track.
      const submapTrack & map_ref = itF->second;

      // We have 2 elements for a track.
      assert(map_ref.size() == 2);
      const IndexT indexI = (map_ref.begin())->second;
      const IndexT indexJ = (++map_ref.begin())->second;

      vec_indexref.push_back(IndMatch(indexI, indexJ));
    }
  }

  /// Return the occurrence of tracks length.
  static void TracksLength(const STLMAPTracks & map_tracks,
    std::map<size_t, size_t> & map_Occurence_TrackLength)
  {
    for (STLMAPTracks::const_iterator iterT = map_tracks.begin();
      iterT != map_tracks.end(); ++iterT)
    {
      const size_t trLength = iterT->second.size();
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
  static void ImageIdInTracks(const STLMAPTracks & map_tracks,
    std::set<size_t> & set_imagesId)
  {
    for (STLMAPTracks::const_iterator iterT = map_tracks.begin();
      iterT != map_tracks.end(); ++iterT)
    {
      const submapTrack & map_ref = iterT->second;
      for (submapTrack::const_iterator iter = map_ref.begin();
        iter != map_ref.end();
        ++iter)
      {
        set_imagesId.insert(iter->first);
      }
    }
  }
};

} // namespace tracks
} // namespace openMVG

#endif // OPENMVG_TRACKS_H_
