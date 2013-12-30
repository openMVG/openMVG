
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
//  typedef std::map< std::pair<size_t, size_t>, std::vector<IndMatch> > map_pairWiseMatches;
//  map_pairWiseMatches map_Matches;
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
using namespace openMVG::matching;

#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>
#include <set>
#include <map>

/// Lightweight copy of the flat_map of BOOST library
/// Use a vector to speed up insertion (preallocated array)
template<typename T1, typename T2>
class flat_pair_map
{
  typedef std::pair<T1, T2> P;
public:
  typedef typename std::vector< P >::iterator iterator;

  typename std::vector< P >::iterator find(const T1 & val)  {
    return lower_bound(m_vec.begin(), m_vec.end(), val, superiorToFirst);
  }

  T2 & operator[](const T1 & val) {
    return lower_bound(m_vec.begin(), m_vec.end(), val, superiorToFirst)->second;
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


namespace openMVG  {
namespace tracks  {
  using namespace std;

// Pairwise matches (indexed matches for a pair <I,J>)
typedef std::map< std::pair<size_t,size_t>, std::vector<IndMatch> > mapPairWiseMatches;

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

  typedef flat_pair_map< lemon::ListDigraph::Node, indexedFeaturePair> MapNodeIndex;
  typedef flat_pair_map< indexedFeaturePair, lemon::ListDigraph::Node > MapIndexNode;

  lemon::ListDigraph g; //Graph container to create the node
  MapNodeIndex reverse_my_Map; //Node to index map
  auto_ptr<IndexMap> index;
  auto_ptr<UnionFindObject> myTracksUF;

  const UnionFindObject & getUnionFindEnum() const {return *myTracksUF; }
  const MapNodeIndex & getReverseMap() const {return reverse_my_Map;}

  /// Build tracks for a given series of pairWise matches
  bool Build( const mapPairWiseMatches &  map_pair_wise_matches)
  {
    typedef std::set<indexedFeaturePair> SetIndexedPair;
    SetIndexedPair myset;
    for (mapPairWiseMatches::const_iterator iter = map_pair_wise_matches.begin();
      iter != map_pair_wise_matches.end();
      ++iter)
    {
      const size_t & I = iter->first.first;
      const size_t & J = iter->first.second;
      const std::vector<IndMatch> & vec_FilteredMatches = iter->second;
      // We have correspondences between I and J image index.

      for( size_t k = 0; k < vec_FilteredMatches.size(); ++k)
      {
        // Look if one of the feature already belong to a track :
        myset.insert(make_pair(I,vec_FilteredMatches[k]._i));
        myset.insert(make_pair(J,vec_FilteredMatches[k]._j));
      }
    }

    // Build the node indirection for each referenced feature
    MapIndexNode my_Map;
    my_Map.reserve(myset.size());
    reverse_my_Map.reserve(myset.size());
    for (SetIndexedPair::const_iterator iter = myset.begin();
      iter != myset.end();
      ++iter)
    {
      lemon::ListDigraph::Node node = g.addNode();
      my_Map.push_back( std::make_pair(*iter, node));
      reverse_my_Map.push_back( std::make_pair(node,*iter));
    }

    // Sort the flat_pair_map
    my_Map.sort();
    reverse_my_Map.sort();

    // Add the element of myset to the UnionFind insert method.
    index = auto_ptr<IndexMap>( new IndexMap(g) );
    myTracksUF = auto_ptr<UnionFindObject>( new UnionFindObject(*index));
    for (ListDigraph::NodeIt it(g); it != INVALID; ++it) {
      myTracksUF->insert(it);
    }

    // Make the union according the pair matches
    for (mapPairWiseMatches::const_iterator iter = map_pair_wise_matches.begin();
      iter != map_pair_wise_matches.end();
      ++iter)
    {
      const size_t & I = iter->first.first;
      const size_t & J = iter->first.second;
      const vector<IndMatch> & vec_FilteredMatches = iter->second;
      // We have correspondences between I and J image index.

      for( size_t k = 0; k < vec_FilteredMatches.size(); ++k)
      {
        indexedFeaturePair pairI = make_pair(I,vec_FilteredMatches[k]._i);
        indexedFeaturePair pairJ = make_pair(J,vec_FilteredMatches[k]._j);
        myTracksUF->join( my_Map[pairI], my_Map[pairJ] );
      }
    }
    return false;
  }

  /// Remove bad tracks, conflict tracks (many times the same image index in a track)
  bool Filter(size_t nLengthSupTo = 2)
  {
    // Remove bad tracks (shorter, conflicts (Many times the same image index)...)

    // Remove tracks that have a conflict (many times the same image index)
    std::set<int> set_classToErase;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*myTracksUF); cit != INVALID; ++cit) {
      size_t cpt = 0;
      std::set<size_t> myset;
      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*myTracksUF, cit); iit != INVALID; ++iit) {
        myset.insert(reverse_my_Map[ iit ].first);
        ++cpt;
      }
      if (myset.size() != cpt || myset.size() < nLengthSupTo)
      {
        set_classToErase.insert(cit.operator int());
      }
    }
    for_each (set_classToErase.begin(), set_classToErase.end(),
      std::bind1st( std::mem_fun( &UnionFindObject::eraseClass ), myTracksUF.get() ));
    return false;
  }

  /// Remove the pair that have too few correspondences.
  bool FilterPairWiseMinimumMatches(size_t minMatchesOccurences, bool bVerbose = false)
  {
    std::vector<size_t> vec_tracksToRemove;
    std::map< size_t, set<size_t> > map_tracksIdPerImages;

    //-- Count the number of track per image Id
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*myTracksUF); cit != INVALID; ++cit) {
      size_t trackId = cit.operator int();
      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*myTracksUF, cit); iit != INVALID; ++iit) {
        const MapNodeIndex::iterator iterTrackValue = reverse_my_Map.find(iit);
        const indexedFeaturePair & currentPair = iterTrackValue->second;
        map_tracksIdPerImages[currentPair.first].insert(trackId);
      }
    }

    //-- Compute cross images matches
    for (std::map<size_t, set<size_t> >::const_iterator iter = map_tracksIdPerImages.begin();
      iter != map_tracksIdPerImages.end();
      ++iter)
    {
      const set<size_t> & setA = iter->second;
      for (std::map<size_t, set<size_t> >::const_iterator iter2 = iter;
        iter2 != map_tracksIdPerImages.end();  ++iter2)
      {
        const set<size_t> & setB = iter2->second;
        vector<size_t> inter;

        set_intersection (setA.begin(), setA.end(), setB.begin(), setB.end(), back_inserter(inter));

        if (inter.size() < minMatchesOccurences)
          copy(inter.begin(), inter.end(), back_inserter(vec_tracksToRemove));
      }
    }
    sort(vec_tracksToRemove.begin(), vec_tracksToRemove.end());
    std::vector<size_t>::iterator it = unique (vec_tracksToRemove.begin(), vec_tracksToRemove.end());
    vec_tracksToRemove.resize( std::distance(vec_tracksToRemove.begin(), it) );
    if (bVerbose)
      cout << endl << endl <<vec_tracksToRemove.size() << " Tracks will be removed"<< endl;
    for_each (vec_tracksToRemove.begin(), vec_tracksToRemove.end(),
      std::bind1st( std::mem_fun( &UnionFindObject::eraseClass ), myTracksUF.get() ));
    return false;
  }

  bool ExportToStream(ostream & os)
  {
    size_t cpt = 0;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*myTracksUF); cit != INVALID; ++cit) {
      os << "Class: " << cpt++ << std::endl;
      size_t cptTrackLength = 0;
      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*myTracksUF, cit); iit != INVALID; ++iit) {
        ++cptTrackLength;
      }
      os << "\t" << "track length: " << cptTrackLength << std::endl;

      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*myTracksUF, cit); iit != INVALID; ++iit) {
        os << reverse_my_Map[ iit ].first << "  " << reverse_my_Map[ iit ].second << std::endl;
      }
    }
    return os.good();
  }

  /// Return the number of connected set in the UnionFind structure (tree forest)
  size_t NbTracks() const
  {
    size_t cpt = 0;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*myTracksUF); cit != INVALID; ++cit)
      ++cpt;
    return cpt;
  }

  /// Export tracks as a map (each entry is a sequence of imageId and featureIndex):
  ///  {TrackIndex => {(imageIndex, featureIndex), ... ,(imageIndex, featureIndex)}
  void ExportToSTL(STLMAPTracks & map_tracks)
  {
    map_tracks.clear();

    size_t cptClass = 0;
    for ( lemon::UnionFindEnum< IndexMap >::ClassIt cit(*myTracksUF); cit != INVALID; ++cit, ++cptClass) {
      pair<STLMAPTracks::iterator, bool> ret =
        map_tracks.insert(std::pair<size_t, submapTrack >(cptClass, submapTrack() ) );
      STLMAPTracks::iterator iterN = ret.first;

      for (lemon::UnionFindEnum< IndexMap >::ItemIt iit(*myTracksUF, cit); iit != INVALID; ++iit) {
        const MapNodeIndex::iterator iterTrackValue = reverse_my_Map.find(iit);
        const indexedFeaturePair & currentPair = iterTrackValue->second;

        iterN->second[currentPair.first] = currentPair.second;
      }
    }
  }
};

struct TracksUtilsMap
{
  /// Return the tracks that are in common to the set_imageIndex indexes.
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
      std::map<size_t, size_t> map_temp;
      for (std::set<size_t>::const_iterator iterIndex = set_imageIndex.begin();
        iterIndex != set_imageIndex.end(); ++iterIndex)
      {
        submapTrack::const_iterator iterSearch = iterT->second.find(*iterIndex);
        if (iterSearch != iterT->second.end())
          map_temp[iterSearch->first] = iterSearch->second;
      }

      if (!map_temp.empty() && map_temp.size() == set_imageIndex.size())
        map_tracksOut.insert(make_pair(iterT->first, map_temp));
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
      size_t trackId = iterT->first;
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

  struct FunctorMapFirstEqual : public std::unary_function <std::map<size_t, std::map<size_t,size_t> >, bool>
  {
    size_t id;
    FunctorMapFirstEqual(size_t val):id(val){};
    bool operator()(const std::pair<size_t, std::map<size_t,size_t> > & val) {
      return ( id == val.first);
    }
  };

  /// Convert a trackId to a vector of indexed Matches.
  /// The input tracks must be compound of only two images index.
  /// Be careful image index are sorted (increasing order)
  /// Only track index contained in the filter vector are kept.
  static void TracksToIndexedMatches(const STLMAPTracks & map_tracks,
    const std::vector<size_t> & vec_filterIndex,
    std::vector<IndMatch> * pvec_index)
  {

    std::vector<IndMatch> & vec_indexref = *pvec_index;
    vec_indexref.clear();
    for (size_t i = 0; i < vec_filterIndex.size(); ++i)
    {
      STLMAPTracks::const_iterator itF =
        find_if(map_tracks.begin(), map_tracks.end(), FunctorMapFirstEqual(vec_filterIndex[i]));
      const submapTrack & map_ref = itF->second;
      submapTrack::const_iterator iter = map_ref.begin();
      size_t indexI = iter->second;
      ++iter;
      size_t indexJ = iter->second;
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
      size_t trLength = iterT->second.size();
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
