// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_FILTERS_HPP
#define OPENMVG_SFM_FILTERS_HPP

#include "openMVG/multiview/rotation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/matching/indMatch.hpp"

namespace openMVG {
namespace sfm {

template<typename IterableIndexTSequence>
inline std::set<IndexT> getIndexes(const IterableIndexTSequence & seq)
{
  std::set<IndexT> setOut;
  for(typename IterableIndexTSequence::const_iterator it = seq.begin(); it != seq.end(); ++it)
    setOut.insert(it->first);
  return setOut;
}

/// Filter the toFilter iterable sequence (keep only the element that share a common index
///  with the provided Ids index list).
template<typename T>
inline void KeepOnlyReferencedElement(
  const std::set<IndexT> & Ids,
  T & toFilter)
{
  std::cout << "Must be specialized for your type" << std::endl;
}

// Specialization for RelativeInfo_Map
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  RelativeInfo_Map& map_relatives)
{
  RelativeInfo_Map map_relatives_infered;
  for (RelativeInfo_Map::const_iterator
    iter = map_relatives.begin();
    iter != map_relatives.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first.first) != set_remainingIds.end() &&
        set_remainingIds.find(iter->first.second) != set_remainingIds.end())
    {
      map_relatives_infered.insert(*iter);
    }
  }
  map_relatives.swap(map_relatives_infered);
}

// Specialization for RelativeInfo_Map
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  rotation_averaging::RelativeRotations& relative_info)
{
  rotation_averaging::RelativeRotations relatives_infered;
  for (rotation_averaging::RelativeRotations::const_iterator
    iter = relative_info.begin();
    iter != relative_info.end(); ++iter)
  {
    if (set_remainingIds.find(iter->i) != set_remainingIds.end() &&
        set_remainingIds.find(iter->j) != set_remainingIds.end())
    {
      relatives_infered.push_back(*iter);
    }
  }
  relative_info.swap(relatives_infered);
}

// Specialization for PairWiseMatches
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  openMVG::matching::PairWiseMatches& map_matches)
{
  openMVG::matching::PairWiseMatches map_matches_E_infered;
  for (openMVG::matching::PairWiseMatches::const_iterator iter = map_matches.begin();
    iter != map_matches.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first.first) != set_remainingIds.end() &&
        set_remainingIds.find(iter->first.second) != set_remainingIds.end())
    {
      map_matches_E_infered.insert(*iter);
    }
  }
  map_matches.swap(map_matches_E_infered);
}

// Specialization for std::map<IndexT,Mat3>
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  std::map<IndexT,Mat3>& map_Mat3)
{
  std::map<IndexT,Mat3> map_infered;
  for (std::map<IndexT,Mat3>::const_iterator iter = map_Mat3.begin();
    iter != map_Mat3.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first) != set_remainingIds.end())
    {
      map_infered.insert(*iter);
    }
  }
  map_Mat3.swap(map_infered);
}

// Specialization for RelativeInfo_Vec
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  RelativeInfo_Vec & relativeInfo_vec)
{
  RelativeInfo_Vec map_infered;
  for (RelativeInfo_Vec::const_iterator iter = relativeInfo_vec.begin();
    iter != relativeInfo_vec.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first.first) != set_remainingIds.end() &&
        set_remainingIds.find(iter->first.second) != set_remainingIds.end())
    {
      map_infered.push_back(*iter);
    }
  }
  relativeInfo_vec.swap(map_infered);
}

// Specialization for std::vector<RelativeInfo_Vec>
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  std::vector<RelativeInfo_Vec> & relative_motion_group)
{
  std::vector<RelativeInfo_Vec> infered_relative_motion;
  for (const openMVG::RelativeInfo_Vec & iter : relative_motion_group)
  {
    RelativeInfo_Vec group;
    for (const relativeInfo & rel : iter)
    {
      if (set_remainingIds.find(rel.first.first) != set_remainingIds.end() &&
          set_remainingIds.find(rel.first.second) != set_remainingIds.end())
      {
        group.push_back(rel);
      }
    }
    if (!group.empty())
    {
      infered_relative_motion.push_back(std::move(group));
    }
  }
  relative_motion_group.swap(infered_relative_motion);
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_FILTERS_HPP
