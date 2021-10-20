// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_FILTERS_HPP
#define OPENMVG_SFM_SFM_FILTERS_HPP

#include <map>
#include <set>
#include <vector>

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/rotation_averaging_common.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/system/logger.hpp"

namespace openMVG {
namespace sfm {

template<typename IterableIndexTSequence>
inline std::set<IndexT> getIndexes(const IterableIndexTSequence & seq)
{
  std::set<IndexT> setOut;
  for (typename IterableIndexTSequence::const_iterator it = seq.begin(); it != seq.end(); ++it)
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
  OPENMVG_LOG_ERROR << "Must be specialized for your type";
}

// Specialization for RelativeInfo_Map
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  RelativeInfo_Map& map_relatives)
{
  RelativeInfo_Map map_relatives_inferred;
  for (RelativeInfo_Map::const_iterator
    iter = map_relatives.begin();
    iter != map_relatives.end(); ++iter)
  {
    if (set_remainingIds.count(iter->first.first) &&
        set_remainingIds.count(iter->first.second))
    {
      map_relatives_inferred.insert(*iter);
    }
  }
  map_relatives.swap(map_relatives_inferred);
}

// Specialization for RelativeInfo_Map
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  rotation_averaging::RelativeRotations& relative_info)
{
  rotation_averaging::RelativeRotations relatives_inferred;
  for (rotation_averaging::RelativeRotations::const_iterator
    iter = relative_info.begin();
    iter != relative_info.end(); ++iter)
  {
    if (set_remainingIds.count(iter->i) &&
        set_remainingIds.count(iter->j))
    {
      relatives_inferred.push_back(*iter);
    }
  }
  relative_info.swap(relatives_inferred);
}

// Specialization for PairWiseMatches
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  openMVG::matching::PairWiseMatches& map_matches)
{
  openMVG::matching::PairWiseMatches map_matches_E_inferred;
  for (openMVG::matching::PairWiseMatches::const_iterator iter = map_matches.begin();
    iter != map_matches.end(); ++iter)
  {
    if (set_remainingIds.count(iter->first.first) &&
        set_remainingIds.count(iter->first.second))
    {
      map_matches_E_inferred.insert(*iter);
    }
  }
  map_matches.swap(map_matches_E_inferred);
}

// Specialization for std::map<IndexT,Mat3>
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  std::map<IndexT,Mat3>& map_Mat3)
{
  std::map<IndexT,Mat3> map_inferred;
  for (std::map<IndexT,Mat3>::const_iterator iter = map_Mat3.begin();
    iter != map_Mat3.end(); ++iter)
  {
    if (set_remainingIds.count(iter->first))
    {
      map_inferred.insert(*iter);
    }
  }
  map_Mat3.swap(map_inferred);
}

// Specialization for RelativeInfo_Vec
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  RelativeInfo_Vec & relativeInfo_vec)
{
  RelativeInfo_Vec map_inferred;
  for (RelativeInfo_Vec::const_iterator iter = relativeInfo_vec.begin();
    iter != relativeInfo_vec.end(); ++iter)
  {
    if (set_remainingIds.count(iter->first.first) &&
        set_remainingIds.count(iter->first.second))
    {
      map_inferred.push_back(*iter);
    }
  }
  relativeInfo_vec.swap(map_inferred);
}

// Specialization for std::vector<RelativeInfo_Vec>
template<>
inline
void KeepOnlyReferencedElement(
  const std::set<IndexT> & set_remainingIds,
  std::vector<RelativeInfo_Vec> & relative_motion_group)
{
  std::vector<RelativeInfo_Vec> inferred_relative_motion;
  for (const openMVG::RelativeInfo_Vec & iter : relative_motion_group)
  {
    RelativeInfo_Vec group;
    for (const relativeInfo & rel : iter)
    {
      if (set_remainingIds.count(rel.first.first) &&
          set_remainingIds.count(rel.first.second))
      {
        group.push_back(rel);
      }
    }
    if (!group.empty())
    {
      inferred_relative_motion.push_back(std::move(group));
    }
  }
  relative_motion_group.swap(inferred_relative_motion);
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_FILTERS_HPP
