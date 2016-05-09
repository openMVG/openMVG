// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_
#define OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_

#include "openMVG/types.hpp"
#include "openMVG/numeric/numeric.h"

#include <utility>
#include <vector>

namespace openMVG {

/// Relative information [Rij|tij] for a pair
typedef std::pair< Pair, std::pair<Mat3,Vec3> > relativeInfo;

typedef std::vector< relativeInfo > RelativeInfo_Vec;
typedef std::map< Pair, std::pair<Mat3, Vec3> > RelativeInfo_Map;

// List the pairs used by the relative motions
inline Pair_Set getPairs(const RelativeInfo_Vec & vec_relative)
{
  Pair_Set pair_set;
  for( const auto & rel : vec_relative ) 
  {
    pair_set.insert(Pair(rel.first.first, rel.first.second));
  }
  return pair_set;
}

inline Pair_Set getPairs(const std::vector<RelativeInfo_Vec> & vec_relative)
{
  Pair_Set pair_set;
  for (const auto & it : vec_relative)
  {
    for (const relativeInfo & iter : it)
    {
      pair_set.insert(Pair(iter.first.first, iter.first.second));
    }
  }

  return pair_set;
}

// List the index used by the relative motions
inline std::set<IndexT> getIndexT(const RelativeInfo_Vec & vec_relative)
{
  std::set<IndexT> indexT_set;
  for ( const auto & rel : vec_relative ) 
  {
    indexT_set.insert(rel.first.first);
    indexT_set.insert(rel.first.second);
  }
  return indexT_set;
}

// List the index used by the relative motions
inline std::set<IndexT> getIndexT(const std::vector<RelativeInfo_Vec> & vec_relative)
{
  std::set<IndexT> indexT_set;
  for (const auto & it : vec_relative)
  {
    for (const relativeInfo & iter : it)
    {
      indexT_set.insert(iter.first.first);
      indexT_set.insert(iter.first.second);
    }
  }
  return indexT_set;
}


} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_TRANSLATION_AVERAGING_COMMON_H_
