// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/translation_averaging_common.hpp"
#include "openMVG/types.hpp"

#include <utility>
#include <vector>

namespace openMVG {

// List the pairs used by the relative motions
Pair_Set getPairs(const RelativeInfo_Vec & vec_relative)
{
  Pair_Set pair_set;
  for (const auto & rel : vec_relative )
  {
    pair_set.insert(Pair(rel.first.first, rel.first.second));
  }
  return pair_set;
}

Pair_Set getPairs(const std::vector<RelativeInfo_Vec> & vec_relative)
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
std::set<IndexT> getIndexT(const RelativeInfo_Vec & vec_relative)
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
std::set<IndexT> getIndexT(const std::vector<RelativeInfo_Vec> & vec_relative)
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

