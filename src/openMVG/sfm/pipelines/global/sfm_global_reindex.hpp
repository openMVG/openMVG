// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_REINDEX_HPP
#define OPENMVG_SFM_GLOBAL_REINDEX_HPP

#include <set>

namespace openMVG {
namespace sfm{

/// Association of Ids to a contiguous set of Ids
template<typename IterablePairs, typename PairValueType>
void reindex
(
  const IterablePairs& pairs,
  Hash_Map<PairValueType, PairValueType> & reindex_forward,
  Hash_Map<PairValueType, PairValueType> & reindex_backward
)
{
  // get a unique set of Ids
  std::set<size_t> unique_id;
  for (typename IterablePairs::const_iterator iter = pairs.begin();
        iter != pairs.end(); ++iter)
  {
    unique_id.insert(iter->first);
    unique_id.insert(iter->second);
  }

  // Build the Forward and Backward mapping
  for (typename IterablePairs::const_iterator iter = pairs.begin();
        iter != pairs.end(); ++iter)
  {
    if (reindex_forward.find(iter->first) == reindex_forward.end())
    {
      const size_t dist = std::distance(unique_id.begin(), unique_id.find(iter->first));
      reindex_forward[iter->first] = dist;
      reindex_backward[dist] = iter->first;
    }
    if (reindex_forward.find(iter->second) == reindex_forward.end())
    {
      const size_t dist = std::distance(unique_id.begin(), unique_id.find(iter->second));
      reindex_forward[iter->second] = dist;
      reindex_backward[dist] = iter->second;
    }
  }
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_REINDEX_HPP
