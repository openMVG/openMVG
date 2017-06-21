// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_TRACKS_UNION_FIND_DISJOINT_SET_HPP
#define OPENMVG_TRACKS_UNION_FIND_DISJOINT_SET_HPP

#include <numeric>
#include <vector>

namespace openMVG  {

// Union-Find/Disjoint-Set data structure
//--
// A disjoint-set data structure also called a union–find data structure
// or merge–find set, is a data structure that keeps track of a set of elements
// partitioned into a number of disjoint (non-overlapping) subsets.
// It supports two operations:
// - Find: Determine which subset a particular element is in.
//   - It returns an item from this set that serves as its "representative";
// - Union: Join two subsets into a single subset.
// Sometime a Connected method is implemented:
// - Connected:
//   - By comparing the result of two Find operations, one can determine whether
//      two elements are in the same subset.
//--
struct UnionFind
{
  // Represent the DS/UF forest thanks to two array:
  // A parent 'pointer tree' where each node holds a reference to its parent node
  std::vector<unsigned int> m_cc_parent;
  // A rank array used for union by rank
  std::vector<unsigned int> m_cc_rank;
  // A 'size array' to know the size of each connected component
  std::vector<unsigned int> m_cc_size;

  // Init the UF structure with num_cc nodes
  void InitSets
  (
    const unsigned int num_cc
  )
  {
    // all set size are 1 (independent nodes)
    m_cc_size.resize(num_cc, 1);
    // Parents id have their own CC id {0,n}
    m_cc_parent.resize(num_cc);
    std::iota(m_cc_parent.begin(), m_cc_parent.end(), 0);
    // Rank array (0)
    m_cc_rank.resize(num_cc, 0);
  }

  // Return the number of nodes that have been initialized in the UF tree
  unsigned int GetNumNodes() const
  {
    return static_cast<unsigned int>(m_cc_size.size());
  }

  // Return the representative set id of I nth component
  unsigned int Find
  (
    unsigned int i
  )
  {
    // Recursively set all branch as children of root (Path compression)
    if (m_cc_parent[i] != i)
        m_cc_parent[i] = Find(m_cc_parent[i]);
     return m_cc_parent[i];
  }

  // Replace sets containing I and J with their union
  void Union
  (
    unsigned int i,
    unsigned int j
  )
  {
    i = Find(i);
    j = Find(j);
    if (i == j)
    { // Already in the same set. Nothing to do
      return;
    }

    // x and y are not already in same set. Merge them.
    // Perform an union by rank:
    //  - always attach the smaller tree to the root of the larger tree
    if (m_cc_rank[i] < m_cc_rank[j])
    {
      m_cc_parent[i] = j;
      m_cc_size[j] += m_cc_size[i];
    }
    else
    {
      m_cc_parent[j] = i;
      m_cc_size[i] += m_cc_size[j];
      if (m_cc_rank[i] == m_cc_rank[j])
        ++m_cc_rank[i];
    }
  }
};

} // namespace openMVG

#endif // OPENMVG_TRACKS_UNION_FIND_DISJOINT_SET_HPP
