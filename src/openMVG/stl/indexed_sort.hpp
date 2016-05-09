
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_STL_INDEXED_SORT_H
#define OPENMVG_STL_INDEXED_SORT_H

#include <vector>

namespace stl
{
namespace indexed_sort
{
  template<typename T1, typename T2>
  struct sort_index_packet_ascend {
    T1 val;
    T2 index;
  };

  template<typename T1, typename T2>
  struct sort_index_packet_descend  {
    T1 val;
    T2 index;
  };

  template<typename T1, typename T2>
  inline
    bool
    operator< (const sort_index_packet_ascend<T1,T2>& A,
                const sort_index_packet_ascend<T1,T2>& B) {
    return A.val < B.val;
  }

  template<typename T1, typename T2>
  inline
    bool
    operator< (const sort_index_packet_descend<T1,T2>& A,
                const sort_index_packet_descend<T1,T2>& B)  {
    return A.val > B.val;
  }

  /// Sort by default all indexed value, else sort only the NN smallest element of the indexed array.
  template<typename packet_type, typename eT>
  inline void sort_index_helper(std::vector<packet_type>& packet_vec,
                      const eT* in_mem, int NN = -1)  {
    const size_t n_elem = packet_vec.size();

    for(size_t i=0; i<n_elem; ++i)  {
      packet_vec[i].val   = in_mem[i];
      packet_vec[i].index = i;
    }

    if (NN == -1)
      std::sort( packet_vec.begin(), packet_vec.end() );
    else
      std::partial_sort(packet_vec.begin(), packet_vec.begin() + NN,
        packet_vec.end());
  }

} // namespace indexed_sort
} // namespace stl

#endif  // OPENMVG_STL_INDEXED_SORT_H
