// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_TYPES_HPP
#define OPENMVG_TYPES_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "openMVG/stl/hash.hpp"

/**
* @brief Main namespace of openMVG API
*/
namespace openMVG
{

namespace detail
{
// Preserve the standard hasher for keys that do not need OpenMVG's pair hash.
template<typename Key>
struct HashSelector
{
  using type = std::hash<Key>;
};

template<typename T1, typename T2>
struct HashSelector<std::pair<T1, T2>>
{
  struct type
  {
    std::size_t operator()(const std::pair<T1, T2> &p) const
    {
      std::size_t seed1(0);
      stl::hash_combine(seed1, p.first, typename HashSelector<T1>::type());
      stl::hash_combine(seed1, p.second, typename HashSelector<T2>::type());

      std::size_t seed2(0);
      stl::hash_combine(seed2, p.second, typename HashSelector<T2>::type());
      stl::hash_combine(seed2, p.first, typename HashSelector<T1>::type());

      return std::min(seed1, seed2);
    }
  };
};
} // namespace detail

/// Portable type used to store an index
using IndexT = uint32_t;

/// Portable value used to save an undefined index value
static const IndexT UndefinedIndexT = std::numeric_limits<IndexT>::max();

/// Standard Pair of IndexT
using Pair = std::pair<IndexT, IndexT>;

/// Set of Pair
using Pair_Set = std::set<Pair>;

/// Vector of Pair
using Pair_Vec = std::vector<Pair>;

/**
* @brief Standard Hash_Map class
* @tparam K type of the keys
* @tparam V type of the values
*/
template<typename Key, typename Value>
using Hash_Map = std::unordered_map<Key, Value, typename detail::HashSelector<Key>::type>;

} // namespace openMVG

#endif  // OPENMVG_TYPES_HPP
