// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_STL_HASH_HPP
#define OPENMVG_STL_HASH_HPP

#include <functional>

namespace stl
{

// Combine hashing value
// http://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine
template <class T, class Hasher>
inline void hash_combine(std::size_t& seed, const T& v, Hasher hasher)
{
  seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
  stl::hash_combine(seed, v, std::hash<T>());
}

} // namespace stl

#endif  // OPENMVG_STL_HASH_HPP
