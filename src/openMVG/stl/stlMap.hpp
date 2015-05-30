
/**
 * @file stlMap.h
 * @brief Functor to select key or value of stl maps
 * @author Pierre MOULON
 *
 * Copyright (c) 2011, 2012, 2013 Pierre MOULON
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef STL_MAP_ADDITION_H
#define STL_MAP_ADDITION_H

// ---------------------------
// Usage example :
// ---------------------------

// std::map<int, double> m;
// m[0] = 2.0;
// m[1] = 1.6;
// std::vector<int> keys;
// // Retrieve all keys
// transform(m.begin(), m.end(), back_inserter(keys), stl::RetrieveKey());
// // Log all keys to console
// copy(keys.begin(), keys.end(), ostream_iterator<int>(cout, "\n"));
// // Retrieve all values
// std::vector<double> values;
// // Retrieve all values
// transform(m.begin(), m.end(), back_inserter(values), stl::RetrieveValue());

#include <map>

namespace stl{

/// Allow to select the Keys of a map.
struct RetrieveKey
{
  template <typename T>
  typename T::first_type operator()(const T & keyValuePair) const
  {
    return keyValuePair.first;
  }
};

/// Allow to select the Values of a map.
struct RetrieveValue
{
  template <typename T>
  typename T::second_type operator()(const T & keyValuePair) const
  {
    return keyValuePair.second;
  }
};

} // namespace stl

#endif // STL_MAP_ADDITION_H

