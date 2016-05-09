// Copyright (c) 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SPLIT_HPP
#define SPLIT_HPP

#include <sstream>
#include <string>
#include <vector>

namespace stl
{
/**
 * Split an input string with a delimiter and fill a string vector
 */
inline bool split
(
  const std::string rhs,
  const char delim,
  std::vector<std::string> & items
)
{
  items.clear();
  std::stringstream ss(rhs);
  std::string item;
  while (std::getline(ss, item, delim))
  {
    items.emplace_back(item);
  }

  // return true if the delimiter is present in the input string
  return ( rhs.find(delim) != std::string::npos );
}
} // namespace stl
#endif // SPLIT_HPP
