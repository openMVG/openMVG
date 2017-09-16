// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/portable_binary.hpp>

#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/indMatch_io.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <map>
#include <vector>

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

namespace openMVG {
namespace matching {

bool Load
(
  PairWiseMatches & matches,
  const std::string & filename
)
{
  matches.clear();
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "txt")
  {
    std::ifstream stream(filename.c_str());
    if (!stream.is_open())
    {
      return false;
    }
    // Read from the text file
    // I J
    // #matches count
    // idx idx
    // ...
    size_t I, J, number;
    while (stream >> I >> J >> number)  {
      std::vector<IndMatch> read_matches(number);
      for (size_t i = 0; i < number; ++i) {
        stream >> read_matches[i];
      }
      matches[{I,J}] = std::move(read_matches);
    }
    stream.close();
    return true;
  }
  else if (ext == "bin")
  {
    std::ifstream stream (filename.c_str(), std::ios::in | std::ios::binary);
    if (stream.is_open())
    {
      cereal::PortableBinaryInputArchive archive(stream);
      archive(matches);
      stream.close();
      return true;
    }
  }
  else
  {
    std::cerr << "Unknown PairWiseMatches input format: " << ext << std::endl;
  }
  return false;
}

bool Save
(
  const PairWiseMatches & matches,
  const std::string & filename
)
{
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "txt")
  {
    std::ofstream stream(filename.c_str());
    if (!stream.is_open())
    {
      return false;
    }
    for ( const auto & cur_match : matches )
    {
      const auto& I = cur_match.first.first;
      const auto& J = cur_match.first.second;

      const std::vector<IndMatch> & pair_matches = cur_match.second;
      stream << I << " " << J << '\n' << pair_matches.size() << '\n';
      copy(pair_matches.begin(), pair_matches.end(),
           std::ostream_iterator<IndMatch>(stream, "\n"));
    }
    stream.close();
    return true;
  }
  else if (ext == "bin")
  {
    std::ofstream stream (filename.c_str(), std::ios::out | std::ios::binary);
    if (stream.is_open())
    {
      cereal::PortableBinaryOutputArchive archive(stream);
      archive(matches);
      stream.close();
      return true;
    }
  }
  else
  {
    std::cerr << "Unknown PairWiseMatches output format: " << ext << std::endl;
  }
  return false;
}
}  // namespace matching
}  // namespace openMVG
