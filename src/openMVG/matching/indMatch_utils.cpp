// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/portable_binary.hpp>

#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/indMatch_io.hpp"
#include "openMVG/system/logger.hpp"

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
  std::ifstream stream;
  const std::string ext = stlplus::extension_part(filename);
  if (ext == "txt")
  {
    stream.open(filename);
    if (stream)
    {
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
      stream.clear(); // necessary since we hit eof with the while
      stream.close();
    }
  }
  else if (ext == "bin")
  {
    stream.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (stream)
    {
      cereal::PortableBinaryInputArchive archive(stream);
      archive(matches);
      stream.close();
    }
  }
  else
  {
    OPENMVG_LOG_ERROR << "Unknown PairWiseMatches file extension: (" << ext << ").";
  }

  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot open the matche file: " << filename << ".";
  }

  return static_cast<bool>(stream);
}

bool Save
(
  const PairWiseMatches & matches,
  const std::string & filename
)
{
  const std::string ext = stlplus::extension_part(filename);
  std::ofstream stream;
  if (ext == "txt")
  {
    stream.open(filename);
    if (stream)
    {
      for ( const auto & cur_match : matches )
      {
        const auto& I = cur_match.first.first;
        const auto& J = cur_match.first.second;

        const std::vector<IndMatch> & pair_matches = cur_match.second;
        stream << I << " " << J << '\n' << pair_matches.size() << '\n';
        copy(pair_matches.cbegin(), pair_matches.cend(),
             std::ostream_iterator<IndMatch>(stream, "\n"));
      }
      stream.close();
    }
  }
  else if (ext == "bin")
  {
    stream.open(filename.c_str(), std::ios::out | std::ios::binary);
    if (stream)
    {
      cereal::PortableBinaryOutputArchive archive(stream);
      archive(matches);
      stream.close();
    }
  }
  else
  {
    OPENMVG_LOG_ERROR << "Unknown PairWiseMatches output file extension: " << filename;
  }

  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot save the matche file: " << filename << ".";
  }
  return static_cast<bool>(stream);
}
}  // namespace matching
}  // namespace openMVG
