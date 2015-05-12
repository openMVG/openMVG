
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_UTILS_H
#define OPENMVG_MATCHING_IND_MATCH_UTILS_H

#include "openMVG/matching/indMatch.hpp"
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace openMVG {
namespace matching {

/// Export vector of IndMatch to a stream
static bool PairedIndMatchToStream(
  const PairWiseMatches & map_indexedMatches,
  std::ostream & os)
{
  for (PairWiseMatches::const_iterator iter = map_indexedMatches.begin();
    iter != map_indexedMatches.end();
    ++iter)
  {
    const size_t I = iter->first.first;
    const size_t J = iter->first.second;
    const std::vector<IndMatch> & vec_matches = iter->second;
    os << I << " " << J << '\n' << vec_matches.size() << '\n';
    copy(vec_matches.begin(), vec_matches.end(),
         std::ostream_iterator<IndMatch>(os, "\n"));
  }
  return os.good();
}

/// Import vector of IndMatch from a file
static bool PairedIndMatchImport(
  const std::string & fileName,
  PairWiseMatches & map_indexedMatches)
{
  bool bOk = false;
  std::ifstream in(fileName.c_str());
  if (in.is_open()) {
    map_indexedMatches.clear();

    size_t I, J, number;
    while (in >> I >> J >> number)  {
      std::vector<IndMatch> matches(number);
      for (size_t i = 0; i < number; ++i) {
        in >> matches[i];
      }
      map_indexedMatches[std::make_pair(I,J)] = matches;
    }
    bOk = true;
  }
  else  {
    std::cout << std::endl << "ERROR indexedMatchesUtils::import(...)" << std::endl
      << "with : " << fileName << std::endl;
    bOk = false;
  }
  return bOk;
}
}  // namespace matching
}  // namespace openMVG

#endif // #define OPENMVG_MATCHING_IND_MATCH_UTILS_H
