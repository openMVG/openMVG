// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_PAIR_BUILDER_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_PAIR_BUILDER_HPP

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "openMVG/types.hpp"
#include "openMVG/stl/split.hpp"

namespace openMVG {

/// Generate all the (I,J) pairs of the upper diagonal of the NxN matrix
inline Pair_Set exhaustivePairs(const size_t N)
{
  Pair_Set pairs;
  for (IndexT I = 0; I < static_cast<IndexT>(N); ++I)
    for (IndexT J = I+1; J < static_cast<IndexT>(N); ++J)
      pairs.insert({I,J});

  return pairs;
}

/// Generate the pairs that have a distance inferior to the overlapSize
/// Usable to match video sequence
inline Pair_Set contiguousWithOverlap(const size_t N, const size_t overlapSize)
{
  Pair_Set pairs;
  for (IndexT I = 0; I < static_cast<IndexT>(N); ++I)
    for (IndexT J = I+1; J < I+1+overlapSize && J < static_cast<IndexT>(N); ++J)
      pairs.insert({I,J});
  return pairs;
}

/// Load a set of Pair_Set from a file
/// I J K L (pair that link I)
inline bool loadPairs(
     const size_t N,  // number of image in the current project (to check index validity)
     const std::string &sFileName, // filename of the list file,
     Pair_Set & pairs)  // output pairs read from the list file
{
  std::ifstream in(sFileName.c_str());
  if (!in.is_open())
  {
    std::cerr << std::endl
      << "loadPairs: Impossible to read the specified file: \"" << sFileName << "\"." << std::endl;
    return false;
  }
  std::string sValue;
  std::vector<std::string> vec_str;
  while (std::getline( in, sValue ) )
  {
    vec_str.clear();
    stl::split(sValue, ' ', vec_str);
    const IndexT str_size (vec_str.size());
    if (str_size < 2)
    {
      std::cerr << "loadPairs: Invalid input file: \"" << sFileName << "\"." << std::endl;
      return false;
    }
    std::stringstream oss;
    oss.clear(); oss.str(vec_str[0]);
    IndexT I, J;
    oss >> I;
    for (IndexT i=1; i<str_size; ++i)
    {
      oss.clear(); oss.str(vec_str[i]);
      oss >> J;
      if ( I > N-1 || J > N-1) //I&J always > 0 since we use unsigned type
      {
        std::cerr << "loadPairs: Invalid input file. Image out of range. "
                << "I: " << I << " J:" << J << " N:" << N << std::endl
                << "File: \"" << sFileName << "\"." << std::endl;
        return false;
      }
      if ( I == J )
      {
        std::cerr << "loadPairs: Invalid input file. Image " << I << " see itself. File: \"" << sFileName << "\"." << std::endl;
        return false;
      }
      // Insert the pair such that .first < .second
      pairs.insert( {std::min(I, J), std::max(I,J)} );
    }
  }
  in.close();
  return true;
}

/// Save a set of Pair_Set to a file (one pair per line)
/// I J
/// I K
/// ...
inline bool savePairs(const std::string &sFileName, const Pair_Set & pairs)
{
  std::ofstream outStream(sFileName.c_str());
  if (!outStream.is_open())  {
    std::cerr << std::endl
      << "savePairs: Impossible to open the output specified file: \"" << sFileName << "\"." << std::endl;
    return false;
  }
  for ( const auto & cur_pair : pairs )
  {
    outStream << cur_pair.first << ' ' << cur_pair.second << '\n';
  }
  const bool bOk = !outStream.bad();
  outStream.close();
  return bOk;
}

} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_PAIR_BUILDER_HPP
