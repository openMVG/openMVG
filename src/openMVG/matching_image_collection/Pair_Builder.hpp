// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/stl/split.hpp"

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

namespace openMVG {

/// Generate all the (I,J) pairs of the upper diagonal of the NxN matrix
static Pair_Set exhaustivePairs(const size_t N)
{
  Pair_Set pairs;
  for(size_t I = 0; I < N; ++I)
    for(size_t J = I+1; J < N; ++J)
      pairs.insert(std::make_pair(I,J));

  return pairs;
}

/// Generate the pairs that have a distance inferior to the overlapSize
/// Usable to match video sequence
static Pair_Set contiguousWithOverlap(const size_t N, const size_t overlapSize)
{
  Pair_Set pairs;
  for(size_t I = 0; I < N; ++I)
    for(size_t J = I+1; J < I+1+overlapSize && J < N; ++J)
      pairs.insert(std::make_pair(I,J));
  return pairs;
}

/// Load a set of Pair_Set from a file
/// I J K L (pair that link I)
static bool loadPairs(
     const size_t N,  // number of image in the current project (to check index validity)
     const std::string &sFileName, // filename of the list file,
     Pair_Set & pairs)  // output pairs read from the list file
{
  std::ifstream in(sFileName.c_str());
  if(!in.is_open())
  {
    std::cerr << std::endl
      << "loadPairs: Impossible to read the specified file: \"" << sFileName << "\"." << std::endl;
    return false;
  }
  std::string sValue;
  std::vector<std::string> vec_str;
  while(std::getline( in, sValue ) )
  {
    vec_str.clear();
    stl::split(sValue, ' ', vec_str);
    const size_t str_size = vec_str.size();
    if (str_size < 2)
    {
      std::cerr << "loadPairs: Invalid input file: \"" << sFileName << "\"." << std::endl;
      return false;
    }
    std::stringstream oss;
    oss.clear(); oss.str(vec_str[0]);
    size_t I, J;
    oss >> I;
    for(size_t i=1; i<str_size ; ++i)
    {
      oss.clear(); oss.str(vec_str[i]);
      oss >> J;
      if( I > N-1 || J > N-1) //I&J always > 0 since we use unsigned type
      {
        std::cerr << "loadPairs: Invalid input file. Image out of range. "
                << "I: " << I << " J:" << J << " N:" << N << std::endl
                << "File: \"" << sFileName << "\"." << std::endl;
        return false;
      }
      if( I == J )
      {
        std::cerr << "loadPairs: Invalid input file. Image " << I << " see itself. File: \"" << sFileName << "\"." << std::endl;
        return false;
      }
      pairs.insert( (I < J) ? std::make_pair(I, J) : std::make_pair(J, I) );
    }
  }
  in.close();
  return true;
}

/// Save a set of Pair_Set to a file (one pair per line)
/// I J
/// I K
/// ...
static bool savePairs(const std::string &sFileName, const Pair_Set & pairs)
{
  std::ofstream outStream(sFileName.c_str());
  if(!outStream.is_open())  {
    std::cerr << std::endl
      << "savePairs: Impossible to open the output specified file: \"" << sFileName << "\"." << std::endl;
    return false;
  }
  for ( const auto & cur_pair : pairs ) 
    {
    outStream << cur_pair.first << ' ' << cur_pair.second << '\n';
  }
  bool bOk = !outStream.bad();
  outStream.close();
  return bOk;
}

}; // namespace openMVG
