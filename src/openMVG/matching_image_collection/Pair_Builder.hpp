// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/split/split.hpp"

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

namespace openMVG {

typedef std::set<std::pair<size_t, size_t> > PairsT;

/// Generate all the (I,J) pairs of the upper diagonal of the NxN matrix
PairsT exhaustivePairs(const size_t N)
{
  PairsT pairs;
  for(size_t I = 0; I < N; ++I)
    for(size_t J = I+1; J < N; ++J)
      pairs.insert(std::make_pair(I,J));

  return pairs;
}

/// Generate the pairs that have a distance inferior to the overlapSize
/// Usable to match video sequence
PairsT contiguousWithOverlap(const size_t N, const size_t overlapSize)
{
  PairsT pairs;
  for(size_t I = 0; I < N; ++I)
    for(size_t J = I+1; J < I+1+overlapSize && J < N; ++J)
      pairs.insert(std::make_pair(I,J));
  return pairs;
}

/// Load a set of Pairs from a file
/// I J K L (pair that link I)
bool loadPairs(
     const size_t N,  // number of image in the current project (to check index validity)
     const std::string &sFileName, // filename of the list file,
     PairsT & pairs)  // output pairs read from the list file
{
  std::ifstream in(sFileName.c_str());
  if(!in.is_open())  {
    std::cerr << std::endl
      << "--loadPairs: Impossible to read the specified file." << std::endl;
    return false;
  }
  std::string sValue;
  std::vector<std::string> vec_str;
  while(std::getline( in, sValue ) )
  {
    vec_str.clear();
    split( sValue, " ", vec_str );
    const size_t str_size = vec_str.size();
    if (str_size < 2)
    {
      std::cerr << "--loadPairs: Invalid input file" << std::endl;
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
      if( I < 0 || I > N-1 || J < 0 || J > N-1)
      {
        std::cerr << "--loadPairs: Invalid input file. Image out of range" << std::endl;
        return false;
      }
      else
      {
        if( I != J )
        {
          pairs.insert( (I < J) ? std::make_pair(I, J) : std::make_pair(J, I) );
        }
        else
        {
          std::cerr << "--loadPairs: Invalid input file. Image see herself " << std::endl;
          return false;
        }
      }
    }
  }
  in.close();
  return true;
}

/// Save a set of Pairs to a file (one pair per line)
/// I J
/// I K
/// ...
bool savePairs(const std::string &sFileName, const PairsT & pairs)
{
  std::ofstream outStream(sFileName.c_str());
  if(!outStream.is_open())  {
    std::cerr << std::endl
      << "--savePairs: Impossible to open the output specified file." << std::endl;
    return false;
  }
  for (PairsT::const_iterator iterP = pairs.begin();
    iterP != pairs.end(); ++iterP)
  {
    outStream << iterP->first << ' ' << iterP->second << '\n';
  }
  bool bOk = !outStream.bad();
  outStream.close();
  return bOk;
}

}; // namespace openMVG
