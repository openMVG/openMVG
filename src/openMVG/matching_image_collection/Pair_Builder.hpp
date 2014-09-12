// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdlib>
#include "openMVG/split/split.hpp"

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
/// Usable to match videos sequence
PairsT contiguousWithOverlap(const size_t N, const size_t overlapSize)
{
  PairsT pairs;
  for(int I = 0; I < N; ++I)
    for(int J = I+1; J < I+1+overlapSize && J < N; ++J)
    {
      pairs.insert(std::make_pair(I,J));
    }
  return pairs;
}

PairsT predefinedPairs(const std::string &sFileName)
{
  PairsT pairs;

  std::ifstream in(sFileName.c_str());
  if(!in.is_open())  {
    std::cerr << std::endl
      << "--pairList: Impossible to read the specified file." << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string sValue;
  std::vector<std::string> vec_str;
  while(getline( in, sValue ) )
  {
    vec_str.clear();
    split( sValue, " ", vec_str );
    size_t str_size=vec_str.size();
    if (str_size < 2)
    {
      std::cerr << "--pairList: Invalid input file" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::stringstream oss;
    oss.clear(); oss.str(vec_str[0]);
    size_t I, J;
    oss >> I;
    for(size_t i=1; i<str_size ; ++i)
    {
      oss.clear(); oss.str(vec_str[i]);
      oss >> J;
      pairs.insert(std::make_pair(I,J));
    }
  }
  in.close();
  return pairs;
}

}; // namespace openMVG
