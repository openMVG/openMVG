// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/stl/split.hpp"
#include "openMVG/sfm/sfm_data.hpp"

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace openMVG {

/// Generate all the (I,J) pairs of the upper diagonal of the NxN matrix
static Pair_Set exhaustivePairs(const sfm::Views& views, int rangeStart=-1, int rangeSize=0)
{
  Pair_Set pairs;
  sfm::Views::const_iterator itA = views.begin();
  sfm::Views::const_iterator itAEnd = views.end();

  // If we have a rangeStart, only compute the matching for (rangeStart, X).
  if(rangeStart != -1 && rangeSize != 0)
  {
    if(rangeStart >= views.size())
      return pairs;
    std::advance(itA, rangeStart);
    itAEnd = views.begin();
    std::advance(itAEnd, std::min(std::size_t(rangeStart+rangeSize), views.size()));
  }
  
  for(; itA != itAEnd; ++itA)
  {
    sfm::Views::const_iterator itB = itA;
    std::advance(itB, 1);
    for(; itB != views.end(); ++itB)
      pairs.insert(std::make_pair(itA->first, itB->first));
  }
  return pairs;
}

/// Generate the pairs that have a distance inferior to the overlapSize
/// Usable to match video sequence
static Pair_Set contiguousWithOverlap(const sfm::Views& views, const size_t overlapSize, int rangeStart=-1, int rangeSize=0)
{
  Pair_Set pairs;
  sfm::Views::const_iterator itA = views.begin();
  sfm::Views::const_iterator itAEnd = views.end();

  // If we have a rangeStart, only compute the matching for (rangeStart, X).
  if(rangeStart != -1 && rangeSize != 0)
  {
    if(rangeStart >= views.size())
      return pairs;
    std::advance(itA, rangeStart);
    itAEnd = views.begin();
    std::advance(itAEnd, std::min(std::size_t(rangeStart+rangeSize), views.size()));
  }
  
  for(; itA != itAEnd; ++itA)
  {
    sfm::Views::const_iterator itB = itA;
    std::advance(itB, 1);
    sfm::Views::const_iterator itBEnd = itA;
    std::size_t advanceSize = overlapSize;
    if(std::distance(itBEnd, views.end()) < 1 + overlapSize)
    {
      advanceSize = std::distance(itBEnd, views.end()) - 1;
    }
    std::advance(itBEnd, 1 + advanceSize);
    
    for(; itB != views.end() && itB != itBEnd; ++itB)
      pairs.insert(std::make_pair(itA->first, itB->first));
  }
  return pairs;
}

/// Load a set of Pair_Set from a file
/// I J K L (pair that link I)
static bool loadPairs(
     const std::string &sFileName, // filename of the list file,
     Pair_Set & pairs,
     bool ordered=true, // output pairs read from the list file
     int rangeStart=-1,
     int rangeSize=0)
{
  std::ifstream in(sFileName.c_str());
  if(!in.is_open())
  {
    std::cerr << std::endl
      << "loadPairs: Impossible to read the specified file: \"" << sFileName << "\"." << std::endl;
    return false;
  }
  std::size_t nbLine = 0;
  std::string sValue;
  std::vector<std::string> vec_str;
  for(; std::getline( in, sValue ); ++nbLine)
  {
    if(rangeStart != -1 && rangeSize != 0)
    {
      if(nbLine < rangeStart)
        continue;
      if(nbLine >= rangeStart + rangeSize)
        break;
    }
    vec_str.clear();
    stl::split(sValue, " ", vec_str);
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
      if( I == J )
      {
        std::cerr << "loadPairs: Invalid input file. Image " << I << " see itself. File: \"" << sFileName << "\"." << std::endl;
        return false;
      }
      Pair pairToInsert;
      if(ordered)
      {
        // Insert pair with I = min(I, J) and J = max(I,J)
        pairToInsert = (( (I < J) ? std::make_pair(I, J) : std::make_pair(J, I) ));  
      }
      else
      {
        // Keep I & J
        Pair_Set::iterator it = pairs.find(std::make_pair(J, I));
        if(it == pairs.end())
          pairToInsert = std::make_pair(I, J);
      }
      pairs.insert(pairToInsert);
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
  for (Pair_Set::const_iterator iterP = pairs.begin();
    iterP != pairs.end(); ++iterP)
  {
    outStream << iterP->first << ' ' << iterP->second << '\n';
  }
  bool bOk = !outStream.bad();
  outStream.close();
  return bOk;
}

}; // namespace openMVG
