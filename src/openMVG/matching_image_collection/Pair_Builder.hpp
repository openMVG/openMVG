// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

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

}; // namespace openMVG
