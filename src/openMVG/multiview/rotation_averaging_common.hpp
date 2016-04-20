
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_COMMON_H_
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_COMMON_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/types.hpp"
#include <vector>
#include <map>

namespace openMVG   {
namespace rotation_averaging  {

/// Representation of weighted relative rotations data between two poses
struct RelativeRotation 
{
  IndexT i, j; // pose's indices
  Mat3 Rij; // pose's relative rotation
  float weight;

  RelativeRotation(IndexT i_=0, IndexT j_=0, const	Mat3 & Rij_=Mat3::Identity(), float weight_=1.0f):
  i(i_), j(j_), Rij(Rij_), weight(weight_)
  {}
};

typedef std::vector<RelativeRotation> RelativeRotations;
typedef std::map<Pair, RelativeRotation> RelativeRotations_map;

/// List the pairs used by the relative rotations
inline Pair_Set getPairs(const RelativeRotations & relRots)
{
  Pair_Set pairs;
  for( const auto & cur_rotation : relRots ) 
    pairs.insert(std::make_pair(cur_rotation.i, cur_rotation.j));
  return pairs;
}

/// Convert a relative motion iterable sequence to RelativeRotation indexed by pairs
inline RelativeRotations_map getMap(const RelativeRotations & relRots)
{
  RelativeRotations_map map_rots;
  for( const auto & cur_rotation : relRots ) 
    map_rots[std::make_pair(cur_rotation.i, cur_rotation.j)] = cur_rotation ;
  return map_rots;
}

} // namespace rotation_averaging
} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_ROTATION_AVERAGING_COMMON_H_

