// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_FILTERS_FRUSTUM_HPP
#define OPENMVG_SFM_SFM_DATA_FILTERS_FRUSTUM_HPP

#include <string>
#include <utility>
#include <vector>

#include "openMVG/geometry/frustum.hpp"
#include "openMVG/types.hpp"

namespace openMVG { namespace sfm { struct SfM_Data; } }


namespace openMVG {
namespace sfm {

struct SfM_Data;

class Frustum_Filter
{
public:
  using FrustumsT = Hash_Map<IndexT, geometry::Frustum>;
  using NearFarPlanesT = Hash_Map<IndexT, std::pair<double, double>>;

  // Constructor
  Frustum_Filter
  (
    const SfM_Data & sfm_data,
    const double zNear = -1.,
    const double zFar = -1.,
    const NearFarPlanesT & z_near_z_far = NearFarPlanesT{}
  );

  // Init a frustum for each valid views of the SfM scene
  void initFrustum(const SfM_Data & sfm_data);

  // Return intersecting View frustum pairs. An optional bounding volume
  // defined as a vector of half-plane objects can also be provided to further
  // limit the intersection area.
  Pair_Set getFrustumIntersectionPairs(
    const std::vector<geometry::halfPlane::HalfPlaneObject>& bounding_volume = {}
  ) const;

  // Export defined frustum in PLY file for viewing
  bool export_Ply(const std::string & filename) const;

private:

  /// Init near and far plane depth from SfM_Data structure or defined value
  void init_z_near_z_far_depth
  (
    const SfM_Data & sfm_data,
    const double zNear = -1.,
    const double zFar = -1.
  );

  //--
  // Data
  //--
  bool _bTruncated; // Tell if we use truncated or infinite frustum

  FrustumsT frustum_perView; // Frustum for the valid view (defined pose+intrinsic)

  NearFarPlanesT z_near_z_far_perView; // Near & Far plane distance per view
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_FILTERS_FRUSTUM_HPP
