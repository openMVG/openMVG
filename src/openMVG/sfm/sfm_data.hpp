// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_HPP
#define OPENMVG_SFM_DATA_HPP

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"

namespace openMVG {

using namespace openMVG::geometry;

/// Define a collection of View
typedef Hash_Map<IndexT, View> Views;

/// Define a collection of Pose (indexed by View::id_pose)
typedef Hash_Map<IndexT, Pose3> Poses;

/// Define a collection of IntrinsicParameter (indexed by  View::id_intrinsic)
typedef Hash_Map<IndexT, std::shared_ptr<IntrinsicBase> > Intrinsics;

/// Define a collection of landmarks are indexed by their TrackId
typedef Hash_Map<IndexT, Landmark> Landmarks;

/// Generic SfM data container
/// Store structure and camera properties:
struct SfM_Data
{
  /// Considered views
  Views views;
  /// Considered poses (indexed by view.id_pose)
  Poses poses;
  /// Considered camera intrinsics (indexed by view.id_cam)
  Intrinsics intrinsics;

  /// Root Views path
  std::string s_root_path;

  //--
  // Accessors
  //--
  const Views & getViews() const {return views;}
  const Poses & getPoses() const {return poses;}
  const Intrinsics & getIntrinsics() const {return intrinsics;}
  const Landmarks & getLandmarks() const {return structure;}

  /// Structure (3D points with their 2D observations)
  Landmarks structure;
};

} // namespace openMVG

#endif // OPENMVG_SFM_DATA_HPP
