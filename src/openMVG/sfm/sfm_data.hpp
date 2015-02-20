// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_HPP
#define OPENMVG_SFM_DATA_HPP

#include "openMVG/geometry/pose3.hpp"
#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include <string>

namespace openMVG {

using namespace openMVG::geometry;

/// A view define an image by a string and unique indexes for the view, the camera intrinsic & the pose
struct View
{
  std::string s_Img_path;

  // Id of the view
  IndexT id_view;

  // Index of intrinsics and the pose
  IndexT id_cam, id_pose;
};

typedef Hash_Map<IndexT, View> Views;
typedef Hash_Map<IndexT, Pose3> Poses;
typedef Hash_Map<IndexT, Intrinsic*> Intrinsics;

/// Define 3D-2D tracking data: 3D landmark with their 2D observations
struct Observation
{
  Observation() {}
  Observation(const Vec2 & p, IndexT idFeat): pos(p), id_feat(idFeat) {}

  Vec2 pos;
  IndexT id_feat;
};
/// Observations are indexed by their View_id
typedef Hash_Map<IndexT, Observation> Observations;

/// Define a landmark (a 3D point, with it's 2d observations)
struct Landmark
{
  Observations obs;
  Vec3 pt3;
};
/// Landmarks are indexed by their TrackId
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

  //--
  // Accessors
  //--
  const Views & getViews() const {return views;}
  const Poses & getPoses() const {return poses;}
  const Intrinsics & getIntrinsics() const {return intrinsics;}

  /// Structure (3D points with their 2D observations)
  Landmarks structure;
};

} // namespace openMVG

#endif // OPENMVG_SFM_DATA_HPP
