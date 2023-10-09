// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_LANDMARK_HPP
#define OPENMVG_SFM_SFM_LANDMARK_HPP

#include "openMVG/types.hpp"

namespace openMVG {
namespace sfm {

/// Define 3D-2D tracking data: 3D landmark with its 2D observations

// This is an observation in one view (a point/feature id)j:w
// It is associated with a 3D feature in Landmark, one per view
struct Observation
{
  Observation():id_feat(UndefinedIndexT) {  }
  Observation(const Vec2 & p, IndexT idFeat): x(p), id_feat(idFeat) {}

  Vec2 x;
  IndexT id_feat;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const;

  // Serialization
  template <class Archive>
  void load( Archive & ar);
};

struct ObservationInfo {
  ObservationInfo() {  }
  ObservationInfo(const Vec2 & tgt): t(tgt) { }
  Vec2 t; // featue orientation / 2D unit tangent / optical flow / etc
};

/// Observations are indexed by their View_id
using Observations = Hash_Map<IndexT, Observation>;

/// ObservationsInfo are indexed by their View_id
/// they go parallel to Observations
using ObservationsInfo = Hash_Map<IndexT, ObservationInfo>;

/// Define a landmark (a 3D point, with its 2d observations)
struct Landmark
{
  Vec3 X;
  Observations obs;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const;

  template <class Archive>
  void load( Archive & ar);
};

// In the future, we can make a generic
// landmarksInfo, and make this a special case.
// But.. keep it simple, keep it maintainable..
struct LandmarkInfo
{
  Vec3 T; // 3D orientation corresponding to 2D feature/SIFT orientations

  double scale; // 3D scale
  /*
  // Serialization
  template <class Archive>
  void save( Archive & ar) const;

  template <class Archive>
  void load( Archive & ar);
  */
  ObservationsInfo obs_info;
};

/// Define a collection of landmarks are indexed by their TrackId
using Landmarks = Hash_Map<IndexT, Landmark>;
// Additional info that may be desired, parallel to Landmarks,
// that is, every addition to Landmark has to go in tandem with this one
using LandmarksInfo = Hash_Map<IndexT, LandmarkInfo>;

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_LANDMARK_HPP
