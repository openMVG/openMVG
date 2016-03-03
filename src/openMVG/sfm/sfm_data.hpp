// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_HPP
#define OPENMVG_SFM_DATA_HPP

#include "openMVG/types.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_metadata.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/cameras/cameras.hpp"

namespace openMVG {
namespace sfm {

/// Define a collection of View
typedef Hash_Map<IndexT, std::shared_ptr<View> > Views;

/// Define a collection of Pose (indexed by View::id_pose)
typedef Hash_Map<IndexT, geometry::Pose3> Poses;

/// Define a collection of IntrinsicParameter (indexed by View::id_intrinsic)
typedef Hash_Map<IndexT, std::shared_ptr<cameras::IntrinsicBase> > Intrinsics;

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
  /// Considered camera intrinsics (indexed by view.id_intrinsic)
  Intrinsics intrinsics;
  /// Structure (3D points with their 2D observations)
  Landmarks structure;
  /// Controls points (stored as Landmarks (id_feat has no meaning here))
  Landmarks control_points;

  /// Root Views path
  std::string s_root_path;

  //--
  // Accessors
  //--
  const Views & GetViews() const {return views;}
  const Poses & GetPoses() const {return poses;}
  const Intrinsics & GetIntrinsics() const {return intrinsics;}
  const Landmarks & GetLandmarks() const {return structure;}
  const Landmarks & GetControl_Points() const {return control_points;}

  std::set<IndexT> GetViewsKeys() const
  {
    std::set<IndexT> viewKeys;
    for(auto v: views)
        viewKeys.insert(v.first);
    return viewKeys;
  }

  /// Check if the View have defined intrinsic and pose
  bool IsPoseAndIntrinsicDefined(const View * view) const
  {
    if (view == NULL) return false;
    return (
      view->id_intrinsic != UndefinedIndexT &&
      view->id_pose != UndefinedIndexT &&
      intrinsics.find(view->id_intrinsic) != intrinsics.end() &&
      poses.find(view->id_pose) != poses.end());
  }

  /// Get the pose associated to a view
  const geometry::Pose3 GetPoseOrDie(const View * view) const
  {
    return poses.at(view->id_pose);
  }

  bool operator==(const SfM_Data& other) const {

    // Views
    if(views.size() != other.views.size())
      return false;
    for(Views::const_iterator it = views.begin(); it != views.end(); ++it)
    {
        const View& view1 = *(it->second.get());
        const View& view2 = *(other.views.at(it->first).get());
        if(!(view1 == view2))
          return false;
        
        // Image paths 
        if(s_root_path + view1.s_Img_path != other.s_root_path + view2.s_Img_path)
          return false;
    }

    // Poses
    if((poses != other.poses))
      return false;

    // Intrinsics
    if(intrinsics.size() != other.intrinsics.size())
      return false;

    Intrinsics::const_iterator it = intrinsics.begin();
    Intrinsics::const_iterator otherIt = other.intrinsics.begin();
    for(; it != intrinsics.end() && otherIt != other.intrinsics.end(); ++it, ++otherIt)
    {
        // Index
        if(it->first != otherIt->first)
          return false;

        // Intrinsic
        cameras::IntrinsicBase& intrinsic1 = *(it->second.get());
        cameras::IntrinsicBase& intrinsic2 = *(otherIt->second.get());
        if(!(intrinsic1 == intrinsic2))
          return false;
    }

    // Points IDs are not preserved
    if(structure.size() != other.structure.size())
      return false;
    Landmarks::const_iterator landMarkIt = structure.begin();
    Landmarks::const_iterator otherLandmarkIt = other.structure.begin();
    for(; landMarkIt != structure.end() && otherLandmarkIt != other.structure.end(); ++landMarkIt, ++otherLandmarkIt)
    {
        // Points IDs are not preserved
        // Landmark
        const Landmark& landmark1 = landMarkIt->second;
        const Landmark& landmark2 = otherLandmarkIt->second;
        if(!(landmark1 == landmark2))
          return false;
    }

    // Control points
    if(control_points != other.control_points)
      return false;

    // Root path can be reseted during exports

    return true;

  }
};

/**
 * @brief ColorizeTracks Add the associated color to each 3D point of
 * the sfm_data, using the track to determine the best view from which
 * to get the color.
 * @param sfm_data The container of the data
 */
void ColorizeTracks( SfM_Data & sfm_data );

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_HPP
