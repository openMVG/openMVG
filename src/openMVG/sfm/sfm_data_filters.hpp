
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_FILTERS_HPP
#define OPENMVG_SFM_DATA_FILTERS_HPP

#include "openMVG/types.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/sfm/sfm_data.hpp"

#include <iterator>
#include <set>

namespace openMVG {
namespace sfm {

/// List the view indexes that have valid camera intrinsic and pose.
static std::set<IndexT> Get_Valid_Views
(
  const SfM_Data & sfm_data
)
{
  std::set<IndexT> valid_idx;
  for (const auto & it : sfm_data.GetViews() ) 
  {
    const View * v = it.second.get();
    if (sfm_data.IsPoseAndIntrinsicDefined(v))
    {
      valid_idx.insert(v->id_view);
    }
  }
  return valid_idx;
}

/// Filter a list of pair: Keep only the pair that are defined in index list
template <typename IterablePairs, typename IterableIndex>
static Pair_Set Pair_filter
(
  const IterablePairs & pairs,
  const IterableIndex & index
)
{
  Pair_Set kept_pairs;
  for (auto& it : pairs)
  {
    if (index.count(it.first) > 0 &&
      index.count(it.second) > 0)
    kept_pairs.insert(it);
  }
  return kept_pairs;
}

// Remove tracks that have a small angle (tracks with tiny angle leads to instable 3D points)
// Return the number of removed tracks
static IndexT RemoveOutliers_PixelResidualError
(
  SfM_Data & sfm_data,
  const double dThresholdPixel,
  const unsigned int minTrackLength = 2
)
{
  IndexT outlier_count = 0;
  auto iterTracks = sfm_data.structure.begin();
  while (iterTracks != sfm_data.structure.end())
  {
    Observations & obs = iterTracks->second.obs;
    auto itObs = obs.begin();
    while (itObs != obs.end())
    {
      const View * view = sfm_data.views.at(itObs->first).get();
      const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
      const cameras::IntrinsicBase * intrinsic = sfm_data.intrinsics.at(view->id_intrinsic).get();
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      if (residual.norm() > dThresholdPixel)
      {
        ++outlier_count;
        itObs = obs.erase(itObs);
      }
      else
        ++itObs;
    }
    if (obs.empty() || obs.size() < minTrackLength)
      iterTracks = sfm_data.structure.erase(iterTracks);
    else
      ++iterTracks;
  }
  return outlier_count;
}

// Remove tracks that have a small angle (tracks with tiny angle leads to instable 3D points)
// Return the number of removed tracks
static IndexT RemoveOutliers_AngleError
(
  SfM_Data & sfm_data,
  const double dMinAcceptedAngle
)
{
  IndexT removedTrack_count = 0;
  Landmarks::iterator iterTracks = sfm_data.structure.begin();
  while (iterTracks != sfm_data.structure.end())
  {
    Observations & obs = iterTracks->second.obs;
    double max_angle = 0.0;
    for (Observations::const_iterator itObs1 = obs.begin();
      itObs1 != obs.end(); ++itObs1)
    {
      const View * view1 = sfm_data.views.at(itObs1->first).get();
      const geometry::Pose3 pose1 = sfm_data.GetPoseOrDie(view1);
      const cameras::IntrinsicBase * intrinsic1 = sfm_data.intrinsics.at(view1->id_intrinsic).get();

      auto itObs2 = itObs1;
      ++itObs2;
      for (; itObs2 != obs.end(); ++itObs2)
      {
        const View * view2 = sfm_data.views.at(itObs2->first).get();
        const geometry::Pose3 pose2 = sfm_data.GetPoseOrDie(view2);
        const cameras::IntrinsicBase * intrinsic2 = sfm_data.intrinsics.at(view2->id_intrinsic).get();

        const double angle = AngleBetweenRay(
          pose1, intrinsic1, pose2, intrinsic2,
          itObs1->second.x, itObs2->second.x);
        max_angle = std::max(angle, max_angle);
      }
    }
    if (max_angle < dMinAcceptedAngle)
    {
      iterTracks = sfm_data.structure.erase(iterTracks);
      ++removedTrack_count;
    }
    else
      ++iterTracks;
  }
  return removedTrack_count;
}

static bool eraseMissingPoses(SfM_Data & sfm_data, const IndexT min_points_per_pose)
{
  IndexT removed_elements = 0;
  const Landmarks & landmarks = sfm_data.structure;

  // Count the observation poses occurence
  Hash_Map<IndexT, IndexT> map_PoseId_Count;
  // Init with 0 count (in order to be able to remove non referenced elements)
  for (const auto & itPoses : sfm_data.GetPoses() )
  {
    map_PoseId_Count[itPoses.first] = 0;
  }

  // Count occurence of the poses in the Landmark observations
  for (const auto & itLandmarks : landmarks )
  {
    const Observations & obs = itLandmarks.second.obs;
    for ( const auto & itObs : obs ) 
    {
      const IndexT ViewId = itObs.first;
      const View * v = sfm_data.GetViews().at(ViewId).get();
      if (map_PoseId_Count.count(v->id_pose))
        map_PoseId_Count.at(v->id_pose) += 1;
      else
        map_PoseId_Count[v->id_pose] = 0;
    }
  }
  // If usage count is smaller than the threshold, remove the Pose
  for (Hash_Map<IndexT, IndexT>::const_iterator it = map_PoseId_Count.begin();
    it != map_PoseId_Count.end(); ++it)
  {
    if (it->second < min_points_per_pose)
    {
      sfm_data.poses.erase(it->first);
      ++removed_elements;
    }
  }
  return removed_elements > 0;
}

static bool eraseObservationsWithMissingPoses(SfM_Data & sfm_data, const IndexT min_points_per_landmark)
{
  IndexT removed_elements = 0;

  std::set<IndexT> pose_Index;
  std::transform(sfm_data.poses.begin(), sfm_data.poses.end(),
    std::inserter(pose_Index, pose_Index.begin()), stl::RetrieveKey());

  // For each landmark:
  //  - Check if we need to keep the observations & the track
  auto itLandmarks = sfm_data.structure.begin();
  while (itLandmarks != sfm_data.structure.end())
  {
    Observations & obs = itLandmarks->second.obs;
    auto itObs = obs.begin();
    while (itObs != obs.end())
    {
      const IndexT ViewId = itObs->first;
      const View * v = sfm_data.GetViews().at(ViewId).get();
      if (pose_Index.count(v->id_pose) == 0)
      {
        itObs = obs.erase(itObs);
        ++removed_elements;
      }
      else
        ++itObs;
    }
    if (obs.empty() || obs.size() < min_points_per_landmark)
      itLandmarks = sfm_data.structure.erase(itLandmarks);
    else
      ++itLandmarks;
  }
  return removed_elements > 0;
}

/// Remove unstable content from analysis of the sfm_data structure
static bool eraseUnstablePosesAndObservations(
  SfM_Data & sfm_data,
  const IndexT min_points_per_pose = 6,
  const IndexT min_points_per_landmark = 2)
{
  IndexT remove_iteration = 0;
  bool bRemovedContent = false;
  do
  {
    bRemovedContent = false;
    if (eraseMissingPoses(sfm_data, min_points_per_pose))
    {
      bRemovedContent = eraseObservationsWithMissingPoses(sfm_data, min_points_per_landmark);
      // Erase some observations can make some Poses index disappear so perform the process in a loop
    }
    remove_iteration += bRemovedContent ? 1 : 0;
  }
  while (bRemovedContent);

  return remove_iteration > 0;
}

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_FILTERS_HPP
