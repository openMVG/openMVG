// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/stellar/stellar_solver.hpp"

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"

#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/stellar/stellar_definitions.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"

#include "openMVG/tracks/tracks.hpp"
#include "openMVG/types.hpp"

#include <algorithm>

#include "ceres/ceres.h"

namespace openMVG {
namespace sfm {

/// Function that estimate the relative scale of a triplet of pose.
/// Since a triplet is made of two edges, we list the 3 view tracks
/// Then compute the median depth for each pair and save the depth ratio
bool EstimateTripletRelativeScale
(
  const Pair_Set & pairs,  // 2 pair of Pose Ids
  const SfM_Data & sfm_data,
  const Features_Provider * features_provider,
  const Matches_Provider * matches_provider,
  const Hash_Map<Pair, geometry::Pose3> & relative_poses,
  Relative_Scale & relative_scale
)
{
  if (pairs.size() != 2)
  {
    return false;
  }
  // Assert that:
  // - at least there is 3 poses => Required condition to compute a scaling factor
  // - that we have valid relative pose for triplet edges.
  Hash_Map<IndexT, IndexT> pose_count;
  std::set<IndexT> set_pose;
  for (const Pair & it : pairs)
  {
    if (relative_poses.count(it) == 0)
      return false;

    ++pose_count[it.first];
    ++pose_count[it.second];
    set_pose.insert(it.first);
    set_pose.insert(it.second);
  }
  if (pose_count.size() != 3)
  {
    return false;
  }
  // Identify the common pose node_id
  const IndexT node_id = [&]
  {
    for (const auto id_count : pose_count)
    {
      if (id_count.second == 2)
        return id_count.first;
    }
    return UndefinedIndexT;
  }();

  if (node_id == UndefinedIndexT)
    return false;

  //
  // Compute tracks/landmarks visibility for the 3-uplet
  //
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  {
    matching::PairWiseMatches matches;
    //-- List all view that shared some content with the used poses
    for (const auto & matches_it : matches_provider->pairWise_matches_)
    {
      const Pair & pair = matches_it.first;

      const View
        * view_I = sfm_data.GetViews().find(pair.first)->second.get(),
        * view_J = sfm_data.GetViews().find(pair.second)->second.get();

      if (set_pose.count(view_I->id_pose)
          && set_pose.count(view_J->id_pose))
      {
        const matching::IndMatches & pair_matches = matches_it.second;
        matches[pair] = pair_matches;
      }
    }

    // Computing tracks
    openMVG::tracks::TracksBuilder tracksBuilder;
    tracksBuilder.Build(matches);
    tracksBuilder.Filter(3);

    // If there is insufficient 3-view track count,
    // We reject this triplet, since the relative scale factor cannot be reliably computed.
    if (tracksBuilder.NbTracks() < 15)
    {
      return false;
    }

    tracksBuilder.ExportToSTL(map_tracksCommon);
  }

  //
  // Triangulate observations for each pose pair
  // Store track depth per pose pair
  // Compute relative the median depth scale ratio factor
  //
  std::map<Pair, std::vector<double>> depths;
  for (const auto & tracks : map_tracksCommon)
  {
    const tracks::submapTrack & track = tracks.second;

    for (const Pair & cu_pair : pairs)
    {
      //
      // Triangulate the observed tracks
      // and store the depth per pair to compute the scale ratio between the pair
      //

      // Check the track is supported at least by 3 poses
      {
        std::set<IndexT> poses_id;
        for (const auto & track_it : track)
        {
          const IndexT view_idx = track_it.first;
          const View * view = sfm_data.GetViews().find(view_idx)->second.get();
          poses_id.insert(view->id_pose);
        }
        if (poses_id.size() < 3)
          continue;
      }

      std::map<IndexT, geometry::Pose3> poses;
      std::vector<Vec3> bearing;
      std::vector<Mat34> vec_poses;

      for (const auto & track_it : track)
      {
        const IndexT view_idx = track_it.first;
        const View * view = sfm_data.GetViews().find(view_idx)->second.get();

        if (view->id_pose == cu_pair.first || view->id_pose == cu_pair.second)
        {
          const geometry::Pose3 pose = (view->id_pose == cu_pair.first) ?
            geometry::Pose3(Mat3::Identity(), Vec3::Zero())
            : relative_poses.at(cu_pair);

          const auto cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();

          poses[view->id_pose] = pose;

          vec_poses.emplace_back(pose.asMatrix());
          const size_t feat_idx = track_it.second;
          const Vec2 feat_pos = features_provider->feats_per_view.at(view_idx)[feat_idx].coords().cast<double>();
          bearing.emplace_back((*cam)(cam->get_ud_pixel(feat_pos)));
        }
      }
      const Eigen::Map<const Mat3X> bearing_matrix(bearing[0].data(), 3, bearing.size());
      Vec4 Xhomogeneous;
      if (TriangulateNViewAlgebraic(bearing_matrix, vec_poses, &Xhomogeneous))
      {
        const Vec3 X = Xhomogeneous.hnormalized();

        const geometry::Pose3 pose = poses[node_id];
        const double depth = (X - pose.center()).norm();
        // Store the depth for this 2-uplet of edges
        depths[cu_pair].push_back(depth);
      }
    }
  }

  if (depths.size() != 2) // Else one of the view did not produce any depth data
  {
    return false;
  }

  // Compute median depths
  std::vector<float> median_depths;
  median_depths.reserve(2);
  for (auto & depth_it : depths)
  {
    if (depth_it.second.empty())
      continue;

    std::nth_element(
      depth_it.second.begin(),
      depth_it.second.begin()+(depth_it.second.size()/2),
      depth_it.second.end());
    const double median = depth_it.second[depth_it.second.size()/2];
    median_depths.push_back(median);
    //std::cout << "(" << depth_it.first.first << "," << depth_it.first.second << ") : " << median << ", "
    //  << depth_it.second.size() << " points." << std::endl;
  }

  const double depth_ratio = median_depths[0] / median_depths[1];

  // Return the computed triplet depth ratio:
  relative_scale = {*pairs.cbegin(), *std::next(pairs.cbegin(), 1), depth_ratio};

  const bool b_refine_triplet = false;
  if (b_refine_triplet)
  {
    // Build a triplet that contained the computed depth ratio:
    const Pair pair01 = relative_scale.pairs[0];
    const Pair pair12 = relative_scale.pairs[1];

    Hash_Map<IndexT, geometry::Pose3> triplet_pose;

    triplet_pose[node_id] = geometry::Pose3(); // Identity
    if (pair01.first == node_id)
    {
      geometry::Pose3 relative_pose = relative_poses.at(pair01);
      relative_pose.center() /= depth_ratio;
      triplet_pose[pair01.second] = relative_pose;
    }
    else
    {
      geometry::Pose3 relative_pose = relative_poses.at(pair01).inverse();
      relative_pose.center() /= depth_ratio;
      triplet_pose[pair01.first] = relative_pose;
    }

    if (pair12.first == node_id)
    {
      geometry::Pose3 relative_pose = relative_poses.at(pair12);
      triplet_pose[pair12.second] = relative_pose;
    }
    else
    {
      geometry::Pose3 relative_pose = relative_poses.at(pair12).inverse();
      triplet_pose[pair12.first] = relative_pose;
    }

    // Create a scene containing the triplet and save it to disk:
    SfM_Data tiny_scene;
    for (const auto &  triplet_pose_it : triplet_pose)
    {
      const IndexT pose_id = triplet_pose_it.first;
      // Add view
      // Add intrinsic
      // Add poses
      const View * view = sfm_data.GetViews().at(pose_id).get();
      tiny_scene.views.insert(*sfm_data.GetViews().find(pose_id));
      tiny_scene.intrinsics.insert(*sfm_data.GetIntrinsics().find(view->id_intrinsic));
      tiny_scene.poses[pose_id] = triplet_pose[pose_id];
    }

    // Add tracks
    Landmarks & landmarks = tiny_scene.structure;
    for (const auto & tracks : map_tracksCommon)
    {
      const tracks::submapTrack & track = tracks.second;
      Observations obs;
      for (const auto & track_it : track)
      {
        const IndexT I = track_it.first;
        const size_t featIndexI = track_it.second;
        const Vec2 x = features_provider->feats_per_view.at(I)[featIndexI].coords().cast<double>();

        obs[I] = std::move(Observation(x, featIndexI));
      }
      landmarks[tracks.first].obs = std::move(obs);
    }

    // Triangulation
    sfm::SfM_Data_Structure_Computation_Blind structure_estimator(false);
    structure_estimator.triangulate(tiny_scene);

    // - refine only Structure and Rotations & translations (keep intrinsic constant)
    Bundle_Adjustment_Ceres::BA_Ceres_options options(true, false);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    const Optimize_Options ba_refine_paramerer_options
    (
      cameras::Intrinsic_Parameter_Type::NONE,
      Extrinsic_Parameter_Type::ADJUST_ALL,
      Structure_Parameter_Type::ADJUST_ALL
    );
    if (!bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_paramerer_options))
    {
      return false;
    }
  } // end -- b_refine_triplet

  return true;
}

Stellar_Solver::Stellar_Solver
(
  const StellarPod & stellar_pod,
  const Hash_Map<Pair, Pose3> & relative_poses,
  const SfM_Data & sfm_data,
  const Matches_Provider * matches_provider,
  const Features_Provider * features_provider,
  const bool use_all_matches,
  const bool use_threading
):stellar_pod_(stellar_pod),
  relative_poses_(relative_poses),
  sfm_data_(sfm_data),
  matches_provider_(matches_provider),
  features_provider_(features_provider),
  use_all_matches_(use_all_matches),
  use_threading_(use_threading)
{
}

bool Stellar_Solver::Solve(Poses & poses)
{
  // Check if the stellar pod is solvable
  // - one common node defined in every pair
  const std::set<IndexT> set_pose = [&] {
    std::set<IndexT> pose_set;
    for (const Pair & it : stellar_pod_)
    {
      pose_set.insert(it.first);
      pose_set.insert(it.second);
    }
    return pose_set;
  }();
  const std::vector<IndexT> pose_count = [&]{
    std::vector<IndexT> vec(set_pose.size(), 0);
    for (const Pair & it : stellar_pod_)
    {
      ++vec[std::distance(set_pose.cbegin(), set_pose.find(it.first))];
      ++vec[std::distance(set_pose.cbegin(), set_pose.find(it.second))];
    }
    return vec;
  }();
  if (pose_count.size() < 3)
    return false;

  const auto iterator_node_id = max_element(pose_count.cbegin(), pose_count.cend());

  if (iterator_node_id == pose_count.cend())
    return false;

  const IndexT central_node_id = *std::next(set_pose.cbegin(), std::distance(pose_count.cbegin(), iterator_node_id));

  const std::vector<Pair_Set> edge_two_uplets = ListEdge2Uplets();

  std::cout
    << "Stellar pod details:\n"
    << "#central pose id: " << central_node_id << "\n"
    << "#pairs: " << stellar_pod_.size() << "\n"
    << "#2-uplets: " << edge_two_uplets.size() << std::endl;

  std::vector<Relative_Scale> relative_scales;
  if (!Solve2UpletsRelativeScales(edge_two_uplets, relative_scales))
    return false;

  Hash_Map<IndexT, geometry::Pose3> stellar_poses;
  if (!SolveStellarPoses(central_node_id, relative_scales, stellar_poses))
    return false;

  SfM_Data sfm_data_stellar_pod;
  if (!Optimize(sfm_data_stellar_pod, stellar_poses, relative_scales))
    return false;
  poses = sfm_data_stellar_pod.poses;

  /*std::ostringstream os;
  os << "Stellar_" << central_node_id;
  Save(sfm_data_stellar_pod,
       stlplus::create_filespec("./",os.str(),".ply"),
       ESfM_Data(ALL));
       */
  return poses.size() >= 2;
}

std::vector<Pair_Set> Stellar_Solver::ListEdge2Uplets()
{
  // List the possible edge 2-uplets (triplet of poses)
  std::vector<Pair_Set> two_uplets;
  for (auto it0 = stellar_pod_.cbegin(); it0 != stellar_pod_.cend(); ++it0)
  {
    for (auto it1 = std::next(it0, 1); it1 != stellar_pod_.cend(); ++it1)
    {
      two_uplets.push_back({*it0, *it1});
    }
  }
  return two_uplets;
}

bool Stellar_Solver::Solve2UpletsRelativeScales
(
  const std::vector<Pair_Set> & edge_two_uplets,
  std::vector<Relative_Scale> & relative_scales
)
{
  // List possible 2-uplet and solve their relative scales
  // some 2-uplet cannot lead to a relative scale if they don't share a sufficient track amount
  for (const auto edge_two_uplet : edge_two_uplets)
  {
    Relative_Scale relative_scale;
    if (EstimateTripletRelativeScale(
          edge_two_uplet,  // Pair of pose Ids that define a triplet of pose
          sfm_data_,
          features_provider_,
          matches_provider_,
          relative_poses_,
          relative_scale))
    {
      relative_scales.emplace_back(relative_scale);
    }
  }
  return !relative_scales.empty();
}

bool Stellar_Solver::SolveStellarPoses
(
  const IndexT & central_node_id,
  const std::vector<Relative_Scale> & relative_scales,
  Hash_Map<IndexT, geometry::Pose3> & triplet_poses
)
{
  // Re-scale edges to a common scale
  if (!Solve_stellar_translation_scales_averaging(
    central_node_id,
    relative_scales,
    relative_poses_,
    triplet_poses,
    Stellar_Translation_Averaging_Solver_Type::SCALING_SOLVER_L2_FULL))
  {
    return false;
  }
  return true;
}

bool Stellar_Solver::Optimize
(
  SfM_Data & stellar_pod_reconstruction,
  const Hash_Map<IndexT, geometry::Pose3> & triplet_poses,
  const std::vector<Relative_Scale> & relative_scales
)
{
  if (triplet_poses.size() < 3)
  {
    return false;
  }

  stellar_pod_reconstruction.poses.clear();
  stellar_pod_reconstruction.intrinsics.clear();
  stellar_pod_reconstruction.structure.clear();
  stellar_pod_reconstruction.views.clear();

  // 1. Fill the sfm_data scene with view, intrinsic, poses
  // 2. Select valid matches pairs, build track and triangulate them
  // 3. Refine the sfm_data scene

  // Add the poses
  {
    for (const auto & pose_it : triplet_poses)
    {
      stellar_pod_reconstruction.poses[pose_it.first] = pose_it.second;
    }
  }

  // Collect the pairs used by this stellar pod
  const Pair_Set & used_pairs = Relative_Scale::Get_pairs(relative_scales);

  // Collect, matches, intrinsics and views data linked to the poses ids
  openMVG::tracks::STLMAPTracks tracks;
  {
    matching::PairWiseMatches matches;
    for (const auto & matches_it : matches_provider_->pairWise_matches_)
    {
      const Pair & pair = matches_it.first;

      const View * view_I = sfm_data_.GetViews().find(pair.first)->second.get();
      const View * view_J = sfm_data_.GetViews().find(pair.second)->second.get();

      if (triplet_poses.count(view_I->id_pose)
          && triplet_poses.count(view_J->id_pose))
      {
        // Collect matches
        const Pair pose_pair(std::min(view_I->id_pose,view_J->id_pose),
                             std::max(view_I->id_pose,view_J->id_pose));
        const matching::IndMatches & pair_matches = matches_it.second;
        if (!use_all_matches_)
        {
          if (used_pairs.count(pose_pair))
            matches[pair] = pair_matches;
          else
            continue; // This pair is ignored
        }
        else
        {
          matches[pair] = pair_matches;
        }

        // Add view information and related intrinsic data
        const std::array<const View*, 2> view_ids = {view_I, view_J};
        for (const View* view : view_ids)
        {
          if (stellar_pod_reconstruction.views.count(view->id_view) == 0)
            stellar_pod_reconstruction.views.insert(*(sfm_data_.GetViews().find(view->id_view)));
          if (stellar_pod_reconstruction.intrinsics.count(view->id_intrinsic) == 0)
            stellar_pod_reconstruction.intrinsics.insert(*(sfm_data_.GetIntrinsics().find(view->id_intrinsic)));
        }
      }
    }
    // Computing tracks
    {
      openMVG::tracks::TracksBuilder tracksBuilder;
      tracksBuilder.Build(matches);
      tracksBuilder.Filter(2); // [2-n] view based matches
      tracksBuilder.ExportToSTL(tracks);
    }
  }

  // Fill the sfm_data landmark observations and 3d initial position (triangulation):
  {
    // Add the track observations to the sfm_data
    Landmarks & landmarks = stellar_pod_reconstruction.structure;
    for (const auto & tracks_it : tracks)
    {
      const tracks::submapTrack & track = tracks_it.second;

      // Check if the track is observed by more than 2 differnt pose ids
      {
        std::set<IndexT> poses_id;
        for (const auto & track_it : track)
        {
          const IndexT view_idx = track_it.first;
          const View * view = sfm_data_.GetViews().find(view_idx)->second.get();
          poses_id.insert(view->id_pose);
          if (poses_id.size() > 2) break; // early exit
        }
        if (poses_id.size() < 2)
          continue;
      }
      // Collect the views and features observing this landmark
      landmarks[tracks_it.first].obs = [&]
      {
        Observations obs;
        for (const auto & track_it : track)
        {
          const IndexT view_idx = track_it.first;
          const IndexT feat_idx = track_it.second;
          const Vec2 x = features_provider_->feats_per_view.at(view_idx)[feat_idx].coords().cast<double>();

          obs[view_idx] = {x, feat_idx};
        }
        return obs;
      }();
    }
    // Triangulate
    sfm::SfM_Data_Structure_Computation_Blind structure_estimator(false);
    structure_estimator.triangulate(stellar_pod_reconstruction);
  }

  // - refine only Structure and Rotations & translations (keep intrinsic constant)
  Bundle_Adjustment_Ceres::BA_Ceres_options options(false, use_threading_);
  options.linear_solver_type_ = ceres::DENSE_SCHUR;
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  const Optimize_Options ba_refine_parameter_options
  (
    cameras::Intrinsic_Parameter_Type::NONE,
    Extrinsic_Parameter_Type::ADJUST_ALL,
    Structure_Parameter_Type::ADJUST_ALL
  );
  if (!bundle_adjustment_obj.Adjust(stellar_pod_reconstruction, ba_refine_parameter_options))
  {
    return false;
  }

  // In order to ensure that the scene is refined correctly and valid we perform
  // a second run of optimization where the matches that most lead to invalid triangulation are removed.
  {
    const IndexT min_point_per_pose = 12; // 6 min -> Keep more since depending of the camera intrinsic DoF we can need more
    const IndexT min_track_length = 3; // 2 min
    const double depth_median_limit_factor =  5.2;    //5.2 * median ~= X84,
    DepthCleaning(
      stellar_pod_reconstruction,
      depth_median_limit_factor,
      min_point_per_pose,
      min_track_length);

    // Remove outliers (max_angle, residual error)
    const size_t pointcount_initial = stellar_pod_reconstruction.structure.size();
    RemoveOutliers_PixelResidualError(stellar_pod_reconstruction, 4.0);
    RemoveOutliers_AngleError(stellar_pod_reconstruction, 2.0);

    // Remove poses that does not cover a sufficient number of observations (some observations are removed too)
    eraseUnstablePosesAndObservations(stellar_pod_reconstruction, min_point_per_pose, min_track_length);

    if (!bundle_adjustment_obj.Adjust(stellar_pod_reconstruction, ba_refine_parameter_options))
    {
      return false;
    }
  }

  return (stellar_pod_reconstruction.GetPoses().size() >= 3);
}

} // namespace sfm
} // namespace openMVG
