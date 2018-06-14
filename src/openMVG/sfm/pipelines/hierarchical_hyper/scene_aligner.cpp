// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner_residuals.hpp"
#include "openMVG/types.hpp"

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres/types.h"
#include "ceres/cost_function.h"


namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

std::shared_ptr<cameras::IntrinsicBase> findBestIntrinsic(const SfM_Data& sfm_data_A, const SfM_Data& sfm_data_B, openMVG::IndexT cam_id)
{
  // both alternatives for the intrinsic
  const auto intrinsicA = sfm_data_A.intrinsics.at(cam_id);
  const auto intrinsicB = sfm_data_B.intrinsics.at(cam_id);

  // quick escape in case intrinsics are equal
  // TODO ? here we only check for pointer equality should we do more ?
  if (intrinsicA == intrinsicB)
    return std::shared_ptr<cameras::IntrinsicBase>(intrinsicA->clone());

  double RMSE_A(0.0), RMSE_B(0.0);
  int n_totalResiduals(0);
  for (const auto & sfm_data : {sfm_data_A, sfm_data_B})
  {
    for (const auto & landmark : sfm_data.GetLandmarks())
    {
      const Observations & observations = landmark.second.obs;
      for (const auto & obs: observations)
      {
        // we have to do the following check because observations of common landmarks are not pruned out when
        // clustering a submap in two...which makes it simpler to merge back together
        const auto & view = sfm_data.GetViews().find(obs.first);
        if (view == sfm_data.GetViews().end())
          continue;

        const IndexT & intrinsic_id = sfm_data.GetViews().at(obs.first)->id_intrinsic;
        if (intrinsic_id != cam_id)
          continue;

        const IndexT & pose_id = sfm_data.GetViews().at(obs.first)->id_pose;
        const auto & pose = sfm_data.GetPoses().at(pose_id);
        const Vec3 X = pose(landmark.second.X);
        const Vec2 residualA = intrinsicA->residual(X, obs.second.x);
        const Vec2 residualB = intrinsicB->residual(X, obs.second.x);
        RMSE_A += residualA(0) * residualA(0);
        RMSE_A += residualA(1) * residualA(1);
        RMSE_B += residualB(0) * residualB(0);
        RMSE_B += residualB(1) * residualB(1);
        ++n_totalResiduals;
      }
    }
  }

  RMSE_A = std::sqrt((RMSE_A)/(n_totalResiduals));
  RMSE_B = std::sqrt((RMSE_B)/(n_totalResiduals));
  std::cout << "RMSE_A : " << RMSE_A << " RMSE_B : "  << RMSE_B << std::endl;

  if (RMSE_A < RMSE_B)
    return std::shared_ptr<cameras::IntrinsicBase>(intrinsicA->clone());
  else
    return std::shared_ptr<cameras::IntrinsicBase>(intrinsicB->clone());
}

std::set<IndexT> getCommonCameraIds(const SfM_Data& sfm_data_1, const SfM_Data& sfm_data_2)
{
  const std::set<IndexT> cam_ids_1 = Get_Valid_Intrinsics_Ids(sfm_data_1);
  const std::set<IndexT> cam_ids_2 = Get_Valid_Intrinsics_Ids(sfm_data_2);

  std::set<IndexT> common_cam_ids;
  std::set_intersection(cam_ids_1.cbegin(), cam_ids_1.cend(),
                        cam_ids_2.cbegin(), cam_ids_2.cend(),
                        std::inserter(common_cam_ids, common_cam_ids.begin()));
  return common_cam_ids;
}

Landmarks getCommonReconstructedLandmarks(const SfM_Data& sfm_data_1, const SfM_Data& sfm_data_2, const std::set<IndexT>& track_ids)
{
  Landmarks common_landmarks;

  for (const auto & track_id : track_ids)
  {
    if (sfm_data_1.structure.count(track_id) > 0
        && sfm_data_2.structure.count(track_id) > 0)
    {
      // we copy from the first scene because we want the landmarks to be in the same referential
      Landmark separator = sfm_data_1.structure.at(track_id);

      // we don't want observations in the separator landmarks, we will only be interested
      // in the position of the landmark. later on observations are added to it via the copySfMDataSceneInto
      // function
      separator.obs.clear();

      common_landmarks[track_id] = separator;
    }
  }

  return common_landmarks;
}

bool MergeScenesUsingCommonTracks
(
  SfM_Data & destination_sfm_data,
  const SfM_Data & sfm_data_fst, // first submap scene
  const SfM_Data & sfm_data_snd, // second submap scene
  const std::set<IndexT> & separator_track_ids,
  SceneAligner *smap_aligner
)
{
  // 1. Find separator landmarks and transformation between two scenes
  Similarity3 sim;

  Landmarks separator_landmarks = getCommonReconstructedLandmarks(sfm_data_fst, sfm_data_snd, separator_track_ids);

  const bool alignment_successful =
      smap_aligner->computeTransformAndCommonLandmarks(
        separator_landmarks, sfm_data_fst, sfm_data_snd, sim);

  // 2. transform second scene and merge it with first scene
  if (alignment_successful)
  {
    destination_sfm_data.intrinsics.clear();
    destination_sfm_data.poses.clear();
    destination_sfm_data.structure.clear();

    // destination sfm data is in the same referential as first submap
    copySfMDataSceneInto(destination_sfm_data, sfm_data_fst, Similarity3());
    copySfMDataSceneInto(destination_sfm_data, sfm_data_snd, sim);

    // find best common intrinsics and update those in destination scene
    const auto common_cam_ids = getCommonCameraIds(sfm_data_fst, sfm_data_snd);
    for (const auto & cam_id : common_cam_ids)
    {
      destination_sfm_data.intrinsics[cam_id] = findBestIntrinsic(sfm_data_fst, sfm_data_snd, cam_id);
    }

    // update position of separator landmarks
    for (const auto & lmk : separator_landmarks)
    {
      // separator landmarks are in the referential of the first scene and the parent scene
      destination_sfm_data.structure.at(lmk.first).X = lmk.second.X;
    }

    return true;
  }
  else
  {
    return false;
  }
}

SceneAligner::SceneAligner(Bundle_Adjustment_Ceres::BA_Ceres_options options)
  : ceres_options_(options)
{}

bool SceneAligner::checkScenesAreAlignable(
    const SfM_Data &sfm_data_first,
    const SfM_Data &sfm_data_second,
    const Landmarks &common_landmarks)
{
  if (common_landmarks.empty())
    return false;

  for (const auto & mapped_lmk : common_landmarks)
  {
    const auto track_id = mapped_lmk.first;
    if (sfm_data_first.structure.count(track_id) == 0
        || sfm_data_second.structure.count(track_id) == 0)
    {
      return false;
    }
  }

  return true;
}

bool SceneAligner::computeTransformAndCommonLandmarks(Landmarks &destination_landmarks,
    const SfM_Data &sfm_data_first,
    const SfM_Data &sfm_data_second,
    Similarity3& sim)
{
  // check precondition for function to run normally (note : careful, this is a virtual method)
  if (!checkScenesAreAlignable(sfm_data_first, sfm_data_second, destination_landmarks))
  {
    return false;
  }

  ceres::Problem problem;

  const Mat3 R = sim.pose_.rotation();
  const Vec3 t = sim.pose_.center();
  double angleAxis[3];
  ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
  std::vector<double> second_base_node_pose =
    {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};
  double scaling_factor = sim.scale_;

  // note : configureProblem is a virtual method !
  configureProblem(
    problem,
    destination_landmarks,
    sfm_data_first,
    sfm_data_second,
    second_base_node_pose,
    scaling_factor
    );

  // Configure a BA engine and run it
  // Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.preconditioner_type =
    static_cast<ceres::PreconditionerType>(ceres_options_.preconditioner_type_);
  ceres_config_options.linear_solver_type =
    static_cast<ceres::LinearSolverType>(ceres_options_.linear_solver_type_);
  ceres_config_options.sparse_linear_algebra_library_type =
    static_cast<ceres::SparseLinearAlgebraLibraryType>(ceres_options_.sparse_linear_algebra_library_type_);
  ceres_config_options.minimizer_progress_to_stdout = false;
  ceres_config_options.logging_type = ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;
  ceres_config_options.max_num_iterations = 100;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  std::cout << summary.FullReport() << std::endl;

  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
    }

    // update similarity
    Vec3 second_base_node_t = Vec3(second_base_node_pose[3],second_base_node_pose[4],second_base_node_pose[5]);
    Mat3 second_base_node_RMat;
    ceres::AngleAxisToRotationMatrix(&second_base_node_pose[0], second_base_node_RMat.data());
    std::cout << " new aligned base node coords : \n"
              << second_base_node_t << std::endl
              << second_base_node_RMat << std::endl
              << " scale : " << scaling_factor << std::endl;

    sim.pose_ = {second_base_node_RMat.transpose(), second_base_node_t};
    sim.scale_ = scaling_factor;

    return true;
  }
}

void SceneAligner::configureProblem(ceres::Problem & problem,
    Landmarks &destination_landmarks,
    const SfM_Data & sfm_data_first,
    const SfM_Data & sfm_data_second,
    std::vector<double> & second_base_node_pose,
    double & scaling_factor
    )
{
  const Landmarks & landmarks_first = sfm_data_first.structure;
  const Landmarks & landmarks_second = sfm_data_second.structure;

  // Add scaling as a parameter
  problem.AddParameterBlock(&scaling_factor, 1);
  problem.SetParameterLowerBound(&scaling_factor, 0, 1e-7); // scaling factor should always be positive

  // Add second base node pose as a parameter
  problem.AddParameterBlock(&second_base_node_pose[0], 6);

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all landmarks that are common to both submaps, add 3d error / residuals
  for (auto & mapped_separator_landmark : destination_landmarks)
  {
    // measurements of the landmark in each submap
    const IndexT track_id = mapped_separator_landmark.first;
    const Landmark & first_measurement = landmarks_first.at(track_id);
    const Landmark & second_measurement = landmarks_second.at(track_id);
    Landmark & separator_landmark = mapped_separator_landmark.second;

    // Build the residual block corresponding to the landmark:
    ceres::CostFunction * cost_function =
      ResidualErrorFunctor_BaseNode_Separators::Create(first_measurement.X, second_measurement.X);

    if (cost_function)
    {
      problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &second_base_node_pose[0],
          separator_landmark.X.data(),// note that the 3d positions will be modified !
          &scaling_factor
          );
    }
  }
}

void copySfMDataSceneInto(
    SfM_Data & destination_sfm_data,
    const SfM_Data & original_sfm_data,
    const Similarity3& sim
    )
{
  SfM_Data transformed_sfm_data = original_sfm_data;
  ApplySimilarity(sim.inverse(), transformed_sfm_data);

  destination_sfm_data.views.insert(transformed_sfm_data.views.cbegin(), transformed_sfm_data.views.cend());
  destination_sfm_data.poses.insert(transformed_sfm_data.poses.cbegin(), transformed_sfm_data.poses.cend());
  destination_sfm_data.control_points.insert(transformed_sfm_data.control_points.cbegin(), transformed_sfm_data.control_points.cend());
  destination_sfm_data.s_root_path = original_sfm_data.s_root_path; // TODO : dangerous...should we check and find common directory + modify local view paths ?...

  // for the intrinsics we do a deep copy
  std::transform(transformed_sfm_data.intrinsics.cbegin(), transformed_sfm_data.intrinsics.cend(),
                 std::inserter(destination_sfm_data.intrinsics, destination_sfm_data.intrinsics.begin()),
                 [](const std::pair<IndexT, std::shared_ptr<cameras::IntrinsicBase>>& intrinsic_pair)
                    {auto deep_copy = std::shared_ptr<cameras::IntrinsicBase>(intrinsic_pair.second->clone());
                      return std::make_pair(intrinsic_pair.first, deep_copy);});

  // find common landmarks (we want to update observations on those)
  for (auto & lmk : destination_sfm_data.structure)
  {
    const auto lmk_id = lmk.first;
    const auto new_lmk = transformed_sfm_data.structure.find(lmk_id);
    if (new_lmk != transformed_sfm_data.structure.end())
    {
      lmk.second.obs.insert(new_lmk->second.obs.cbegin(), new_lmk->second.obs.cend());
    }
  }
  // the rest of the landmarks we insert directly
  destination_sfm_data.structure.insert(transformed_sfm_data.structure.cbegin(), transformed_sfm_data.structure.cend());
}

} // namespace sfm
} // namespace openMVG
