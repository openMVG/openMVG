// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner_residuals.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/types.hpp"

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres/types.h"
#include "ceres/cost_function.h"


namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

bool MergeScenesUsingCommonTracks
(
  SfM_Data & destination_sfm_data,
  const SfM_Data & sfm_data_first, // first submap scene
  const SfM_Data & sfm_data_second, // second submap scene
  const std::set<IndexT> & common_track_ids,
  SceneAligner *smap_aligner
)
{
  // 1. Find separator landmarks and transformation between two scenes
  Similarity3 sim;

  Landmarks separator_landmarks;
  std::copy_if(sfm_data_first.structure.cbegin(), sfm_data_first.structure.cend(),
               std::inserter(separator_landmarks, separator_landmarks.begin()),
               [&common_track_ids](const std::pair<IndexT, Landmark>& landmark_pair)
                  {return common_track_ids.count(landmark_pair.first) > 0;});

  const bool alignment_successful =
      smap_aligner->computeTransformAndCommonLandmarks(
        separator_landmarks, sfm_data_first, sfm_data_second,
        sim, common_track_ids);

  // 2. transform second scene and merge it with first scene
  if (alignment_successful)
  {
    destination_sfm_data.intrinsics.clear();
    destination_sfm_data.poses.clear();
    destination_sfm_data.structure.clear();

    // initialize destination sfm data
    destination_sfm_data.intrinsics = sfm_data_first.intrinsics;
    for (const auto & intrinsic : sfm_data_second.GetIntrinsics())
    {
      if (destination_sfm_data.intrinsics.find(intrinsic.first) == destination_sfm_data.intrinsics.end())
        destination_sfm_data.intrinsics[intrinsic.first] = intrinsic.second;
    }
    destination_sfm_data.views = sfm_data_first.views; // TODO : should we really deal with views here ?
    destination_sfm_data.poses = sfm_data_first.poses;
    destination_sfm_data.structure = sfm_data_first.structure;
    // overwrite separators with newly computed data
    destination_sfm_data.structure.insert(separator_landmarks.cbegin(), separator_landmarks.cend());

    // update the landmarks, poses and views from the second
    // scene transformed into the destination scene.
    transformSfMDataSceneInto(destination_sfm_data, sfm_data_second, sim);
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
    const std::set<IndexT> &common_track_ids)
{
  if (common_track_ids.empty())
    return false;

  for (const auto & track_id : common_track_ids)
  {
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
    Similarity3& sim,
    const std::set<IndexT> & common_track_ids)
{
  // check precondition for function to run normally
  if (!checkScenesAreAlignable(sfm_data_first, sfm_data_second, common_track_ids))
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
    scaling_factor,
    common_track_ids
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
    double & scaling_factor,
    const std::set<IndexT> & common_track_ids
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
  for (const auto & track_id : common_track_ids)
  {
    // measurements of the landmark in each submap
    const Landmark & first_measurement = landmarks_first.at(track_id);
    const Landmark & second_measurement = landmarks_second.at(track_id);
    Landmark & separator_landmark = destination_landmarks.at(track_id);

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

/**
 * @brief takes a SfM_Data scene and transforms it into a destination scene, with a geometric transformation.
 * @note mostly used for merging submaps together
 * @warning landmarks, poses and views from original scene will overwrite conflicting data
 * (i.e. map elements with same keys) in the destination scene ! Use with caution
 * @param destination_sfm_data : the destination scene, can already contain views, landmarks ... etc
 * @param original_sfm_data : the scene to be transformed into the destination scene
 * @param rotation
 * @param translation
 * @param scaling_factor
 */
void transformSfMDataSceneInto(
    SfM_Data & destination_sfm_data,
    const SfM_Data & original_sfm_data,
    const Similarity3& sim
    )
{
  SfM_Data transformed_sfm_data = original_sfm_data;
  ApplySimilarity(sim.inverse(), transformed_sfm_data);

  destination_sfm_data.views.insert(transformed_sfm_data.views.cbegin(), transformed_sfm_data.views.cend());
  destination_sfm_data.poses.insert(transformed_sfm_data.poses.cbegin(), transformed_sfm_data.poses.cend());
  destination_sfm_data.structure.insert(transformed_sfm_data.structure.cbegin(), transformed_sfm_data.structure.cend());
  destination_sfm_data.control_points.insert(transformed_sfm_data.control_points.cbegin(), transformed_sfm_data.control_points.cend());
}

} // namespace sfm
} // namespace openMVG
