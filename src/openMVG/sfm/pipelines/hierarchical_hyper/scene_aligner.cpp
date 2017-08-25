// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner.hpp"
#include "openMVG/sfm/pipelines/hierarchical_hyper/scene_aligner_residuals.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres/types.h"
#include "ceres/cost_function.h"


namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

bool scenesAreAlignable(const SfM_Data & sfm_data_first, const SfM_Data & sfm_data_second, const std::vector<size_t> & common_track_ids)
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

SceneAligner::SceneAligner(Bundle_Adjustment_Ceres::BA_Ceres_options options)
  : ceres_options_(options)
{}

bool SceneAligner::computeTransformAndDestinationSeparators(
    SfM_Data &destination_sfm_data,
    const SfM_Data &sfm_data_first,
    const SfM_Data &sfm_data_second,
    std::vector<double> &second_base_node_pose,
    double &scaling_factor,
    const std::vector<size_t> & common_track_ids)
{
  // check precondition for function to run normally
  if (!scenesAreAlignable(sfm_data_first, sfm_data_second, common_track_ids))
  {
    return false;
  }

  ceres::Problem problem;

  // note : configureProblem is a virtual method !
  configureProblem(
    problem,
    destination_sfm_data,
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
  }
}

void SceneAligner::configureProblem(ceres::Problem & problem,
    SfM_Data &destination_sfm_data,
    const SfM_Data & sfm_data_first,
    const SfM_Data & sfm_data_second,
    std::vector<double> & second_base_node_pose,
    double & scaling_factor,
    const std::vector<size_t> & common_track_ids
    )
{
  const Landmarks & landmarks_first = sfm_data_first.structure;
  const Landmarks & landmarks_second = sfm_data_second.structure;
  Landmarks & destination_landmarks = destination_sfm_data.structure;

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

bool MergeScenesUsingCommonTracks
(SfM_Data & destination_sfm_data,
  const SfM_Data & sfm_data_first, // first submap scene
  const SfM_Data & sfm_data_second, // second submap scene
  const std::vector<size_t> & common_track_ids,
  SceneAligner *smap_aligner
)
{
  // initialize destination sfm data
  destination_sfm_data.intrinsics = sfm_data_first.intrinsics;
  for (const auto & intrinsic : sfm_data_second.GetIntrinsics())
  {
    if (destination_sfm_data.intrinsics.find(intrinsic.first) == destination_sfm_data.intrinsics.end())
      destination_sfm_data.intrinsics[intrinsic.first] = intrinsic.second;
  }
  destination_sfm_data.poses = sfm_data_first.poses;
  destination_sfm_data.structure = sfm_data_first.structure;

  // initialize the transformation parameters
  std::vector<double> second_base_node_pose(6, 0.0);
  double scaling_factor(1.0);

  const bool alignment_successful =
      smap_aligner->computeTransformAndDestinationSeparators(
        destination_sfm_data, sfm_data_first, sfm_data_second,
        second_base_node_pose, scaling_factor, common_track_ids);

  if (alignment_successful)
  {
    Vec3 second_base_node_t = Vec3(second_base_node_pose[3],second_base_node_pose[4],second_base_node_pose[5]);
    Mat3 second_base_node_RMat;
    ceres::AngleAxisToRotationMatrix(&second_base_node_pose[0], second_base_node_RMat.data());

    std::cout << " new aligned base node coords : \n"
    << second_base_node_t << std::endl
    << second_base_node_RMat << std::endl
    << " scale : " << scaling_factor << std::endl;

    // update the landmarks and poses from the second
    // scene transformed into the destination scene.
    transformSfMDataScene(destination_sfm_data, sfm_data_second,
                          second_base_node_RMat, second_base_node_t, scaling_factor);
    return true;
  }
  else
  {
    return false;
  }
}

/**
 * @brief transforms the content of a sfm data scene into another using a
 * rotation, a translation, and a scale factor
 * @param the original sfm data scene
 * @return the destination sfm data scene
 * @param the rotation matrix
 * @param a translation vector
 * @param a scale factor
 * @note both poses and landmarks are transformed, note that they will be overwritten
 * into the destination scene !
 */
void transformSfMDataScene(
    SfM_Data & destination_sfm_data,
    const SfM_Data & original_sfm_data,
    const Mat3 & rotation,
    const Vec3 & translation,
    const double scaling_factor)
{
  // update camera poses
  const Poses & original_poses = original_sfm_data.poses;
  Poses & destination_poses = destination_sfm_data.poses;
  for (const auto & pose : original_poses)
  {
    geometry::Pose3 new_pose = pose.second;
    const Vec3 original_position = new_pose.center();
    const Mat3 original_rotation = new_pose.rotation();

    new_pose = geometry::Pose3(original_rotation * rotation.transpose(),
        (rotation * ((1.0/scaling_factor)*original_position) + translation));

    // push to parent sfm data
    destination_poses[pose.first] = new_pose;
  }

  // update landmarks
  const Landmarks & original_landmarks = original_sfm_data.structure;
  Landmarks & destination_landmarks = destination_sfm_data.structure;
  for (const auto & landmark : original_landmarks)
  {
    const Vec3 original_position = landmark.second.X;

    Landmark new_landmark = landmark.second;
    new_landmark.X = rotation * ((1.0/scaling_factor) * original_position) + translation;

    // push to destination sfm data
    destination_landmarks[landmark.first] = new_landmark;
  }
}

} // namespace sfm
} // namespace openMVG
