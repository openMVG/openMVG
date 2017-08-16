#include "openMVG/sfm/pipelines/hierarchical_hyper/submap_BA_fixed_separators.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace openMVG{
namespace sfm{

  Bundle_Adjustment_Fixed_Separators::Bundle_Adjustment_Fixed_Separators
(
  Bundle_Adjustment_Ceres::BA_Ceres_options options
)
: ceres_options_(options)
{}

void Bundle_Adjustment_Fixed_Separators::configureProblem(ceres::Problem & problem, SfM_Data& sfm_data, const std::set<size_t>& separator_tracks_ids, Hash_Map<IndexT, std::vector<double> > & map_intrinsics, Hash_Map<IndexT, std::vector<double> > & map_poses)
{
  for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
  {
    const IndexT indexPose = itPose->first;

    const Pose3 & pose = itPose->second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
  }

  // Setup Intrinsics data
  // TODO : do we want to modify the intrinsics during this phase ???
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const IndexT indexCam = itIntrinsic->first;

    if (isValid(itIntrinsic->second->getType()))
    {
      map_intrinsics[indexCam] = itIntrinsic->second->getParams();
      double * parameter_block = &map_intrinsics[indexCam][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
      problem.SetParameterBlockConstant(parameter_block);
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(itObs->first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

      if (cost_function)
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[view->id_intrinsic][0],
          &map_poses[view->id_pose][0],
          iterTracks->second.X.data());
    }
    // if this is a separator track, we keep it fixed !
    if (separator_tracks_ids.count(iterTracks->first) != 0)
      problem.SetParameterBlockConstant(iterTracks->second.X.data());
  }
}

// Modified version of Bundle_Adjustment_Ceres::Adjust where we
// don't take any options as input (TODO should we though ??)
// and where we fix the position of the separator landmarks
bool Bundle_Adjustment_Fixed_Separators::Adjust(
   SfM_Data & sfm_data, // scene to refine
   const std::set<size_t> & separator_tracks_ids)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------
  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  configureProblem(problem, sfm_data, separator_tracks_ids, map_intrinsics, map_poses);

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
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

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
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
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
    }

    // Update camera poses with refined data
    for (Poses::iterator itPose = sfm_data.poses.begin();
      itPose != sfm_data.poses.end(); ++itPose)
    {
      const IndexT indexPose = itPose->first;

      Mat3 R_refined;
      ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
      Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
      // Update the pose
      Pose3 & pose = itPose->second;
      pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
    }

    // Update camera intrinsics with refined data
    for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
    {
      const IndexT indexCam = itIntrinsic->first;

      const std::vector<double> & vec_params = map_intrinsics[indexCam];
      itIntrinsic->second.get()->updateFromParams(vec_params);
    }
    // Structure is already updated directly if needed (no data wrapping)
    return true;
  }
}


}
}
