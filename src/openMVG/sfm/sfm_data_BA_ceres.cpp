// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#include "ceres/rotation.h"

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

/// Create the appropriate cost functor according the provided input camera intrinsic model
ceres::CostFunction * IntrinsicsToCostFunction(IntrinsicBase * intrinsic, const Vec2 & observation)
{
  switch(intrinsic->getType())
  {
    case PINHOLE_CAMERA:
      return new ceres::AutoDiffCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic, 2, 3, 6, 3>(
        new ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()));
    break;
    case PINHOLE_CAMERA_RADIAL1:
      return new ceres::AutoDiffCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1, 2, 4, 6, 3>(
        new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data()));
    break;
    case PINHOLE_CAMERA_RADIAL3:
      return new ceres::AutoDiffCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 3>(
        new ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data()));
    break;
    case PINHOLE_CAMERA_BROWN:
      return new ceres::AutoDiffCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2, 2, 8, 6, 3>(
        new ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data()));
    default:
      return NULL;
  }
}


Bundle_Adjustment_Ceres::BA_options::BA_options(const bool bVerbose, bool bmultithreaded)
  :_bVerbose(bVerbose),
   _nbThreads(1)
{
  #ifdef OPENMVG_USE_OPENMP
    _nbThreads = omp_get_max_threads();
  #endif // OPENMVG_USE_OPENMP
  if (!bmultithreaded)
    _nbThreads = 1;

  _bCeres_Summary = false;

  // Default configuration use a DENSE representation
  _linear_solver_type = ceres::DENSE_SCHUR;
  _preconditioner_type = ceres::JACOBI;
  // If Sparse linear solver are available
  // Descending priority order by efficiency (SUITE_SPARSE > CX_SPARSE > EIGEN_SPARSE)
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
  {
    _sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
    _linear_solver_type = ceres::SPARSE_SCHUR;
  }
  else
  {
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
    {
      _sparse_linear_algebra_library_type = ceres::CX_SPARSE;
      _linear_solver_type = ceres::SPARSE_SCHUR;
    }
    else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
    {
      _sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
      _linear_solver_type = ceres::SPARSE_SCHUR;
    }
  }
}


Bundle_Adjustment_Ceres::Bundle_Adjustment_Ceres(
  Bundle_Adjustment_Ceres::BA_options options)
  : _openMVG_options(options)
{}

bool Bundle_Adjustment_Ceres::Adjust(
  SfM_Data & sfm_data,     // the SfM scene to refine
  bool bRefineRotations,   // tell if pose rotations will be refined
  bool bRefineTranslations,// tell if the pose translation will be refined
  bool bRefineIntrinsics,  // tell if the camera intrinsic will be refined
  bool bRefineStructure)   // tell if the structure will be refined
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

  // Setup Poses data & subparametrization
  for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
  {
    const IndexT indexPose = itPose->first;

    const Pose3 & pose = itPose->second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    map_poses[indexPose].reserve(6); //angleAxis + translation
    map_poses[indexPose].push_back(angleAxis[0]);
    map_poses[indexPose].push_back(angleAxis[1]);
    map_poses[indexPose].push_back(angleAxis[2]);
    map_poses[indexPose].push_back(t(0));
    map_poses[indexPose].push_back(t(1));
    map_poses[indexPose].push_back(t(2));

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (!bRefineTranslations && !bRefineRotations)
    {
      //set the whole parameter block as constant for best performance.
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  {
      // Subset parametrization
      std::vector<int> vec_constant_extrinsic;
      if(!bRefineRotations)
      {
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      if(!bRefineTranslations)
      {
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const IndexT indexCam = itIntrinsic->first;

    if (isValid(itIntrinsic->second->getType()))
    {
      map_intrinsics[indexCam] = itIntrinsic->second->getParams();

      double * parameter_block = &map_intrinsics[indexCam][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
      if (!bRefineIntrinsics)
      {
        //set the whole parameter block as constant for best performance.
        problem.SetParameterBlockConstant(parameter_block);
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));
  // TODO: make the LOSS function and the parameter an option

  // For all visibility add reprojections errors:
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views[itObs->first].get();

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
          iterTracks->second.X.data()); //Do we need to copy 3D point to avoid false motion, if failure ?
    }
    if (!bRefineStructure)
      problem.SetParameterBlockConstant(iterTracks->second.X.data());
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options options;
  options.preconditioner_type = _openMVG_options._preconditioner_type;
  options.linear_solver_type = _openMVG_options._linear_solver_type;
  options.sparse_linear_algebra_library_type = _openMVG_options._sparse_linear_algebra_library_type;
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;
  options.num_threads = _openMVG_options._nbThreads;
  options.num_linear_solver_threads = _openMVG_options._nbThreads;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  if (_openMVG_options._bCeres_Summary)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (_openMVG_options._bVerbose)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (_openMVG_options._bVerbose)
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
        << std::endl;
    }

    // Update camera poses with refined data
    if (bRefineRotations || bRefineTranslations)
    {
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
    }

    // Update camera intrinsics with refined data
    if (bRefineIntrinsics)
    {
      for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
        itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
      {
        const IndexT indexCam = itIntrinsic->first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        itIntrinsic->second.get()->updateFromParams(vec_params);
      }
    }
    return true;
  }
}

} // namespace sfm
} // namespace openMVG

