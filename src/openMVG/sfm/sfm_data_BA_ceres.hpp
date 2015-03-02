// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_CERES_HPP
#define OPENMVG_SFM_DATA_BA_CERES_HPP

#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"
#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace openMVG {

class Bundle_Adjustment_Ceres : public Bundle_Adjustment
{
  public:
  bool bAdjust(
    SfM_Data & sfm_data, // the SfM scene to refine
    bool bRefineRotations = true, // tell if pose rotations will be refined
    bool bRefineTranslations = true, // tell if the pose translation will be refined
    bool bRefineIntrinsics = true) // tell if the camera intrinsic will be refined)
  {
    // Add camera parameters
    // - intrinsics
    // - poses [R|t]

    using namespace openMVG::bundle_adjustment;

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
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
      if (!bRefineTranslations && !bRefineIntrinsics)
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

      switch (itIntrinsic->second->getType())
      {
        case PINHOLE_CAMERA:
        {
          const std::vector<double> vec_params = itIntrinsic->second->getParams();
          map_intrinsics[indexCam].reserve(3);
          map_intrinsics[indexCam].push_back(vec_params[0]);
          map_intrinsics[indexCam].push_back(vec_params[1]);
          map_intrinsics[indexCam].push_back(vec_params[2]);

          double * parameter_block = &map_intrinsics[indexCam][0];
          problem.AddParameterBlock(parameter_block, 3);
          if (!bRefineIntrinsics)
          {
            //set the whole parameter block as constant for best performance.
            problem.SetParameterBlockConstant(parameter_block);
          }
        }
        break;
      }
    }

    // Set a LossFunction to be less penalized by false measurements
    //  - set it to NULL if you don't want use a lossFunction.
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(2.0));
    // TODO: make the LOSS function and the parameter an option

    // For all visibility add reprojections errors:
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
      iterTracks!= sfm_data.structure.end(); ++iterTracks)
    {
      const Observations & obs = iterTracks->second.obs;
      for(Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View & view = sfm_data.views[itObs->first];

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function = NULL;
        switch(sfm_data.intrinsics[view.id_intrinsic]->getType())
        {
          case PINHOLE_CAMERA:
          {
            cost_function =
             new ceres::AutoDiffCostFunction<pinhole_reprojectionError::ErrorFunc_Refine_Intrinsic_Motion_3DPoints, 2, 3, 6, 3>(
              new pinhole_reprojectionError::ErrorFunc_Refine_Intrinsic_Motion_3DPoints(itObs->second.x.data()));
          }
          break;
          default:
            std::cout << "Unsupported camera type: please create the appropriate ceres functor." << std::endl;
        }
        if (cost_function)
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[view.id_intrinsic][0],
            &map_poses[view.id_pose][0],
            iterTracks->second.X.data()); //Do we need to copy 3D point to avoid false motion, if failure ?
      }
    }

    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options options;
    options.preconditioner_type = ceres::JACOBI;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
    else
    {
      if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
        options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
      else
      {
        // No sparse backend for Ceres.
        // Use dense solving
        options.linear_solver_type = ceres::DENSE_SCHUR;
      }
    }
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;
  #ifdef OPENMVG_USE_OPENMP
    options.num_threads = omp_get_max_threads();
    options.num_linear_solver_threads = omp_get_max_threads();
  #endif // OPENMVG_USE_OPENMP

    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    //std::cout << summary.FullReport() << std::endl;

    // If no error, get back refined parameters
    if (!summary.IsSolutionUsable())
    {
      std::cout << "Bundle Adjustment failed." << std::endl;
      return false;
    }
    else // Solution is usable
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << std::endl;

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
      return true;
    }
  }
};
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_CERES_HPP
