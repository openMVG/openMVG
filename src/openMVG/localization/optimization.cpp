/* 
 * File:   optimization.cpp
 * Author: sgaspari
 * 
 * Created on October 23, 2015, 12:02 AM
 */

#include "optimization.hpp"
#include <openMVG/sfm/sfm_data_BA_ceres.hpp>
#include <openMVG/rig/rig_BA_ceres.hpp>
#include <openMVG/logger.hpp>

namespace openMVG{
namespace localization{

bool refineSequence(cameras::Pinhole_Intrinsic_Radial_K3 *intrinsics,
                    std::vector<LocalizationResult> & vec_localizationResult,
                    bool b_refine_pose /*= true*/,
                    bool b_refine_intrinsic /*= true*/,
                    bool b_refine_structure /*= false*/)
{
  std::vector<cameras::Pinhole_Intrinsic_Radial_K3* > vec_intrinsics;
  vec_intrinsics.push_back(intrinsics);
  return refineSequence(vec_intrinsics, vec_localizationResult);
}

bool refineSequence(std::vector<cameras::Pinhole_Intrinsic_Radial_K3* > vec_intrinsics,
                    std::vector<LocalizationResult> & vec_localizationResult,
                    bool b_refine_pose /*= true*/,
                    bool b_refine_intrinsic /*= true*/,
                    bool b_refine_structure /*= false*/)
{
   
  // flags for BA
//  const bool b_refine_pose = true;
//  const bool b_refine_intrinsic = true;
//  const bool b_refine_structure = false;
  
  // vec_intrinsics must be either of the same size of localization result (a 
  // camera for each found pose) or it must contain only 1 element, meaning that 
  // it's the same camera for the whole sequence
  assert(vec_intrinsics.size()==vec_localizationResult.size() || vec_intrinsics.size()==1);
  
  const size_t numViews = vec_localizationResult.size();
  const bool singleCamera = vec_intrinsics.size()==1;
  
  // the id for the instrinsic group
  IndexT intrinsicID = 0;
    
  // Setup a tiny SfM scene with the corresponding 2D-3D data
  sfm::SfM_Data tinyScene;
  
  // if we have only one camera just set the intrinsics group once for all
  if(singleCamera)
  {
    // intrinsic (the shared_ptr does not take the ownership, will not release the input pointer)
    tinyScene.intrinsics[intrinsicID] = std::shared_ptr<cameras::Pinhole_Intrinsic_Radial_K3>(vec_intrinsics[0], [](cameras::Pinhole_Intrinsic_Radial_K3*){});
  }
  
  for(size_t viewID = 0; viewID < numViews; ++viewID)
  {
    const LocalizationResult &currResult = vec_localizationResult[viewID];
    cameras::Pinhole_Intrinsic_Radial_K3* currIntrinsics = vec_intrinsics[viewID];
    // skip invalid poses
    if(!currResult.isValid())
    {
      std::cout << "\n*****\nskipping invalid View " << viewID << std::endl;
      continue;
    }
    
    std::cout << "\n*****\nView " << viewID << std::endl;
    // view
    tinyScene.views.insert( std::make_pair(viewID, std::make_shared<sfm::View>("",viewID, intrinsicID, viewID)));
    // pose
    tinyScene.poses[viewID] = currResult.getPose();
    
    if(!singleCamera)
    {
       // intrinsic (the shared_ptr does not take the ownership, will not release the input pointer)
      tinyScene.intrinsics[intrinsicID] = std::shared_ptr<cameras::Pinhole_Intrinsic_Radial_K3>(currIntrinsics, [](cameras::Pinhole_Intrinsic_Radial_K3*){});
      ++intrinsicID;
    }
    
    // structure data (2D-3D correspondences)
    const vector<pair<IndexT, IndexT> > &currentIDs = currResult.getIndMatch3D2D();
    
    for(const size_t idx : currResult.getInliers() )
    {
      // the idx should be in the size range of the data
      assert(idx < currResult.getPt3D().cols());
      // get the corresponding 3D point (landmark) ID
      const IndexT landmarkID = currentIDs[idx].first;
      // get the corresponding 2D point ID
      const IndexT featID = currentIDs[idx].second;
      std::cout << "inlier " << idx << " is land " << landmarkID << " and feat " << featID << std::endl;
      // get the corresponding feature
      const Vec2 &feature = currResult.getPt2D().col(idx);
      // check if the point exists already
      if(tinyScene.structure.count(landmarkID))
      {
        // normally there should be no other features already associated to this
        // 3D point in this view
//        assert(tinyScene.structure[landmarkID].obs.count(viewID) == 0);
        if(tinyScene.structure[landmarkID].obs.count(viewID) != 0)
        {
          // this is weird but it could happen when two features are really close to each other (?)
          std::cout << "Point 3D " << landmarkID << " has multiple features " << tinyScene.structure[landmarkID].obs.size() << " in the same view " << viewID << " size "  << std::endl; 
          continue;
        }
        
        // the 3D point exists already, add the observation
        tinyScene.structure[landmarkID].obs[viewID] =  sfm::Observation(feature, featID);
      }
      else
      {
        // create a new 3D point
        sfm::Landmark landmark;
        landmark.X = currResult.getPt3D().col(idx);
        landmark.obs[viewID] = sfm::Observation(feature, featID);
        tinyScene.structure[landmarkID] = std::move(landmark);
      }
    }
  }
  POPART_COUT("Number of 3D-2D associations " << tinyScene.structure.size());
  
  if(singleCamera)
  {
    // just debugging stuff
    const cameras::Pinhole_Intrinsic_Radial_K3* intrinsics = vec_intrinsics[0];
    std::vector<double> params = intrinsics->getParams();
    POPART_COUT("K before bundle:" << params[0] << " " << params[1] << " "<< params[2]);
    POPART_COUT("Distortion before bundle" << params[3] << " " << params[4] << " "<< params[5]);
  }

  sfm::Bundle_Adjustment_Ceres bundle_adjustment_obj;
  const bool b_BA_Status = bundle_adjustment_obj.Adjust(tinyScene, b_refine_pose, b_refine_pose, b_refine_intrinsic, b_refine_structure);
  if(b_BA_Status)
  {
    // get back the results
    for(const auto &pose : tinyScene.poses)
    {
      const IndexT idPose = pose.first;
      vec_localizationResult[idPose].setPose(pose.second);
    }
  }
  
  if(singleCamera)
  {
    // just debugging stuff
    const cameras::Pinhole_Intrinsic_Radial_K3* intrinsics = vec_intrinsics[0];
    std::vector<double> params = intrinsics->getParams();
    POPART_COUT("K after bundle:" << params[0] << " " << params[1] << " "<< params[2]);
    POPART_COUT("Distortion after bundle" << params[3] << " " << params[4] << " "<< params[5]);
  }
  
  return b_BA_Status;
}

bool refineRigPose(const std::vector<geometry::Pose3 > &vec_subPoses,
                   const std::vector<localization::LocalizationResult> vec_localizationResults,
                   geometry::Pose3 & rigPose)
{
  ceres::Problem problem;
  
  const openMVG::Mat3 & R = rigPose.rotation();
  const openMVG::Vec3 & t = rigPose.translation();

  double mainPose[6];
  ceres::RotationMatrixToAngleAxis((const double*)R.data(), mainPose);

  mainPose[3] = t(0);
  mainPose[4] = t(1);
  mainPose[5] = t(2);
  problem.AddParameterBlock(mainPose, 6);


  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = nullptr;//new ceres::HuberLoss(Square(4.0));
  // todo: make the LOSS function and the parameter an option

  // For all visibility add reprojections errors:
  for(int iLocalizer = 0; iLocalizer < vec_localizationResults.size(); ++iLocalizer)
  {
    const localization::LocalizationResult & localizationResult = vec_localizationResults[iLocalizer];

    if(!localizationResult.isValid())
    {
      POPART_COUT("Skipping camera " << iLocalizer << " as it has not been localized");
      continue;
    }
    // Get the inliers 3D points
    const Mat & points3D = localizationResult.getPt3D();
    // Get their image locations (also referred as observations)
    const Mat & points2D = localizationResult.getPt2D();

    // Add a residual block for all inliers
    for(const IndexT iPoint : localizationResult.getInliers())
    {
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observations
      // and the 3D point and compares the reprojection against the observation.
      ceres::CostFunction* cost_function;

      // Vector-2 residual, pose of the rig parameterized by 6 parameters
      //                  + relative pose of the secondary camera parameterized by 6 parameters
      
      geometry::Pose3 subPose;
      // if it is not the main camera (whose subpose is the identity)
      if(iLocalizer != 0)
      {
        subPose = vec_subPoses[iLocalizer - 1];
      }

      cost_function = new ceres::AutoDiffCostFunction<rig::ResidualErrorSecondaryCameraFixedRelativeFunctor, 2, 6>(
              new rig::ResidualErrorSecondaryCameraFixedRelativeFunctor(localizationResult.getIntrinsics(),
                                                                   points2D.col(iPoint),
                                                                   points3D.col(iPoint),
                                                                   subPose));

      if(!cost_function)
      {
        problem.AddResidualBlock(cost_function,
                                 p_LossFunction,
                                 mainPose);
      }
      else
      {
        POPART_CERR("Fail in adding residual block for the " << iLocalizer 
                << " camera while adding point id " << iPoint);
      }
    }
  }

  // Configure a BA engine and run it
  // todo: Set the most appropriate options
  openMVG::sfm::Bundle_Adjustment_Ceres::BA_options openMVG_options; // Set all
  // the options field in our owm struct - unnecessary dependancy to openMVG here
  
  ceres::Solver::Options options;
  
  options.preconditioner_type = openMVG_options._preconditioner_type;
  options.linear_solver_type = openMVG_options._linear_solver_type;
  options.sparse_linear_algebra_library_type = openMVG_options._sparse_linear_algebra_library_type;
  options.minimizer_progress_to_stdout = true;
  //options.logging_type = ceres::SILENT;
  options.num_threads = 1;//openMVG_options._nbThreads;
  options.num_linear_solver_threads = 1;//openMVG_options._nbThreads;
  
  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  
  if (openMVG_options._bCeres_Summary)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (openMVG_options._bVerbose)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }

  if(openMVG_options._bVerbose)
  {
    // Display statistics about the minimization
    std::cout << std::endl
            << "Bundle Adjustment statistics (approximated RMSE):\n"
            << " #localizers: " << vec_localizationResults.size() << "\n"
            << " #residuals: " << summary.num_residuals << "\n"
            << " Initial RMSE: " << std::sqrt(summary.initial_cost / summary.num_residuals) << "\n"
            << " Final RMSE: " << std::sqrt(summary.final_cost / summary.num_residuals) << "\n"
            << std::endl;
  }

  // update the rigPose 
  openMVG::Mat3 R_refined;
  ceres::AngleAxisToRotationMatrix(mainPose, R_refined.data());
  openMVG::Vec3 t_refined(mainPose[3], mainPose[4], mainPose[5]);
  // Push the optimized pose
  rigPose = geometry::Pose3(R_refined, -R_refined.transpose() * t_refined);

//  displayRelativePoseReprojection(geometry::Pose3(openMVG::Mat3::Identity(), openMVG::Vec3::Zero()), 0);

  // @todo do we want to update pose inside the LocalizationResults 

  return true;
}


} //namespace localization
} //namespace openMVG

