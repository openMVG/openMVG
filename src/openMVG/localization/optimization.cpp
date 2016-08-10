/* 
 * File:   optimization.cpp
 * Author: sgaspari
 * 
 * Created on October 23, 2015, 12:02 AM
 */

#include "optimization.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include <openMVG/sfm/sfm_data_BA_ceres.hpp>
#include <openMVG/rig/rig_BA_ceres.hpp>
#include <openMVG/logger.hpp>


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum.hpp>

namespace openMVG{
namespace localization{

bool refineSequence(std::vector<LocalizationResult> & vec_localizationResult,
                    bool allTheSameIntrinsics /*= true*/,
                    bool b_refine_intrinsic /*= true*/,
                    bool b_no_distortion /*= false*/,
                    bool b_refine_pose /*= true*/,
                    bool b_refine_structure /*= false*/,
                    const std::string outputFilename /*= ""*/,
                    std::size_t minPointVisibility /*=0*/)
{
  
  const std::size_t numViews = vec_localizationResult.size();
  assert(numViews > 0 );
   
  // the id for the instrinsic group
  IndexT intrinsicID = 0;
    
  // Setup a tiny SfM scene with the corresponding 2D-3D data
  sfm::SfM_Data tinyScene;
  
  // if we have only one camera just set the intrinsics group once for all
  if(allTheSameIntrinsics)
  {
    // find the first valid localization result and use its intrinsics
    std::size_t intrinsicIndex = 0;
    for(std::size_t viewID = 0; viewID < numViews; ++viewID, ++intrinsicIndex)
    {
      if(vec_localizationResult[viewID].isValid())
        break;
    }
    // it may be the case that all the localization result are invalid...
    if(!vec_localizationResult[intrinsicIndex].isValid())
    {
      POPART_CERR("Apparently all the vec_localizationResult are invalid! Aborting...");
      return false;
    }
    
    POPART_CERR("allTheSameIntrinsics mode: using the intrinsics of the " << intrinsicIndex << " result");
    
    cameras::Pinhole_Intrinsic_Radial_K3* currIntrinsics = &vec_localizationResult[intrinsicIndex].getIntrinsics();
    
    if(b_no_distortion)
    {
      // no distortion refinement
      POPART_COUT("Optical distortion won't be considered");
      // just add a simple pinhole camera with the same K as the input camera
      Vec2 pp = currIntrinsics->principal_point();
      tinyScene.intrinsics[intrinsicID] = std::make_shared<cameras::Pinhole_Intrinsic>(currIntrinsics->_w, currIntrinsics->_h, currIntrinsics->focal(), pp(0), pp(1));
    }
    else
    {
      // intrinsic (the shared_ptr does not take the ownership, will not release the input pointer)
      tinyScene.intrinsics[intrinsicID] = std::shared_ptr<cameras::Pinhole_Intrinsic_Radial_K3>(currIntrinsics, [](cameras::Pinhole_Intrinsic_Radial_K3*){});
      POPART_COUT("Type of intrinsics " <<tinyScene.intrinsics[0].get()->getType());
    }
  }
  
  {
    // just debugging -- print reprojection errors -- this block can be safely removed/commented out
    for(size_t viewID = 0; viewID < numViews; ++viewID)
    {
      const LocalizationResult &currResult = vec_localizationResult[viewID];
      if(!currResult.isValid())
      {
        continue;
      }
      
      const Mat2X residuals = currResult.computeResiduals();
      
      const auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
      POPART_COUT("View " << viewID << " RMSE = " << std::sqrt(sqrErrors.mean()) 
              << " min = " << std::sqrt(sqrErrors.minCoeff()) 
              << " mean = " << std::sqrt(sqrErrors.mean())
              << " max = " << std::sqrt(sqrErrors.maxCoeff()));
    }
  }
  
  for(size_t viewID = 0; viewID < numViews; ++viewID)
  {
    LocalizationResult &currResult = vec_localizationResult[viewID];
    // skip invalid poses
    if(!currResult.isValid())
    {
      POPART_COUT("\n*****\nskipping invalid View " << viewID);
      continue;
    }
    
//    POPART_COUT("\n*****\nView " << viewID);
    // view
    tinyScene.views.insert( std::make_pair(viewID, std::make_shared<sfm::View>("",viewID, intrinsicID, viewID)));
    // pose
    tinyScene.poses[viewID] = currResult.getPose();
    
    if(!allTheSameIntrinsics)
    {
      cameras::Pinhole_Intrinsic_Radial_K3* currIntrinsics = &currResult.getIntrinsics();
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
//      POPART_COUT("inlier " << idx << " is land " << landmarkID << " and feat " << featID);
      // get the corresponding feature
      const Vec2 &feature = currResult.getPt2D().col(idx);
      // check if the point exists already
      if(tinyScene.structure.count(landmarkID))
      {
        // normally there should be no other features already associated to this
        // 3D point in this view
        if(tinyScene.structure[landmarkID].obs.count(viewID) != 0)
        {
          // this is weird but it could happen when two features are really close to each other (?)
          POPART_COUT("Point 3D " << landmarkID << " has multiple features " 
                  << " in the same view " << viewID << " , size of obs: " 
                  << tinyScene.structure[landmarkID].obs.size() );
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

//  {
//    POPART_COUT("Number of 3D-2D associations before filtering " << tinyScene.structure.size());
//    sfm::Landmarks &landmarks = tinyScene.structure;
//    for(sfm::Landmarks::iterator it = landmarks.begin(), ite = landmarks.end(); it != ite;)
//    {
//      if(it->second.obs.size() < 5)
//      {
//         it = landmarks.erase(it);
//      }
//      else
//        ++it;
//    }
//  }
  
  {
    // just debugging some stats -- this block can be safely removed/commented out
    
    POPART_COUT("Number of 3D-2D associations " << tinyScene.structure.size());
    
    std::size_t maxObs = 0;
    for(const auto landmark : tinyScene.GetLandmarks() )
    {
      if(landmark.second.obs.size() > maxObs)
        maxObs = landmark.second.obs.size();
    }
    namespace bacc = boost::accumulators;
    bacc::accumulator_set<std::size_t, bacc::stats<bacc::tag::mean, bacc::tag::min, bacc::tag::max, bacc::tag::sum > > stats;
    std::vector<std::size_t> hist(maxObs+1, 0);
    for(const auto landmark : tinyScene.GetLandmarks() )
    {
      const std::size_t nobs = landmark.second.obs.size();
      assert(nobs < hist.size());
      stats(nobs);
      hist[nobs]++;
    }
    POPART_COUT("Min number of observations per point:   " << bacc::min(stats) );
    POPART_COUT("Mean number of observations per point:   " << bacc::mean(stats) );
    POPART_COUT("Max number of observations per point:   " << bacc::max(stats) );
    
    for( int i = 0; i < hist.size(); i++ ) 
    {
      POPART_COUT("Points with " << i << " observations: " << hist[i]); 
    }

    // just debugging stuff
    if(allTheSameIntrinsics)
    {
      std::vector<double> params = tinyScene.intrinsics[0].get()->getParams();
      POPART_COUT("K before bundle: " << params[0] << " " << params[1] << " "<< params[2]);
      if(params.size() == 6)
        POPART_COUT("Distortion before bundle: " << params[3] << " " << params[4] << " "<< params[5]);
    }
  }

  // filter out the 3D points having too few observations.
  if(minPointVisibility > 0)
  {
    auto &landmarks = tinyScene.structure;
    auto iter = landmarks.begin();
    auto endIter = landmarks.end();

    for(; iter != endIter;)
    {
      if(iter->second.obs.size() < minPointVisibility)
      {
        iter = landmarks.erase(iter);
      }
      else
      {
        ++iter;
      }
    }
  }
  
  if(!outputFilename.empty())
  {
    const std::string outfile = outputFilename+".BEFORE.json";
    if(!sfm::Save(tinyScene, outfile, sfm::ESfM_Data::ALL))
      POPART_CERR("Could not save " << outfile);
  }

  sfm::Bundle_Adjustment_Ceres bundle_adjustment_obj;
  sfm::BA_Refine refineOptions = sfm::BA_REFINE_NONE;
  if(b_refine_pose)
    refineOptions |= sfm::BA_REFINE_ROTATION | sfm::BA_REFINE_TRANSLATION;
  if(b_refine_intrinsic)
    refineOptions |= sfm::BA_REFINE_INTRINSICS_ALL;
  if(b_refine_structure)
    refineOptions |= sfm::BA_REFINE_STRUCTURE;

  const bool b_BA_Status = bundle_adjustment_obj.Adjust(tinyScene, refineOptions);
  if(b_BA_Status)
  {
    // get back the results and update the localization result with the refined pose
    for(const auto &pose : tinyScene.poses)
    {
      const IndexT idPose = pose.first;
      vec_localizationResult[idPose].setPose(pose.second);
    }

    if(!outputFilename.empty())
    {
      const std::string outfile = outputFilename+".AFTER.json";
      if(!sfm::Save(tinyScene, outfile, sfm::ESfM_Data::ALL))
        POPART_CERR("Could not save " << outfile);
    }
  }
  
  if(allTheSameIntrinsics)
  {
    // if we used the same intrinsics for all the localization results we need to
    // update the intrinsics of each localization result
    
    // get its optimized parameters
    std::vector<double> params = tinyScene.intrinsics[0].get()->getParams();
    POPART_COUT("Type of intrinsics " <<tinyScene.intrinsics[0].get()->getType());
    if(params.size() == 3)
    {
      // this means that the b_no_distortion has been passed
      // set distortion to 0
      params.push_back(0);
      params.push_back(0);
      params.push_back(0);
    }
    assert(params.size() == 6);
    POPART_COUT("K after bundle: " << params[0] << " " << params[1] << " "<< params[2]);
    POPART_COUT("Distortion after bundle " << params[3] << " " << params[4] << " "<< params[5]);

    // update the intrinsics of the each localization result
    for(size_t viewID = 0; viewID < numViews; ++viewID)
    {
      LocalizationResult &currResult = vec_localizationResult[viewID];
      if(!currResult.isValid())
      {
        continue;
      }
      currResult.updateIntrinsics(params);
      
      // just debugging -- print reprojection errors
      const Mat2X residuals = currResult.computeResiduals();
      
      const auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
      POPART_COUT("View " << viewID << " RMSE = " << std::sqrt(sqrErrors.mean()) 
              << " min = " << std::sqrt(sqrErrors.minCoeff()) 
              << " mean = " << std::sqrt(sqrErrors.mean())
              << " max = " << std::sqrt(sqrErrors.maxCoeff()));
    }
    
  }
  
  return b_BA_Status;
}

bool refineRigPose(const std::vector<geometry::Pose3 > &vec_subPoses,
                   const std::vector<localization::LocalizationResult> vec_localizationResults,
                   geometry::Pose3 & rigPose)
{
  const std::size_t numCameras = vec_localizationResults.size();
  assert(vec_subPoses.size() == numCameras - 1);
  
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
  for(std::size_t iLocalizer = 0; iLocalizer < numCameras; ++iLocalizer)
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

      if(cost_function)
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
  options.minimizer_progress_to_stdout = openMVG_options._bVerbose;
  options.logging_type = ceres::SILENT;
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

bool refineRigPose(const std::vector<Mat> &pts2d,
                   const std::vector<Mat> &pts3d,
                   const std::vector<std::vector<std::size_t> > &inliers,
                   const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                   const std::vector<geometry::Pose3 > &vec_subPoses,
                   geometry::Pose3 &rigPose)
{
  const size_t numCameras = pts2d.size();
  assert(pts3d.size() == numCameras);
  assert(vec_queryIntrinsics.size() == numCameras);
  assert(inliers.size() == numCameras);
  assert(vec_subPoses.size() == numCameras - 1);
  
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
  for(size_t cam = 0; cam < numCameras; ++cam)
  {

    // Get the inliers 3D points
    const Mat & points3D = pts3d[cam];
    assert(points3D.rows() == 3);
    // Get their image locations (also referred as observations)
    const Mat & points2D = pts2d[cam];
    assert(points2D.rows() == 2);

    // Add a residual block for all inliers
    for(const IndexT iPoint : inliers[cam])
    {
      assert(iPoint < points2D.cols());
      assert(iPoint < points3D.cols());
      
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observations
      // and the 3D point and compares the reprojection against the observation.
      ceres::CostFunction* cost_function;

      // Vector-2 residual, pose of the rig parameterized by 6 parameters
      //                  + relative pose of the secondary camera parameterized by 6 parameters
      
      geometry::Pose3 subPose;
      // if it is not the main camera (whose subpose is the identity)
      if(cam != 0)
      {
        subPose = vec_subPoses[cam - 1];
      }

      cost_function = new ceres::AutoDiffCostFunction<rig::ResidualErrorSecondaryCameraFixedRelativeFunctor, 2, 6>(
              new rig::ResidualErrorSecondaryCameraFixedRelativeFunctor(vec_queryIntrinsics[cam],
                                                                   points2D.col(iPoint),
                                                                   points3D.col(iPoint),
                                                                   subPose));

      if(cost_function)
      {
        problem.AddResidualBlock(cost_function,
                                 p_LossFunction,
                                 mainPose);
      }
      else
      {
        POPART_CERR("Fail in adding residual block for the " << cam 
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
            << " #cameras: " << numCameras << "\n"
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

