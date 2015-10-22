#include "Rig.hpp"
//#include "bundleAdjustmentCeresFunctor.hpp"

//#include <vision/cameraTracking/debug/visualDebug.hpp>
//#include <eigen3/Eigen/src/Core/MatrixBase.h>

#include <openMVG/sfm/sfm_data_BA_ceres.hpp>
#include <ceres/rotation.h>

namespace openMVG {
namespace rig {

Rig::~Rig()
{
}

void Rig::setTrackingResult(
        std::vector<localization::LocalizationResult> vLocalizationResults,
        std::size_t i)
{
  _vLocalizationResults.emplace(i, vLocalizationResults);
}

bool Rig::initializeCalibration()
{
  const std::size_t nCams = _vLocalizationResults.size();
  
  // Tracker of the main cameras
  std::vector<localization::LocalizationResult> & resMainCamera = _vLocalizationResults[0];
  
  // Clear all relative poses
  _vRelativePoses.clear();
  _vRelativePoses.reserve(nCams-1);
  
  // Loop over all witness cameras
  for (int i=1 ; i < nCams ; ++i)
  {
    // Perform the pose averaging over all relative pose between the main camera
    // (index 0) and the witness camera (index i)
    std::vector<localization::LocalizationResult> & resWitnessCamera = _vLocalizationResults[i];
    
    // vRelativePoses will store all the relative poses overall frames where both
    // the pose computation of the main camera and witness camera succeed
    std::vector<geometry::Pose3> vRelativePoses;
    vRelativePoses.reserve(resWitnessCamera.size());
    
    for(int iView=0 ; iView < resWitnessCamera.size() ; ++iView )
    {
      // Check that both pose computations succeed 
      if ( resMainCamera[iView].isValid() && resWitnessCamera[iView].isValid() )
      {
        const openMVG::Mat3 R1 = resMainCamera[iView].getPose().rotation();
        const openMVG::Vec3 t1 = resMainCamera[iView].getPose().translation();
        const openMVG::Mat3 R2 = resWitnessCamera[iView].getPose().rotation();
        const openMVG::Vec3 t2 = resWitnessCamera[iView].getPose().translation();
        
        const openMVG::Mat3 R12 = R2 * R1.transpose();
        const openMVG::Vec3 t12 = t2 - R12 * t1;
        
        const geometry::Pose3 relativePose( R12 , -R12.transpose()*t12 );
        
        vRelativePoses.push_back(relativePose);
      }
    }
    geometry::Pose3 optimalRelativePose;
    findOptimalPose(vRelativePoses, i, optimalRelativePose );
    //poseAveraging(vRelativePoses, averageRelativePose);
    _vRelativePoses.push_back(optimalRelativePose);
  }
  
  // Update all poses in all localization
  for (int iRelativePose = 0 ; iRelativePose < _vRelativePoses.size() ; ++iRelativePose )
  {
    std::size_t iRes = iRelativePose+1;
    for (int iView = 0 ; iView < _vLocalizationResults[iRes].size() ; ++iView )
    {
      if( _vLocalizationResults[iRes][iView].isValid() && _vLocalizationResults[0][iView].isValid() )
      {
        const geometry::Pose3 & relativePose = _vRelativePoses[iRelativePose];
        
        const openMVG::Mat3 R1 = _vLocalizationResults[0][iView].getPose().rotation();
        const openMVG::Vec3 t1 = _vLocalizationResults[0][iView].getPose().translation();
        
        const openMVG::Mat3 R12 = relativePose.rotation();
        const openMVG::Vec3 t12 = relativePose.translation();
        
        const openMVG::Mat3 R2 = R12 * R1;
        const openMVG::Vec3 t2 = R12 * t1 + t12 ;
        
        _vLocalizationResults[iRes][iView].setPose(geometry::Pose3( R2 , -R2.transpose() * t2 ));
      }
    }
  }
}

#if 0
geometry::Pose3 Rig::product(const geometry::Pose3 & poseA, const geometry::Pose3 & poseB)
{
  openMVG::Mat3 R = poseA.rotation()*poseB.rotation();
  openMVG::Vec3 t = poseA.rotation()*poseB.translation()+poseA.translation();
  return geometry::Pose3(R, -R.transpose()*t);
}

geometry::Pose3 Rig::productInv(const geometry::Pose3 & poseA, const geometry::Pose3 & poseB)
{
  openMVG::Vec3 t = -poseA.rotation().transpose()*poseA.translation();
  openMVG::Mat3 R = poseA.rotation().transpose();
  
  geometry::Pose3 poseC(R, -R.transpose()*t);
  return product(poseC,poseB);
}



double Rig::distance(openMVG::Vec3 va, openMVG::Vec3 vb)
{
  double d1 = va(0)-vb(0);
  double d2 = va(1)-vb(1);
  double d3 = va(2)-vb(2);
  return sqrt(d1*d1+d2*d2+d3*d3);
}

#endif

// From a set of relative pose, find the optimal one for a given tracker iTraker which
// minimize the reprojection errors over all images
void Rig::findOptimalPose(
        const std::vector<geometry::Pose3> & vPoses,
        std::size_t iRes,
        geometry::Pose3 & result )
{
  std::vector<localization::LocalizationResult> & resMainCamera = _vLocalizationResults[0];
  std::vector<localization::LocalizationResult> & resWitnessCamera = _vLocalizationResults[iRes];
  
  double minReprojError = std::numeric_limits<double>::max();
  double iMin = 0;
  
  for(int i=0 ; i < vPoses.size() ; ++i)
  {
    const geometry::Pose3 & relativePose = vPoses[i];

    double error = 0;
    for(int j=0 ; j < resWitnessCamera.size() ; ++j )
    {
      // Check that both pose computations succeed 
      if ( ( resMainCamera[j].isValid() ) && ( resWitnessCamera[j].isValid() ) )
      {
        const openMVG::Mat3 R1 = resMainCamera[j].getPose().rotation();
        const openMVG::Vec3 t1 = resMainCamera[j].getPose().translation();
        
        const openMVG::Mat3 R12 = relativePose.rotation();
        const openMVG::Vec3 t12 = relativePose.translation();
        
        const openMVG::Mat3 R2 = R12 * R1;
        const openMVG::Vec3 t2 = R12 * t1 + t12 ;
        
        const geometry::Pose3 pose( R2 , -R2.transpose() * t2 );
        
        // todo: uncomment:
        // error += reprojectionError(toOMVG(resWitnessCamera.intrinsics()[j].getIntrinsics().getK()), pose, tracker.imgPts()[j], tracker.pts()[j]);
      }
    }
    if ( error < minReprojError )
    {
      iMin = i;
      minReprojError = error;
    }           
  }
  result = vPoses[iMin];
  
  //displayRelativePoseReprojection(geometry::Pose3(openMVG::Mat3::Identity(), openMVG::Vec3::Zero()), 0);
  //displayRelativePoseReprojection(result, iTracker);
}

#if 0

// Display reprojection error based on a relative pose
void Rig::displayRelativePoseReprojection(const geometry::Pose3 & relativePose, std::size_t iTracker)
{
#ifdef VISUAL_DEBUG_MODE
  PonctualMarkerTracker & mainTracker = _vLocalizationResults[0];
  PonctualMarkerTracker & tracker = _vLocalizationResults[iTracker];

  // Set the marker size
  std::size_t semiWidth = 3.0;
  
  for(int j=0 ; j < tracker.poses().size() ; ++j )
  {
    // Window to display reprojection errors
    cv::Mat imgRes(tracker._height, tracker._width, CV_8UC3);// todo
    imgRes = cv::Scalar(255,255,255);
    cvNamedWindow("Reprojection", CV_WINDOW_NORMAL);
    cv::moveWindow("Reprojection",0,0);
    cv::resizeWindow("Reprojection", tracker._width/2, tracker._height/2);
    cv::imshow("Reprojection", imgRes);

    // Check that both pose computations succeed 
    if ( tracker.poses()[j].second )
    {
      //geometry::Pose3 pose = mainTracker.poses()[j].first*relativePose;
      
      const openMVG::Mat3 R1 = mainTracker.poses()[j].first.rotation();
      const openMVG::Vec3 t1 = mainTracker.poses()[j].first.translation();

      const openMVG::Mat3 R12 = relativePose.rotation();
      const openMVG::Vec3 t12 = relativePose.translation();

      const openMVG::Mat3 R2 = R12 * R1;
      const openMVG::Vec3 t2 = R12 * t1 + t12 ;

      const geometry::Pose3 pose( R2 , -R2.transpose() * t2 );
        
      //
      for(int k=0 ; k < tracker.pts()[j].size() ; ++k)
      {
        // Reprojections
        openMVG::Vec2 calibReprojPt;
        reproject( toOMVG(tracker.intrinsics()[j].getIntrinsics().getK()), pose, tracker.pts()[j][k], calibReprojPt);
        
        // Observations
        openMVG::Vec2 calibimgPt;
        popart::vision::projectiveTransform(toOMVG(tracker.intrinsics()[j].getIntrinsics().getK()), tracker.imgPts()[j][k], calibimgPt);
        
        // Display reprojections and observations
        cv::rectangle(imgRes, 
                cvPoint(calibimgPt(0)-semiWidth,calibimgPt(1)-semiWidth),
                cvPoint(calibimgPt(0)+semiWidth,calibimgPt(1)+semiWidth),
                cv::Scalar(0,255,0));
        
        cv::rectangle(imgRes,
                cvPoint(calibReprojPt(0)-semiWidth,calibReprojPt(1)-semiWidth),
                cvPoint(calibReprojPt(0)+semiWidth,calibReprojPt(1)+semiWidth),
                cv::Scalar(255,0,0));
      }
      cv::imshow("Reprojection", imgRes);
      cvpause();
    }
  }
#endif
}

bool Rig::optimizeCalibration()
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------

  //google::InitGoogleLogging("/home/lilian/toto.txt");
  
  ceres::Problem problem;

  // Add relative pose as a parameter block for each witness cameras (i.e. each Tracker)
  std::vector<std::vector<double> > vRelativePoses;
  for (int iRelativePose = 0 ; iRelativePose < _vRelativePoses.size() ; ++iRelativePose )
  {
    geometry::Pose3 & pose = _vRelativePoses[iRelativePose];
    
    const openMVG::Mat3 R = pose.rotation();
    const openMVG::Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    
    std::vector<double> relativePose;
    relativePose.reserve(6); //angleAxis + translation
    relativePose.push_back(angleAxis[0]);
    relativePose.push_back(angleAxis[1]);
    relativePose.push_back(angleAxis[2]);
    relativePose.push_back(t(0));
    relativePose.push_back(t(1));
    relativePose.push_back(t(2));

    vRelativePoses.push_back(relativePose);
    
    double * parameter_block = &vRelativePoses.back()[0];
    problem.AddParameterBlock(parameter_block, 6);
  }
  
  
  std::vector<std::vector<double> > vMainPoses;
  for (int iView = 0 ; iView < _vLocalizationResults[0].poses().size() ; ++iView )
  {
    if ( _vLocalizationResults[0].poses()[iView].second )
    {
      geometry::Pose3 & pose = _vLocalizationResults[0].poses()[iView].first;

      POPART_COUT_VAR(pose.rotation());
      const openMVG::Mat3 R = pose.rotation();
      const openMVG::Vec3 t = pose.translation();

      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);

      std::vector<double> mainPose;
      mainPose.reserve(6); //angleAxis + translation
      mainPose.push_back(angleAxis[0]);
      mainPose.push_back(angleAxis[1]);
      mainPose.push_back(angleAxis[2]);
      mainPose.push_back(t(0));
      mainPose.push_back(t(1));
      mainPose.push_back(t(2));

      vMainPoses.push_back(mainPose);
      
      double * parameter_block = &vMainPoses.back()[0];
      problem.AddParameterBlock(parameter_block, 6);
    }
  }
  
  POPART_COUT("Second init loop");

// The following code will be used if the intrinsics have to be refined in the bundle adjustment
  
//  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
//    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
//  {
//    const IndexT indexCam = itIntrinsic->first;
//
//    switch (itIntrinsic->second->getType())
//    {
//      case PINHOLE_CAMERA:
//      {
//        std::vector<double> vec_params = itIntrinsic->second->getParams();
//        map_intrinsics[indexCam].reserve(3);
//        map_intrinsics[indexCam].push_back(vec_params[0]);
//        map_intrinsics[indexCam].push_back(vec_params[1]);
//        map_intrinsics[indexCam].push_back(vec_params[2]);
//
//        double * parameter_block = &map_intrinsics[indexCam][0];
//        problem.AddParameterBlock(parameter_block, 3);
//        if (!bRefineIntrinsics)
//        {
//          //set the whole parameter block as constant for best performance.
//          problem.SetParameterBlockConstant(parameter_block);
//        }
//      }
//      break;
//    }
//  }


//#if 0
   // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = NULL;//new ceres::HuberLoss(Square(4.0));
  // TODO: make the LOSS function and the parameter an option

  // For all visibility add reprojections errors:
  for (int iTracker = 0 ; iTracker < _vLocalizationResults.size() ; ++iTracker)
  {
    const PonctualMarkerTracker & tracker = _vLocalizationResults[iTracker];
    
    for (int iView = 0 ; iView < tracker.poses().size() ; ++iView)
    {
      const TrackingCamera & cameraInfo = tracker.intrinsics()[iView];
      
      const std::vector<openMVG::Vec3> & points = tracker.pts()[iView];
      const std::vector<openMVG::Vec2> & observations = tracker.imgPts()[iView];
      
      for (int iPoint = 0 ; iPoint < points.size() ; ++iPoint)
      {
        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observations
        // and the 3D point and compares the reprojection against the observation.
        ceres::CostFunction* cost_function;

        // Add a residual block for the main camera
        if ( iTracker == 0 )
        {
          // Add the residual block if the resection (of the main camera) succeeded
          if ( _vLocalizationResults[iTracker].poses()[iView].second )
          {
            // Vector-2 residual, pose of the rig parameterized by 6 parameters
            cost_function = new ceres::AutoDiffCostFunction<ResidualErrorMainCameraFunctor, 2, 6>(
            new ResidualErrorMainCameraFunctor(toOMVG(cameraInfo.getIntrinsics().getK()), observations[iPoint], points[iPoint]));

            if (cost_function)
            {
              problem.AddResidualBlock( cost_function,
                                        p_LossFunction,
                                        &vMainPoses[iView][0]);
            }else
            {
              POPART_COUT("Fail in adding residual block for the main camera");
            }
          }

        }else
        // Add a residual block for a secondary camera
        {
          // Add the residual block if the resection (of the secondary camera) succeeded
          if ( _vLocalizationResults[iTracker].poses()[iView].second )
          {
                    //POPART_COUT_VAR(observations[iPoint]);
                    //POPART_COUT_VAR(points[iPoint]);
                   // POPART_COUT_VAR(toOMVG(cameraInfo.getIntrinsics().getK()));
            // Vector-2 residual, pose of the rig parametrised by 6 parameters
            //                  + relative pose of the secondary camera parameterized by 6 parameters
            cost_function = new ceres::AutoDiffCostFunction<ResidualErrorSecondaryCameraFunctor, 2, 6, 6>(
            new ResidualErrorSecondaryCameraFunctor(toOMVG(cameraInfo.getIntrinsics().getK()), observations[iPoint], points[iPoint]));

            if (cost_function)
            {
              //POPART_COUT("Second: Add residual block");
              problem.AddResidualBlock( cost_function,
                                        p_LossFunction,
                                        &vMainPoses[iView][0],
                                        &vRelativePoses[iTracker-1][0]);
            }else
            {
              POPART_COUT("Fail in adding residual block for a secondary camera");
            }
          }
        }
      }
    }
  }

  // Configure a BA engine and run it
  // todo @L set the most suitable options
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
  
#if 1
  if (openMVG_options._bCeres_Summary)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (openMVG_options._bVerbose)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (openMVG_options._bVerbose)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
      //  << " #views: " << sfm_data.views.size() << "\n"
      //  << " #poses: " << sfm_data.poses.size() << "\n"
      //  << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
      //  << " #tracks: " << sfm_data.structure.size() << "\n"        // Set with the correct sfm data derived from the rig
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << std::endl;
    }
    
    // Update relative pose after optimization
    for (int iRelativePose = 0 ; iRelativePose < _vRelativePoses.size() ; ++iRelativePose )
    {
      openMVG::Mat3 R_refined;
      std::vector<double> vPose;
      vPose.reserve(6);
      ceres::AngleAxisToRotationMatrix(&vRelativePoses[iRelativePose][0], R_refined.data());
      openMVG::Vec3 t_refined(vRelativePoses[iRelativePose][3], vRelativePoses[iRelativePose][4], vRelativePoses[iRelativePose][5]);
      // Update the pose
      geometry::Pose3 & pose = _vRelativePoses[iRelativePose];
      pose = geometry::Pose3(R_refined, -R_refined.transpose() * t_refined);
    }
    
    // Update relative poses along with rig poses after optimization
    for (int iPose = 0 ; iPose < _vLocalizationResults[0].poses().size() ; ++iPose )
    {
      if( _vLocalizationResults[0].poses()[iPose].second )
      {
        openMVG::Mat3 R_refined;
        std::vector<double> vPose;
        vPose.reserve(6);
        ceres::AngleAxisToRotationMatrix(&vMainPoses[iPose][0], R_refined.data());
        openMVG::Vec3 t_refined(vMainPoses[iPose][3], vMainPoses[iPose][4], vMainPoses[iPose][5]);
        // Push the optimized pose
        geometry::Pose3 pose = geometry::Pose3(R_refined, -R_refined.transpose() * t_refined);
        _vLocalizationResults[0].poses()[iPose].first = pose;
        _vPoses.push_back(pose);
      }else{
        // todo@L garbage...
        _vPoses.push_back(geometry::Pose3());
      }
    }
    
    // Update all poses in all trackers
    for (int iRelativePose = 0 ; iRelativePose < _vRelativePoses.size() ; ++iRelativePose )
    {
      std::size_t iTracker = iRelativePose+1;
      for (int iPose = 0 ; iPose < _vLocalizationResults[iTracker].poses().size() ; ++iPose )
      {
        if( _vLocalizationResults[iTracker].poses()[iPose].second && _vLocalizationResults[0].poses()[iPose].second )
        {
          const geometry::Pose3 & relativePose = _vRelativePoses[iRelativePose];
          
          const openMVG::Mat3 R1 = _vLocalizationResults[0].poses()[iPose].first.rotation();
          const openMVG::Vec3 t1 = _vLocalizationResults[0].poses()[iPose].first.translation();

          const openMVG::Mat3 R12 = relativePose.rotation();
          const openMVG::Vec3 t12 = relativePose.translation();

          const openMVG::Mat3 R2 = R12 * R1;
          const openMVG::Vec3 t2 = R12 * t1 + t12 ;
          
          _vLocalizationResults[iTracker].poses()[iPose].first = geometry::Pose3( R2 , -R2.transpose() * t2 );
          
          
        }
      }
    }
    
    displayRelativePoseReprojection(geometry::Pose3(openMVG::Mat3::Identity(), openMVG::Vec3::Zero()), 0);
    displayRelativePoseReprojection(_vRelativePoses[0], 1);
// Possibility to update the intrinsics here
    
// Update camera intrinsics with refined data
//    for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
//      itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
//    {
//      const IndexT indexCam = itIntrinsic->first;
//
//      const std::vector<double> & vec_params = map_intrinsics[indexCam];
//      itIntrinsic->second.get()->updateFromParams(vec_params);
//    }
    return true;
  }
  #endif
}

// Compute an average pose
void poseAveraging(const std::vector<geometry::Pose3> & vPoses, geometry::Pose3 & result)
{
  // todo
}

double reprojectionError(const openMVG::Mat3 & K, const geometry::Pose3 & pose, const std::vector<openMVG::Vec2> & imgPts, const std::vector<openMVG::Vec3> & pts)
{
  double res = 0;
  
  //assert( imgPts.size - pts.size() == 0 );
  
  for(int i=0 ; i < imgPts.size() ; ++i)
  {
    res += residual(K, pose, pts[i], imgPts[i]);
  }
  return res;
}

// Compute the residual between the 3D projected point X and an image observation x
double residual( const openMVG::Mat3 & K, const geometry::Pose3 & pose, const openMVG::Vec3 & X, const openMVG::Vec2 & x)
{
  openMVG::Vec2 reproj;
  reproject( K, pose, X, reproj );
  openMVG::Vec2 calibImgPt;
  projectiveTransform(K, x, calibImgPt);
  double diffX, diffY;
  diffX = reproj(0) - calibImgPt(0);
  diffY = reproj(1) - calibImgPt(1);
  return ( diffX*diffX + diffY*diffY );
}

// Compute the reprojection of X through the camera whose pose pose
void reproject( const openMVG::Mat3 & K, const geometry::Pose3 & pose, const openMVG::Vec3 & X, openMVG::Vec2 & x)
{
  openMVG::Vec3 proj = K * pose(X);
  x(0) = proj(0)/proj(2);
  x(1) = proj(1)/proj(2);
}

// Apply a 2D projective transformation
void projectiveTransform(const openMVG::Mat3 & H, const openMVG::Vec2 & x, openMVG::Vec2 & res)
{
  double u = H(0,0)*x(0)+H(0,1)*x(1)+H(0,2);
  double v = H(1,0)*x(0)+H(1,1)*x(1)+H(1,2);
  double w = H(2,0)*x(0)+H(2,1)*x(1)+H(2,2);
  res(0) = u/w;
  res(1) = v/w;
}
#endif

} // namespace rig
} // namespace openMVG
