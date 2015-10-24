#include "Rig.hpp"
#include "bundleAdjustmentCeresFunctor.hpp"

//#include <vision/cameraTracking/debug/visualDebug.hpp>
//#include <eigen3/Eigen/src/Core/MatrixBase.h>

#include <openMVG/sfm/sfm_data_BA_ceres.hpp>

#include <ceres/rotation.h>
#include <opencv2/opencv.hpp>

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
  for (int iLocalizer=1 ; iLocalizer < nCams ; ++iLocalizer)
  {
    // Perform the pose averaging over all relative pose between the main camera
    // (index 0) and the witness camera (index i)
    std::vector<localization::LocalizationResult> & resWitnessCamera = _vLocalizationResults[iLocalizer];
    
    // vRelativePoses will store all the relative poses overall frames where both
    // the pose computation of the main camera and witness camera succeed
    std::vector<geometry::Pose3> vRelativePoses;
    vRelativePoses.reserve(resWitnessCamera.size());
    
    for(int iView=0 ; iView < resWitnessCamera.size() ; ++iView )
    {
      // Check that both pose computations succeed 
      if ( resMainCamera[iView].isValid() && resWitnessCamera[iView].isValid() )
      {
        vRelativePoses.push_back(
                computeRelativePose(resMainCamera[iView].getPose(), resWitnessCamera[iView].getPose())
                );
      }
    }
    geometry::Pose3 optimalRelativePose;
    findBestRelativePose(vRelativePoses, iLocalizer, optimalRelativePose );
  
    //poseAveraging(vRelativePoses, averageRelativePose);
    _vRelativePoses.push_back(optimalRelativePose);
  }
  
  // Update all poses in all localization results
  for (int iRelativePose = 0 ; iRelativePose < _vRelativePoses.size() ; ++iRelativePose )
  {
    std::size_t iRes = iRelativePose+1;
    for (int iView = 0 ; iView < _vLocalizationResults[iRes].size() ; ++iView )
    {
      if(  _vLocalizationResults[0][iView].isValid() )
      {
        const geometry::Pose3 & relativePose = _vRelativePoses[iRelativePose];
        
        const geometry::Pose3 poseWitnessCamera = poseFromMainToWitness(_vLocalizationResults[0][iView].getPose(), relativePose);
        _vLocalizationResults[iRes][iView].setPose(poseWitnessCamera);
        
      }
    }
  }
}

// From a set of relative pose, find the optimal one for a given tracker iTraker which
// minimize the reprojection errors over all images
void Rig::findBestRelativePose(
        const std::vector<geometry::Pose3> & vPoses,
        std::size_t iLocalizer,
        geometry::Pose3 & result )
{
  std::vector<localization::LocalizationResult> & resMainCamera = _vLocalizationResults[0];
  std::vector<localization::LocalizationResult> & resWitnessCamera = _vLocalizationResults[iLocalizer];
  
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
        const geometry::Pose3 poseWitnessCamera = poseFromMainToWitness(resMainCamera[j].getPose(), relativePose);
        error += reprojectionError(resWitnessCamera[j], poseWitnessCamera);
      }
    }
    if ( error < minReprojError )
    {
      iMin = i;
      minReprojError = error;
    }           
  }
  result = vPoses[iMin];
  
  displayRelativePoseReprojection(geometry::Pose3(openMVG::Mat3::Identity(), openMVG::Vec3::Zero()), 0);
  displayRelativePoseReprojection(result, iLocalizer);
}

geometry::Pose3 computeRelativePose(geometry::Pose3 poseMainCamera, geometry::Pose3 poseWitnessCamera)
{
  const openMVG::Mat3 & R1 = poseMainCamera.rotation();
  const openMVG::Vec3 & t1 = poseMainCamera.translation();
  const openMVG::Mat3 & R2 = poseWitnessCamera.rotation();
  const openMVG::Vec3 & t2 = poseWitnessCamera.translation();

  const openMVG::Mat3 R12 = R2 * R1.transpose();
  const openMVG::Vec3 t12 = t2 - R12 * t1;

  return geometry::Pose3( R12 , -R12.transpose()*t12 );
}
        

geometry::Pose3 poseFromMainToWitness(geometry::Pose3 poseMainCamera, geometry::Pose3 relativePose)
{
  const openMVG::Mat3 & R1 = poseMainCamera.rotation();
  const openMVG::Vec3 & t1 = poseMainCamera.translation();

  const openMVG::Mat3 & R12 = relativePose.rotation();
  const openMVG::Vec3 & t12 = relativePose.translation();

  const openMVG::Mat3 R2 = R12 * R1;
  const openMVG::Vec3 t2 = R12 * t1 + t12 ;
  
  return geometry::Pose3( R2 , -R2.transpose() * t2 );
}

double reprojectionError(const localization::LocalizationResult & localizationResult, const geometry::Pose3 & pose)
{
  double residual = 0;
  
  for(const IndexT iInliers : localizationResult.getMatchData().vec_inliers)
  {
    // Inlier 3D point
    const Vec3 & point3D = localizationResult.getMatchData().pt3D.col(iInliers);
    // Its reprojection
    Vec2 itsReprojection = localizationResult.getIntrinsics().project(pose, point3D);
    // Its associated observation location
    const Vec2 & point2D = localizationResult.getMatchData().pt2D.col(iInliers);
    // Residual
    residual += (point2D(0) - itsReprojection(0))*(point2D(0) - itsReprojection(0));
    residual += (point2D(1) - itsReprojection(1))*(point2D(1) - itsReprojection(1));
  }
  return residual;
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

// Display reprojection error based on a relative pose
void Rig::displayRelativePoseReprojection(const geometry::Pose3 & relativePose, std::size_t iLocalizer)
{
#ifdef VISUAL_DEBUG_MODE
  std::vector<localization::LocalizationResult> & mainLocalizerResults = _vLocalizationResults[0];
  std::vector<localization::LocalizationResult> & witnessLocalizerResults = _vLocalizationResults[iLocalizer];

  // Set the marker size
  std::size_t semiWidth = 3.0;
  
  for(int iView=0 ; iView < witnessLocalizerResults.size() ; ++iView )
  {
        // Check that both pose computations succeed 
    if ( witnessLocalizerResults[iView].isValid() )
    {
      const std::size_t width = witnessLocalizerResults[iView].getIntrinsics()._w;
      const std::size_t height = witnessLocalizerResults[iView].getIntrinsics()._h;

      // Window to display reprojection errors
      cv::Mat imgRes(height, width, CV_8UC3);
      imgRes = cv::Scalar(255,255,255);
      cvNamedWindow("Reprojection", CV_WINDOW_NORMAL);
      cv::moveWindow("Reprojection",0,0);
      cv::resizeWindow("Reprojection", width/2, height/2);
      cv::imshow("Reprojection", imgRes);
      
      const geometry::Pose3 poseWitnessCamera = poseFromMainToWitness(mainLocalizerResults[iView].getPose(), relativePose);
        
      const openMVG::Mat points2D = witnessLocalizerResults[iView].getMatchData().pt2D;
      const openMVG::Mat points3D = witnessLocalizerResults[iView].getMatchData().pt3D;
      
      for(const IndexT iInlier : witnessLocalizerResults[iView].getMatchData().vec_inliers)
      {
        // Reprojections
        const Vec3 & point3D = points3D.col(iInlier);
        // Its reprojection
        Vec2 itsReprojection = witnessLocalizerResults[iView].getIntrinsics().project(poseWitnessCamera, point3D);
        // Its associated observation location
        const Vec2 & point2D = points2D.col(iInlier);
 
        // Display reprojections and observations
        cv::rectangle(imgRes, 
                cvPoint(point2D(0)-semiWidth,point2D(1)-semiWidth),
                cvPoint(point2D(0)+semiWidth,point2D(1)+semiWidth),
                cv::Scalar(0,255,0));
        
        cv::rectangle(imgRes,
                cvPoint(itsReprojection(0)-semiWidth,itsReprojection(1)-semiWidth),
                cvPoint(itsReprojection(0)+semiWidth,itsReprojection(1)+semiWidth),
                cv::Scalar(255,0,0));
      }
      cv::imshow("Reprojection", imgRes);
      cvpause();
    }
  }
#endif
}

/* Pause function while using cv::namedwindows*/
void cvpause(){
#ifdef VISUAL_DEBUG_MODE
  int keyboard;
  while( !(char) keyboard ){ // ASCII code for 'CR'
    keyboard = cv::waitKey( 0 );
  }

  if ( (char) keyboard == 'q' )
  {
    std::cerr << "The program has been manually stopped" << std::endl;
    std::exit(0);
  }
#endif
}

bool Rig::optimizeCalibration()
{
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
  for (int iView = 0 ; iView < _vLocalizationResults[0].size() ; ++iView )
  {
    if ( _vLocalizationResults[0][iView].isValid() )
    {
      const geometry::Pose3 & pose = _vLocalizationResults[0][iView].getPose();

      //POPART_COUT_VAR(pose.rotation());
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

// The following code can be used if the intrinsics have to be refined in the bundle adjustment
  
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

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = NULL;//new ceres::HuberLoss(Square(4.0));
  // TODO: make the LOSS function and the parameter an option

  // For all visibility add reprojections errors:
  for (int iLocalizer = 0 ; iLocalizer < _vLocalizationResults.size() ; ++iLocalizer)
  {
    const std::vector<localization::LocalizationResult> & localizationResults = _vLocalizationResults[iLocalizer];
    
    for (int iView = 0 ; iView < localizationResults.size() ; ++iView)
    {
      const cameras::Pinhole_Intrinsic & cameraInfo = localizationResults[iView].getIntrinsics();
      
      // Get the inliers 3D points
      const Mat & points3D = localizationResults[iView].getMatchData().pt3D;
      // Get their image locations
      const Mat & points2D = localizationResults[iView].getMatchData().pt2D;
      
      // Add a residual block for all inliers
      for (const IndexT iPoint : localizationResults[iView].getMatchData().vec_inliers )
      {
        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observations
        // and the 3D point and compares the reprojection against the observation.
        ceres::CostFunction* cost_function;

        // Add a residual block for the main camera
        if ( iLocalizer == 0 )
        {
          // Add the residual block if the resection (of the main camera) succeeded
          if ( _vLocalizationResults[iLocalizer][iView].isValid() )
          {
            // Vector-2 residual, pose of the rig parameterized by 6 parameters
            cost_function = new ceres::AutoDiffCostFunction<ResidualErrorMainCameraFunctor, 2, 6>(
            new ResidualErrorMainCameraFunctor(_vLocalizationResults[iLocalizer][iView].getIntrinsics().K(), points2D.col(iPoint), points3D.col(iPoint) ));

            if (cost_function)
            {
              problem.AddResidualBlock( cost_function,
                                        p_LossFunction,
                                        &vMainPoses[iView][0]);
            }else
            {
              std::cout << "Fail in adding residual block for the main camera" << std::endl;
            }
          }

        }else
        // Add a residual block for a secondary camera
        {
          // Add the residual block if the resection (of the secondary camera) succeeded
          if ( _vLocalizationResults[iLocalizer][iView].isValid() )
          {
                    //POPART_COUT_VAR(observations[iPoint]);
                    //POPART_COUT_VAR(points[iPoint]);
                   // POPART_COUT_VAR(toOMVG(cameraInfo.getIntrinsics().getK()));
            // Vector-2 residual, pose of the rig parametrised by 6 parameters
            //                  + relative pose of the secondary camera parameterized by 6 parameters
            cost_function = new ceres::AutoDiffCostFunction<ResidualErrorSecondaryCameraFunctor, 2, 6, 6>(
            new ResidualErrorSecondaryCameraFunctor(_vLocalizationResults[iLocalizer][iView].getIntrinsics().K(), points2D.col(iPoint), points3D.col(iPoint)));

            if (cost_function)
            {
              //POPART_COUT("Second: Add residual block");
              problem.AddResidualBlock( cost_function,
                                        p_LossFunction,
                                        &vMainPoses[iView][0],
                                        &vRelativePoses[iLocalizer-1][0]);
            }else
            {
              std::cout << "Fail in adding residual block for a secondary camera" << std::endl;
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
    
    // Update the main camera pose after optimization
    for (int iView = 0 ; iView < _vLocalizationResults[0].size() ; ++iView )
    {
      if( _vLocalizationResults[0][iView].isValid() )
      {
        openMVG::Mat3 R_refined;
        std::vector<double> vPose;
        vPose.reserve(6);
        ceres::AngleAxisToRotationMatrix(&vMainPoses[iView][0], R_refined.data());
        openMVG::Vec3 t_refined(vMainPoses[iView][3], vMainPoses[iView][4], vMainPoses[iView][5]);
        // Push the optimized pose
        geometry::Pose3 pose = geometry::Pose3(R_refined, -R_refined.transpose() * t_refined);
        _vLocalizationResults[0][iView].setPose(pose);
        _vPoses.push_back(pose);
      }
    }
    
    // Update all poses over all witness cameras after optimization
    for (int iRelativePose = 0 ; iRelativePose < _vRelativePoses.size() ; ++iRelativePose )
    {
      std::size_t iLocalizer = iRelativePose+1;
      for (int iView = 0 ; iView < _vLocalizationResults[iLocalizer].size() ; ++iView )
      {
        if( _vLocalizationResults[iLocalizer][iView].isValid() && _vLocalizationResults[0][iView].isValid() )
        {
          // Retrieve the witness camera pose from the main camera one.
          const geometry::Pose3 poseWitnessCamera = poseFromMainToWitness(_vLocalizationResults[0][iView].getPose(), _vRelativePoses[iRelativePose]);
          _vLocalizationResults[iLocalizer][iView].setPose(poseWitnessCamera);
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
}

#if 0

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
