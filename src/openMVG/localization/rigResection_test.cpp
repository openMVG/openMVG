/* 
 * File:   rigResection_test.cpp
 * Author: sgaspari
 *
 * Created on January 2, 2016, 11:38 PM
 */

#include "rigResection.hpp"
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>
#include <openMVG/geometry/pose3.hpp>

#include "testing/testing.h"
#include "optimization.hpp"

#include <vector>
#include <math.h>
#include <chrono>
#include <random>
#include <random>

using namespace openMVG;

Mat3 generateRotation(double x, double y, double z)
{
//std::cout << "generateRotation" << std::endl; 
  Mat3 R1 = Mat3::Identity();
  R1(1,1) = cos(x);
  R1(1,2) = -sin(x);
  R1(2,1) = -R1(1,2);
  R1(2,2) = R1(1,1);

  Mat3 R2 = Mat3::Identity();
  R2(0,0) = cos(y);
  R2(0,2) = sin(y);
  R2(2,0) = -R2(0,2);
  R2(2,2) = R2(0,0);

  Mat3 R3 = Mat3::Identity();
  R3(0,0) = cos(z);
  R3(0,1) = -sin(z);
  R3(1,0) =-R3(0,1);
  R3(1,1) = R3(0,0);
  return R3 * R2 * R1;
}

Mat3 generateRotation(const Vec3 &angles)
{
  return generateRotation(angles(0), angles(1), angles(2));
}

Mat3 generateRandomRotation(const Vec3 &maxAngles = Vec3::Constant(2*M_PI))
{
  const Vec3 angles = Vec3::Random().cwiseProduct(maxAngles);
//  std::cout << "generateRandomRotation" << std::endl; 
  return generateRotation(angles);
}

Vec3 generateRandomTranslation(double maxNorm)
{
//  std::cout << "generateRandomTranslation" << std::endl; 
  Vec3 translation = Vec3::Random();
  return maxNorm * (translation / translation.norm());
}

Vec3 generateRandomPoint(double thetaMin, double thetaMax, double depthMax, double depthMin)
{
  // variables contains the spherical parameters (r, theta, phi) for the point to generate its 
  // spherical coordinatrs
  Vec3 variables = Vec3::Random();
  variables(0) = (variables(0)+1)/2;    // put the random in [0,1]
  variables(1) = (variables(1)+1)/2;    // put the random in [0,1]
   
  // get random value for r
  const double &r = variables(0) = depthMin*variables(0) + (1-variables(0))*depthMax;
  // get random value for theta
  const double &theta = variables(1) = thetaMin*variables(1) + (1-variables(1))*thetaMax;
  // get random value for phi
  const double &phi = variables(2) = M_PI*variables(2);
  
  return Vec3(r*std::sin(theta)*std::cos(phi), 
              r*std::sin(theta)*std::sin(phi), 
              r*std::cos(theta));
}

Mat3X generateRandomPoints(std::size_t numPts, double thetaMin, double thetaMax, double depthMax, double depthMin)
{
  Mat3X points = Mat(3, numPts);
  for(std::size_t i = 0; i < numPts; ++i)
  {
    points.col(i) = generateRandomPoint(thetaMin, thetaMax, depthMax, depthMin);
  }
  return points;
}

geometry::Pose3 generateRandomPose(const Vec3 &maxAngles = Vec3::Constant(2*M_PI), double maxNorm = 1)
{
//  std::cout << "generateRandomPose" << std::endl; 
  return geometry::Pose3(generateRandomRotation(maxAngles), generateRandomTranslation(maxNorm));
}

void generateRandomExperiment(std::size_t numCameras, 
                              std::size_t numPoints,
                              double outliersPercentage,
                              double noise,
                              geometry::Pose3 &rigPoseGT,
                              Mat3X &pointsGT,
                              std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                              std::vector<geometry::Pose3 > &vec_subPoses,
                              std::vector<Mat> &vec_pts3d,
                              std::vector<Mat> &vec_pts2d)
{
      // generate random pose for the rig
    rigPoseGT = generateRandomPose(Vec3::Constant(M_PI/10), 5);

    // generate random 3D points
    const double thetaMin = 0;
    const double thetaMax = M_PI/3;
    const double depthMax = 50;
    const double depthMin = 10;
    pointsGT = generateRandomPoints(numPoints, thetaMin, thetaMax, depthMax, depthMin);

    // apply random rig pose to points
    const Mat3X points = rigPoseGT(pointsGT);
    std::cout << "rigPoseGT\n" << rigPoseGT.rotation() << "\n" << rigPoseGT.center()<< std::endl;

    // generate numCameras random poses and intrinsics
    for(std::size_t cam = 0; cam < numCameras; ++cam)
    {
      // first camera is in I 0
      if(cam != 0)
        vec_subPoses.push_back(generateRandomPose(Vec3::Constant(M_PI/10), 1.5));

      // let's keep it simple
      vec_queryIntrinsics.push_back(cameras::Pinhole_Intrinsic_Radial_K3(640, 480, 500, 320, 240));
    }
    assert(vec_subPoses.size() == numCameras-1);
  //  for(std::size_t i = 0; i < vec_subPoses.size(); ++i)
  //    std::cout << "vec_subPoses\n" << vec_subPoses[i].rotation() << "\n" << vec_subPoses[i].center()<< std::endl;;

    // for each camera generate the features (if 3D point is "in front" of the camera)
    vec_pts3d.reserve(numCameras);
    vec_pts2d.reserve(numCameras);
    
    const std::size_t numOutliers = (std::size_t)numPoints*outliersPercentage;

    for(std::size_t cam = 0; cam < numCameras; ++cam)
    {
      Mat3X localPts;
      if(cam != 0)
        localPts = vec_subPoses[cam-1](points);
      else
        localPts = points;

      // count the number of points in front of the camera
      // ie take the 3rd coordinate and check it is > 0
      const std::size_t validPoints = (localPts.row(2).array() > 0).count();
      Mat3X pts3d = Mat(3, validPoints + numOutliers);
      Mat2X pts2d = Mat(2, validPoints + numOutliers);

      // for each 3D point
      std::size_t idx = 0;
      for(std::size_t i = 0; i < numPoints; ++i)
      {
        // if it is in front of the camera
        if(localPts(2,i) > 0)
        {
          // project it
          Vec2 feat = vec_queryIntrinsics[cam].project(geometry::Pose3(), localPts.col(i));
          
          if(noise > 0.0)
          {
            feat = feat + noise*Vec2::Random();
          }

          // add the 3d and 2d point
          pts3d.col(idx) = pointsGT.col(i);
          pts2d.col(idx) = feat;
          ++idx;
        }
      }
      assert(idx == validPoints);
      
      if(numOutliers)
      {
        //add some other random associations
        for(std::size_t i = 0; i < numOutliers; ++i)
        {
          pts3d.col(idx) = 10*Vec3::Random();
          pts2d.col(idx) = 100*Vec2::Random();
          ++idx;
        }
      }

//      std::cout << "Cam " << cam << std::endl;
//      std::cout << "pts2d\n" << pts2d << std::endl;
//      std::cout << "pts3d\n" << pts3d << std::endl;
//      std::cout << "pts3dTRA\n" << localPts << std::endl;
//      
//      auto residuals = vec_queryIntrinsics[cam].residuals(geometry::Pose3(), localPts, pts2d);
//      auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
//      
//      std::cout << "residuals\n" << sqrErrors << std::endl;

  //    if(cam!=0)
  //      residuals = vec_queryIntrinsics[cam].residuals(vec_subPoses[cam-1], points, pts2d);
  //    else
  //      residuals = vec_queryIntrinsics[cam].residuals(geometry::Pose3(), points, pts2d);
  //    
  //    auto sqrErrors2 = (residuals.cwiseProduct(residuals)).colwise().sum();
  //    
  //    std::cout << "residuals2\n" << sqrErrors2 << std::endl;

      vec_pts3d.push_back(pts3d);
      vec_pts2d.push_back(pts2d);
    }
}


TEST(rigResection, simpleNoNoiseNoOutliers)
{
  const std::size_t numCameras = 3;
  const std::size_t numPoints = 10;
  const std::size_t numTrials = 10;
  const double threshold = 1e-3;
  
  for(std::size_t trial = 0; trial < numTrials; ++trial)
  {

    // generate random pose for the rig
    geometry::Pose3 rigPoseGT;
    std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > vec_queryIntrinsics;
    std::vector<geometry::Pose3 > vec_subPoses;
    std::vector<Mat> vec_pts3d;
    std::vector<Mat> vec_pts2d;
    Mat3X pointsGT;

    generateRandomExperiment(numCameras, 
                             numPoints,
                             0,
                             0,
                             rigPoseGT,
                             pointsGT,
                             vec_queryIntrinsics,
                             vec_subPoses,
                             vec_pts3d,
                             vec_pts2d);

    // call the GPNP
    std::vector<std::vector<std::size_t> > inliers;
    geometry::Pose3 rigPose;
    EXPECT_TRUE(localization::rigResection(vec_pts2d,
                                          vec_pts3d,
                                          vec_queryIntrinsics,
                                          vec_subPoses,
                                          rigPose,
                                          inliers));

    std::cout << "rigPose\n" << rigPose.rotation() << "\n" << rigPose.center()<< std::endl;

    // check result for the pose
    const auto poseDiff = rigPose*rigPoseGT.inverse();
    std::cout << "diffR\n" << poseDiff.rotation() << std::endl;
    std::cout << "diffC\n" << poseDiff.center().norm() << std::endl;
    
    // check result for the posediff, how close it is to the identity
    const Mat3 &rot = poseDiff.rotation();
    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
      {
        if(i == j)
        {
          EXPECT_NEAR(rot(i,j), 1.0, threshold);
        }
        else
        {
          EXPECT_NEAR(rot(i,j), 0.0, threshold);
        }
      }
    }

    const Vec3 &center = poseDiff.center();
    for(std::size_t i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(center(i), 0.0, 10*threshold);
    }

    // check inliers
    EXPECT_TRUE(inliers.size() == numCameras);
    for(std::size_t i = 0; i < numCameras; ++i)
    {
      EXPECT_TRUE(inliers[i].size() == vec_pts2d[i].cols());
    }

    EXPECT_TRUE(localization::refineRigPose(vec_pts2d,
                                            vec_pts3d,
                                            inliers,
                                            vec_queryIntrinsics,
                                            vec_subPoses,
                                            rigPose));

    // THIS TEST CAN FAIL AS THE REPROJECTION ERROR USED INSIDE THE GP3P IS BASED
    // ON THE ANGLE DIFFERENCE RATHER THAN THE REPROJECTED POINT DIFFERENCE
    // check reprojection errors
//    for(std::size_t cam = 0; cam < numCameras; ++cam)
//    {
//      const std::size_t numPts = vec_pts2d[cam].cols();
//      const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera = vec_queryIntrinsics[cam];
//      Mat2X residuals;
//      if(cam!=0)
//        residuals = currCamera.residuals(vec_subPoses[cam-1]*rigPose, vec_pts3d[cam], vec_pts2d[cam]);
//      else
//        residuals = currCamera.residuals(geometry::Pose3()*rigPose, vec_pts3d[cam], vec_pts2d[cam]);
//
//      auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
//      std::cout << sqrErrors << std::endl;
//      for(std::size_t j = 0; j < numPts; ++j)
//      {
//        EXPECT_TRUE(sqrErrors(j) <= threshold);
//      }
//    }
  }
}

TEST(rigResection, simpleWithNoiseNoOutliers)
{
  const std::size_t numCameras = 3;
  const std::size_t numPoints = 10;
  const std::size_t numTrials = 10;
  const double noise = 0.1;
  const double threshold = 0.1;
  
  for(std::size_t trial = 0; trial < numTrials; ++trial)
  {

    // generate random pose for the rig
    geometry::Pose3 rigPoseGT;
    std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > vec_queryIntrinsics;
    std::vector<geometry::Pose3 > vec_subPoses;
    std::vector<Mat> vec_pts3d;
    std::vector<Mat> vec_pts2d;
    Mat3X pointsGT;

    generateRandomExperiment(numCameras, 
                             numPoints,
                             0,
                             noise,
                             rigPoseGT,
                             pointsGT,
                             vec_queryIntrinsics,
                             vec_subPoses,
                             vec_pts3d,
                             vec_pts2d);

    // call the GPNP
    std::vector<std::vector<std::size_t> > inliers;
    geometry::Pose3 rigPose;
    EXPECT_TRUE(localization::rigResection(vec_pts2d,
                                           vec_pts3d,
                                           vec_queryIntrinsics,
                                           vec_subPoses,
                                           rigPose,
                                           inliers,
                                           D2R(0.008)));

    std::cout << "rigPose\n" << rigPose.rotation() << "\n" << rigPose.center()<< std::endl;
    
    const auto poseDiff = rigPose*rigPoseGT.inverse();
    std::cout << "diffR\n" << poseDiff.rotation() << std::endl;
    std::cout << "diffC\n" << poseDiff.center().norm() << std::endl;
    
    // check result for the posediff, how close it is to the identity
    const Mat3 &rot = poseDiff.rotation();
    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
      {
        if(i == j)
        {
          EXPECT_NEAR(rot(i,j), 1.0, threshold);
        }
        else
        {
          EXPECT_NEAR(rot(i,j), 0.0, threshold);
        }
      }
    }

    const Vec3 &center = poseDiff.center();
    for(std::size_t i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(center(i), 0.0, threshold);
    }

    EXPECT_TRUE(localization::refineRigPose(vec_pts2d,
                                            vec_pts3d,
                                            inliers,
                                            vec_queryIntrinsics,
                                            vec_subPoses,
                                            rigPose));
  }
}

TEST(rigResection, simpleNoNoiseWithOutliers)
{
  const std::size_t numCameras = 3;
  const std::size_t numPoints = 10;
  const std::size_t numTrials = 10;
  const double noise = 0;
  const double threshold = 1e-3;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 0.4);
  
  for(std::size_t trial = 0; trial < numTrials; ++trial)
  {

    // generate random pose for the rig
    geometry::Pose3 rigPoseGT;
    const double outlierPercentage = 0.4*dis(gen);
    std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > vec_queryIntrinsics;
    std::vector<geometry::Pose3 > vec_subPoses;
    std::vector<Mat> vec_pts3d;
    std::vector<Mat> vec_pts2d;
    Mat3X pointsGT;

    generateRandomExperiment(numCameras, 
                             numPoints,
                             outlierPercentage,
                             noise,
                             rigPoseGT,
                             pointsGT,
                             vec_queryIntrinsics,
                             vec_subPoses,
                             vec_pts3d,
                             vec_pts2d);

    // call the GPNP
    std::vector<std::vector<std::size_t> > inliers;
    geometry::Pose3 rigPose;
    EXPECT_TRUE(localization::rigResection(vec_pts2d,
                                          vec_pts3d,
                                          vec_queryIntrinsics,
                                          vec_subPoses,
                                          rigPose,
                                          inliers,
                                          D2R(0.008)));

    std::cout << "rigPose\n" << rigPose.rotation() << "\n" << rigPose.center()<< std::endl;

    EXPECT_TRUE(localization::refineRigPose(vec_pts2d,
                                            vec_pts3d,
                                            inliers,
                                            vec_queryIntrinsics,
                                            vec_subPoses,
                                            rigPose));

    // THIS TEST CAN FAIL AS THE REPROJECTION ERROR USED INSIDE THE GP3P IS BASED
    // ON THE ANGLE DIFFERENCE RATHER THAN THE REPROJECTED POINT DIFFERENCE
//    // check reprojection errors
//    for(std::size_t cam = 0; cam < numCameras; ++cam)
//    {
//      const std::size_t numPts = vec_pts2d[cam].cols();
//      const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera = vec_queryIntrinsics[cam];
//      Mat2X residuals;
//      if(cam!=0)
//        residuals = currCamera.residuals(vec_subPoses[cam-1]*rigPose, vec_pts3d[cam], vec_pts2d[cam]);
//      else
//        residuals = currCamera.residuals(geometry::Pose3()*rigPose, vec_pts3d[cam], vec_pts2d[cam]);
//  
//      auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
//
////      std::cout << sqrErrors << std::endl;
//
//      const auto &currInliers = inliers[cam];
//      for(std::size_t j = 0; j < currInliers.size(); ++j)
//      {
////        std::cout << sqrErrors(currInliers[j]) << std::endl;
//        EXPECT_TRUE(sqrErrors(currInliers[j]) <= threshold);
//      }
//    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

