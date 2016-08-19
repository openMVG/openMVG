/* 
 * File:   rigResection.cpp
 * Author: sgaspari
 * 
 * Created on January 2, 2016, 5:51 PM
 */

#include "rigResection.hpp"

#if HAVE_OPENGV
#include <Eigen/Eigen>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#endif

#include <stdio.h>
#include <iostream>
#include <vector>
#include <chrono>

namespace openMVG{
namespace localization{

#if HAVE_OPENGV

bool rigResection(const std::vector<Mat> &pts2d, 
                  const std::vector<Mat> &pts3d,
                  const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                  const std::vector<geometry::Pose3 > &vec_subPoses,
                  geometry::Pose3 &rigPose,
                  std::vector<std::vector<std::size_t> > &inliers,
                  double threshold /*= 1e-6*/,
                  size_t maxIterations /*= 100*/,
                  bool verbosity /*= true*/)
{
  const std::size_t numCameras = pts2d.size();
  
  assert(pts3d.size() == numCameras);
  assert(vec_queryIntrinsics.size() == numCameras);
  // the subposes are always wrt the main camera
  assert(vec_subPoses.size() == numCameras-1);
  
  
  // Copy the subposes into separate vectors for rotations and centers
  opengv::translations_t camOffsets(numCameras); //centers
  opengv::rotations_t camRotations(numCameras);
  
  // the first is the identity (the main camera)
  camOffsets[0] = Vec3::Zero();
  camRotations[0] = Mat3::Identity();
  // now copy the others
  for(std::size_t i = 1; i < numCameras; ++i)
  {
    camOffsets[i] = vec_subPoses[i-1].center();
    // it takes as input a real pose, ie the rotation wrt the world frame
    camRotations[i] = vec_subPoses[i-1].rotation().inverse();
  }
  
  // count the total number of associations
  std::size_t numTotalPoints = 0;
  for(std::size_t i = 0; i < numCameras; ++i)
  {
    const std::size_t numPts = pts2d[i].cols();
    assert(pts3d[i].cols() == numPts);
    assert(pts3d[i].rows() == 3);
    assert(pts2d[i].rows() == 2);
    numTotalPoints += numPts;
  }
  // transform the points into bearingVectors
  opengv::bearingVectors_t bearingVectors;  // this contains the bearing vector associated to each image feature
  opengv::points_t points;                  // this contains the 3d points
  points.reserve(numTotalPoints);
  std::vector<int> camCorrespondences;      // this is to keep track, for each bearing vector, whose camera it belongs
  
  // this map is used to remap the indices of the flatten structure containing
  // the points (points, bearingvector) to the original local index of each pts2d/pts3d: 
  // eg 
  // localIdx = absoluteToLocal[absoluteID]; such that
  // pts3d[localIdx.first].col(localIdx.second) == points[absoluteID];
  std::map<std::size_t, std::pair<std::size_t,std::size_t>> absoluteToLocal;   
  
  
  for(std::size_t cam = 0; cam < numCameras; ++cam)
  {
    const std::size_t numPts = pts2d[cam].cols();
    const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera = vec_queryIntrinsics[cam];
    
    for(std::size_t i = 0; i < numPts; ++i)
    {
      // store the 3D point
      points.push_back(pts3d[cam].col(i));
      
      // we first remove the distortion and then we transform the undistorted point in
      // normalized camera coordinates (inv(K)*undistortedPoint)
      auto pt = currCamera.ima2cam(currCamera.get_ud_pixel(pts2d[cam].col(i)));
      
      opengv::bearingVector_t bearing(pt(0), pt(1), 1.0);
      
      // normalize the bearing-vector to 1
      bearing = bearing / bearing.norm();

      bearingVectors.push_back(bearing);
      
      camCorrespondences.push_back(cam);
      
      // fill the indices map
      absoluteToLocal[points.size()-1] = std::make_pair(cam, i);
    }
  }

  assert(points.size() == numTotalPoints);
  assert(bearingVectors.size() == numTotalPoints);
  assert(camCorrespondences.size() == numTotalPoints);

  using namespace opengv;
  //create a non-central absolute adapter
  absolute_pose::NoncentralAbsoluteAdapter adapter(bearingVectors,
                                                   camCorrespondences,
                                                   points,
                                                   camOffsets,
                                                   camRotations);

  //Create a AbsolutePoseSacProblem and Ransac
  //The method is set to GP3P
  sac::Ransac<sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
  std::shared_ptr<
          sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
              new sac_problems::absolute_pose::AbsolutePoseSacProblem(
                  adapter,
                  sac_problems::absolute_pose::AbsolutePoseSacProblem::GP3P));  
  ransac.sac_model_ = absposeproblem_ptr;
  ransac.threshold_ = threshold;
  ransac.max_iterations_ = maxIterations;

  auto detect_start = std::chrono::steady_clock::now();
  const bool success = ransac.computeModel();
  auto detect_end = std::chrono::steady_clock::now();
  const auto detect_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(detect_end - detect_start);
  
  if(!success)
  {
    if(verbosity)
      std::cout << "Resection Failed" << std::endl;
    
    return false;
  }
  
  const auto &sol = ransac.model_coefficients_;
  // again, the output is a real pose, with the center and the rotation expressed
  // wrt the world frame, ie we need to take the transpose/inverse of the rotation
  rigPose = geometry::Pose3(sol.block<3,3>(0,0).transpose(), sol.col(3));

  // copy the inliers
  const std::size_t numInliers = ransac.inliers_.size();
  inliers.clear();
  inliers.resize(numCameras);
  
  if(verbosity)
  {
    std::cout << "\n"
    << "-------------------------------" << "\n"
    << "-- Robust Resection " << "\n"
    << "-- Resection status: " << success << "\n"
    << "-- #Points used for Resection: " << numTotalPoints << "\n"
    << "-- #Points validated by robust Resection: " << numInliers << "\n"
    << "-- #Iterations needed: " << ransac.iterations_ << "\n"
    << "-- #Thresehold used: " << ransac.threshold_ << "\n"
    << "-- Time spent in ransac [ms]: " << detect_elapsed.count() << "\n"
    << "-------------------------------" << std::endl;
  }
  
  // remap the inliers
  for(std::size_t i = 0; i < numInliers; i++)
  {
    const auto idx = absoluteToLocal.at(ransac.inliers_[i]);
    inliers[idx.first].emplace_back(idx.second);
  }

//  for(size_t i = 0; i < ransac.inliers_.size(); i++)
//    std::cout << ransac.inliers_[i] << " ";
//  std::cout << std::endl;
//  
//  for(size_t i = 0; i < inliers.size(); ++i)
//  {
//    std::cout << "Inliers cam " << i << ": ";
//    for(size_t j = 0; j < inliers[i].size(); ++j)
//      std::cout << inliers[i][j] << " ";
//    std::cout << std::endl;
//  }
  
  return success;
}

#endif // #if HAVE_OPENGV

}
}