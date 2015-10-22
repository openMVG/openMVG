#pragma once

//#include <vision/cameraTracking/markerBasedTracker/PonctualMarkerTracker.hpp>

#include <openMVG/localization/CCTagLocalizer.hpp>
#include <openMVG/geometry/pose3.hpp>
#include <openMVG/numeric/numeric.h>

#include <Eigen/Dense>

#include <vector>

namespace openMVG {
namespace rig {
  
class Rig {
public:
  
  Rig(){};
 
  virtual ~Rig();
  
#if 0
  std::size_t nCams() const { return _vTrackers.size(); }
  
  const PonctualMarkerTracker& tracker(std::size_t i) const { return _vTrackers[i]; }
  PonctualMarkerTracker& tracker(std::size_t i) { return _vTrackers[i]; }
  
  const Pose& relativePose(std::size_t i) const { return _vRelativePoses[i-1]; }
  Pose& relativePose(std::size_t i) { return _vRelativePoses[i-1]; } 
  
  const std::vector<Pose> & poses( ) const { return _vPoses; }
  std::vector<Pose>& poses( ) { return _vPoses; } 
  
  // Compute the intial relative positions between all witness cameras and the 
  // main one
  // ouput: true if the initial calibration succeed, false otherwise
  bool initializeCalibration();
  
  double distance(openMVG::Vec3 va, openMVG::Vec3 vb);
  
  Pose product(const Pose & poseA, const Pose & poseB);
  Pose productInv(const Pose & poseA, const Pose & poseB);
  
  // Optimize the initial solutions over all images
  // output: mean of the reprojection errors
  bool optimizeCalibration();
  
  // Direct method to compute an optimal solution to the relative pose problem for
  // the traker iTracker
  void findOptimalPose(const std::vector<Pose> & vPoses, std::size_t iTracker, Pose & result );
  
  void displayRelativePoseReprojection(const Pose & relativePose, std::size_t iTracker);
 #endif 
  
private:
  
  // Set of camera trackers (all of them associated to a camera)
  // THE FIRST index is associated to the main camera for which the pose
  // corresponds to the whole system pose.
  std::vector<localization::CCTagLocalizer> _vTrackers; // todo@L: replace by Tracker interface
  
  // (_vRelativePoses.size == nCams-1) where nCams represents the number of cameras including the main one
  // _vRelativePoses contains the relative poses of all witness cameras
  std::vector<geometry::Pose3> _vRelativePoses;
  // todo@L: discuss the "degree" of dependancy with openMVG
  // (BA aside in cameraLoc or used mainly the openMVG implementation?)
  
  std::vector<geometry::Pose3> _vPoses; // Pose of the rig, (== pose of the main camera == pose of _vTrackers[0])
};

#if 0

// Depreciated
// ugly -> required until type in camera localization change.
Pose toOMVG(const DeprePose & pose);

// Depreciated
// ugly -> required until type in camera localization change.
openMVG::Mat3 toOMVG(const numerical::BoundedMatrix3x3d & mat);

// Compute an average pose
void poseAveraging(const std::vector<Pose> & vPoses, Pose & result);

// Compute the sum of the reprojections error for imaged points imgPts of a set of points
// pts viewed through a camera whose pose Pose
double reprojectionError(const openMVG::Mat3 & K, const Pose & pose, const std::vector<openMVG::Vec2> & imgPts, const std::vector<openMVG::Vec3> & pts);

 // Compute the residual between the 3D projected point X and an image observation x
double residual( const openMVG::Mat3 & K, const Pose & pose, const openMVG::Vec3 & X, const openMVG::Vec2 & x);

// Compute the reprojection of X through the camera whose pose pose
void reproject( const openMVG::Mat3 & K, const Pose & pose, const openMVG::Vec3 & X, openMVG::Vec2 & res);

/* Apply a 2D projective transformation */
void projectiveTransform(const openMVG::Mat3 & H, const openMVG::Vec2 & x, openMVG::Vec2 & res);

#endif

} // namespace rig
} // namespace openMVG

