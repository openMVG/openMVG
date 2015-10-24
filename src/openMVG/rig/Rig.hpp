#pragma once

//#include <vision/cameraTracking/markerBasedTracker/PonctualMarkerTracker.hpp>

#include <openMVG/localization/LocalizationResult.hpp>
#include <openMVG/geometry/pose3.hpp>
#include <openMVG/numeric/numeric.h>

#include <Eigen/Dense>

#include <vector>
#include <map>

namespace openMVG {
namespace rig {
  
class Rig {
public:
  
  Rig(){};
 
  virtual ~Rig();
  
  std::size_t nCams() const { return _vLocalizationResults.size(); }
  
  std::vector<localization::LocalizationResult> & getLocalizationResults(IndexT i)
  { 
    return _vLocalizationResults[i];
  }
  //localization::LocalizationResult& getLocalizationResults(std::size_t i) { return _vLocalizationResults[i]; }
  
  const geometry::Pose3& getRelativePose(std::size_t i) const { return _vRelativePoses[i-1]; }
  //geometry::Pose3& getRelativePose(std::size_t i) { return _vRelativePoses[i-1]; } 
  
  const std::vector<geometry::Pose3> & getPoses( ) const { return _vPoses; }
  //std::vector<geometry::Pose3>& getPoses( ) { return _vPoses; }
  
  /*
   * @brief Compute the initial guess related to the relative positions between all witness 
   * cameras and the main camera
   * @return True if the initial calibration succeed, false otherwise.
   */
  bool initializeCalibration();
  
  /*
   * @brief 
   */
  void setTrackingResult(
          std::vector<localization::LocalizationResult> vLocalizationResults,
          std::size_t i);
  
  //double distance(openMVG::Vec3 va, openMVG::Vec3 vb);
  
  //Pose product(const Pose & poseA, const Pose & poseB);
 // Pose productInv(const Pose & poseA, const Pose & poseB);
  
  /*
   * @brief From a set of relative poses, find the one for a given 
   *        localization result which minimizes the reprojection errors over
   *        all views.
   * 
   * @param[in] vPoses Relative poses to test
   * @param[in] Index of the considered witness camera
   * @param[out] Best relative pose belonging in vPoses
   */
  void findBestRelativePose(
        const std::vector<geometry::Pose3> & vPoses,
        std::size_t iRes,
        geometry::Pose3 & result );

  // Optimize the initial solutions over all images
  // output: mean of the reprojection errors
  //bool optimizeCalibration();
  
  //void displayRelativePoseReprojection(const Pose & relativePose, std::size_t iTracker);
  
private:
  // Set of localization results (all of them associated to a camera)
  // The FIRST INDEX is associated to the MAIN CAMERA for which the pose
  // corresponds to the entire system pose.
  std::map<IndexT,std::vector<localization::LocalizationResult> > _vLocalizationResults;
  
  // (_vRelativePoses.size == nCams-1) where nCams represents the number of cameras
  // including the main one
  // _vRelativePoses contains the relative poses of all witness cameras
  std::vector<geometry::Pose3> _vRelativePoses;
  
  // Rig pose
  std::vector<geometry::Pose3> _vPoses; // (i.e. pose of the main camera)
};

/*
 * @brief For a given localization result, compute the sum of the reprojection errors
 * (over all points) related to another pose than the one previously computed and store
 * in the provided localizationResult instance.
 * 
 * @param[in] localizationResult The localization result
 * @param[in] pose The pose
 * @return The reprojection error over all inliers store in localizationResult
 */
double reprojectionError(const localization::LocalizationResult & localizationResult, const geometry::Pose3 & pose);

/*
 *
 */
geometry::Pose3 computeRelativePose(geometry::Pose3 poseMainCamera, geometry::Pose3 poseWitnessCamera);
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

