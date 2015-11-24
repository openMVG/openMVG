#pragma once

#include <openMVG/localization/LocalizationResult.hpp>
#include <openMVG/geometry/pose3.hpp>
#include <openMVG/numeric/numeric.h>

#include <Eigen/Dense>

#include <vector>
#include <map>

#ifdef HAVE_OPENCV
//#define VISUAL_DEBUG_MODE
#endif

namespace openMVG {
namespace rig {
  
class Rig {
public:
  
  Rig(){};
 
  virtual ~Rig();
  
  // Accessors
  std::size_t nCams() const { return _vLocalizationResults.size(); }
  
  std::vector<localization::LocalizationResult> & getLocalizationResults(IndexT i)
  { 
    return _vLocalizationResults[i];
  }
  
  const geometry::Pose3& getRelativePose(std::size_t i) const { return _vRelativePoses[i-1]; }
  
  const std::vector<geometry::Pose3> & getPoses( ) const { return _vPoses; }
  
  /*
   * @brief Compute the initial guess related to the relative positions between all witness 
   * cameras and the main camera
   * @return True if the initial calibration succeed, false otherwise.
   */
  bool initializeCalibration();
  
  void setTrackingResult(
          std::vector<localization::LocalizationResult> vLocalizationResults,
          std::size_t i);
  
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

  /*
   * @brief Perform the rig bundle adjustment
   */
  bool optimizeCalibration();
  
  /**
   * @brief Save the calibrated poses to a text file.
   * @param filename The filename for the calibration file.
   * @return true if everything went ok.
   * @see saveRigCalibration()
   */
  bool saveCalibration(std::string &filename);
  
  /*
   * @brief Visual debug function displaying the reprojected 3D points and their
   * associated observation.
   */
  void displayRelativePoseReprojection(const geometry::Pose3 & relativePose, std::size_t iTracker);
  
private:
  // Set of localization results (each of them associated to a camera)
  // The FIRST INDEX is associated to the MAIN CAMERA for which the pose
  // corresponds to the entire system pose.
  std::map<IndexT,std::vector<localization::LocalizationResult> > _vLocalizationResults;
  
  // (_vRelativePoses.size == nCams-1) where nCams represents the number of cameras
  // including the main one
  // _vRelativePoses contains the relative poses of all witness cameras
  std::vector<geometry::Pose3> _vRelativePoses;
  
  // Rig pose
  std::vector<geometry::Pose3> _vPoses; // (i.e., by convention, pose of the main camera)
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
 * @brief Compute the witness camera from the main camera pose and the relative pose 
 * from the main camera to the witness camera
 * 
 * @param[in] poseMainCamera Pose of the main camera
 * @param[in] relativePose Relative pose from the main camera to the witness camera
 * @return The absolute pose of the witness camera
 */
geometry::Pose3 poseFromMainToWitness(geometry::Pose3 poseMainCamera, geometry::Pose3 relativePose);

geometry::Pose3 computeRelativePose(geometry::Pose3 poseMainCamera, geometry::Pose3 poseWitnessCamera);

/*
 * @brief Visual debug function doing a pause during the program execution.
 */
void cvpause();

/**
 * @brief Load the set of subposes from a simple text file.
 * @param[in] filename The file from which to load the subposes.
 * @param[out] subposes The loaded subposes.
 * @return true if everything went ok.
 * 
 * The format of the file is the following
 * numCam
 * R[0][0] // first camera rotation
 * R[0][1]
 * ...
 * t[0] // first camera translation
 * t[1]
 * t[2]
 * R[0][0] // second camera rotation
 * ...
 */
bool loadRigCalibration(const std::string &filename, std::vector<geometry::Pose3> &subposes);

/**
 * @brief Save the set of subposes from a simple text file.
 * @param[in] filename The file to which the subposes are saved.
 * @param[in] subposes The subposes to write.
 * @return true if everything went ok.
 * 
 * The format of the file is the following
 * numCam
 * R[0][0] // first camera rotation
 * R[0][1]
 * ...
 * t[0] // first camera translation
 * t[1]
 * t[2]
 * R[0][0] // second camera rotation
 * ...
 */
bool saveRigCalibration(const std::string &filename, const std::vector<geometry::Pose3> &subposes);

} // namespace rig
} // namespace openMVG

