/* 
 * File:   optimization.hpp
 * Author: sgaspari
 *
 * Created on October 23, 2015, 12:02 AM
 */

#pragma once

#include "LocalizationResult.hpp"
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>
#include <openMVG/geometry/pose3.hpp>

#include <vector>
#include <string>
#include <tuple>

namespace openMVG{
namespace localization{

/**
 * @brief Refine a sequence of camera positions and their camera intrinsics.
 * The camera parameters can change during the sequence and they are refined
 * in terms of K (principal point and focal length) and the distortion coefficients.
 * The camera poses are refined as well. If \p allTheSameIntrinsics is provided  
 * a single camera having constant internal parameters during the whole sequence is assumed.
 *  
 * @param[in,out] vec_localizationResult The series of camera poses and point correspondences. 
 * @param[in] allTheSameIntrinsics If true all the intrinsics of the sequence are
 * assumed to be the same, ie a sequence in which the camrera parameters do not change.
 * @param[in] b_refine_intrinsic Whether to refine the camera parameters.
 * @param[in] b_no_distortion If b_refine_intrinsic is true, this allow to not consider
 * the optical distortion, setting it to 0. 
 * @param[in] b_refine_pose Whether to refine the camera poses.
 * @param[in] b_refine_structure Whether to refine the 3D points.
 * @param[in] outputFilename If not empty, a filename (possibly with path) without 
 * extension where to save the scene before and after the bundle adjustment. The 
 * file will be named outputFilename.BEFORE.json and outputFilename.AFTER.json.
 * @param[in] minPointVisibility if > 0 it allows to use only the 3D points that 
 * are seen in at least \p minPointVisibility views/frames, all the other 
 * points (and associated 2D features) will be discarded.
 * @return true if the bundle adjustment has success.
 */
bool refineSequence(std::vector<LocalizationResult> & vec_localizationResult,
                    bool allTheSameIntrinsics = true,
                    bool b_refine_intrinsic = true,
                    bool b_no_distortion = false,
                    bool b_refine_pose = true,
                    bool b_refine_structure = false,
                    const std::string outputFilename = "",
                    std::size_t minPointVisibility = 0);

/**
 * @brief refine the pose of a camera rig by minimizing the reprojection error in
 * each camera with the bundle adjustment.
 * 
 * @param[in] vec_subPoses The rig calibration, ie the pose of the cameras wrt the main
 * camera; for a rig of N cameras the size of the vector is N-1 as the main camera has
 * identity as pose.
 * @param[in] vec_localizationResults A vector of N element containing the localization
 * results for each camera.
 * @param[in,out] rigPose The current rig pose and the refined rig pose if the bundle
 * adjustment succeeds.
 * @return true if the bundle adjustment succeeds.
 */
bool refineRigPose(const std::vector<geometry::Pose3 > &vec_subPoses,
                   const std::vector<localization::LocalizationResult> vec_localizationResults,
                   geometry::Pose3 & rigPose);

/**
 * @brief refine the pose of a camera rig of N cameras by minimizing the reprojection error in
 * each camera with the bundle adjustment.
 * 
 * @param[in] pts2d A vector of N 2xM matrices containing the image points of the 
 * 2D-3D associations.
 * @param[in] pts3d A vector of N 3xM matrices containing the image points of the 
 * 2D-3D associations.
 * @param[in] inliers A vector of N vectors, each contining the inliers for the 2D-3D
 * associations.
 * @param[in] vec_queryIntrinsics A vector of N intrinsics, one for each camera of the rig.
 * @param[in] vec_subPoses The rig calibration, ie the pose of the cameras wrt the main
 * camera; for a rig of N cameras the size of the vector is N-1 as the main camera has
 * identity as pose.
 * @param[in,out] rigPose rigPose The current rig pose and the refined rig pose if the bundle
 * adjustment succeeds.
 * @return true if the bundle adjustment succeeds.
 */
bool refineRigPose(const std::vector<Mat> &pts2d,
                   const std::vector<Mat> &pts3d,
                   const std::vector<std::vector<std::size_t> > &inliers,
                   const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                   const std::vector<geometry::Pose3 > &vec_subPoses,
                   geometry::Pose3 &rigPose);

/**
 * 
 * @param pts2d
 * @param pts3d
 * @param vec_queryIntrinsics
 * @param vec_subPoses
 * @param maxReprojectionError
 * @param minNumPoints
 * @param vec_inliers
 * @param rigPose
 * @param maxIterationNumber
 * @return 
 */
bool iterativeRefineRigPose(const std::vector<Mat> &pts2d,
                            const std::vector<Mat> &pts3d,
                            const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                            const std::vector<geometry::Pose3 > &vec_subPoses,
                            double maxReprojectionError,
                            std::size_t minNumPoints,
                            std::vector<std::vector<std::size_t> > &vec_inliers,
                            geometry::Pose3 &rigPose,
                            std::size_t maxIterationNumber = 10);

std::pair<double, bool> computeInliers(const std::vector<Mat> &vec_pts2d,
                                       const std::vector<Mat> &vec_pts3d,
                                       const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                       const std::vector<geometry::Pose3 > &vec_subPoses,
                                       const double maxReprojectionError,
                                       const geometry::Pose3 &rigPose,
                                       std::vector<std::vector<std::size_t> > &vec_inliers);

void printRigRMSEStats(const std::vector<Mat> &vec_pts2D,
                       const std::vector<Mat> &vec_pts3D,
                       const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                       const std::vector<geometry::Pose3 > &vec_subPoses,
                       const geometry::Pose3 &rigPose,
                       const std::vector<std::vector<std::size_t> > &vec_inliers);

/**
 * 
 * @param pts2d
 * @param pts3d
 * @param currCamera
 * @param currInliers
 * @param subPoses
 * @param rigPose
 * @return 
 */
std::tuple<double, double, double> computeStatistics(const Mat &pts2d, 
                                                     const Mat &pts3d,
                                                     const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera,
                                                     const std::vector<std::size_t> &currInliers,
                                                     const geometry::Pose3 &subPoses,
                                                     const geometry::Pose3 &rigPose);

} //namespace localization
} //namespace openMVG