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


bool refineRigPose(const std::vector<geometry::Pose3 > &vec_subPoses,
                   const std::vector<localization::LocalizationResult> vec_localizationResults,
                   geometry::Pose3 & rigPose);

} //namespace localization
} //namespace openMVG