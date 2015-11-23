/* 
 * File:   optimization.hpp
 * Author: sgaspari
 *
 * Created on October 23, 2015, 12:02 AM
 */

#pragma once

#include "LocalizationResult.hpp"
#include <openMVG/cameras/Camera_Pinhole_Radial.hpp>

namespace openMVG{
namespace localization{

/**
 * @brief Refine a sequence of camera positions and the camera intrinsics.
 * The camera parameters are not changing during the sequence and they are refined
 * in terms of K (principal point and focal length) and the distortion coefficients.
 * The camera poses are refined as well.
 *  
 * @param[in,out] intrinsics The camera intrinsics.
 * @param[in,out] localizationResult The series of camera poses and point correspondences.
 * @param[in] b_refine_pose Whether to refine the camera poses.
 * @param[in] b_refine_intrinsic Whether to refine the camera parameters.
 * @param[in] b_refine_structure Whether to refine the 3D points.
 * @return true if the bundle adjustment has success.
 */
bool refineSequence(cameras::Pinhole_Intrinsic_Radial_K3 *intrinsics,
                    std::vector<LocalizationResult> & localizationResult,
                    bool b_refine_pose = true,
                    bool b_refine_intrinsic = true,
                    bool b_refine_structure = false);

/**
 * @brief Refine a sequence of camera positions and the camera intrinsics.
 * The camera parameters can change during the sequence and they are refined
 * in terms of K (principal point and focal length) and the distortion coefficients.
 * The camera poses are refined as well. If only one element is provided in the 
 * vector of intrinsics a single camera having constant internal parameters 
 * during the whole sequence is assumed.
 *  
 * @param[in,out] vec_intrinsics The vector containing the camera intrinsics. 
 * It can have the same size as \p vec_localizationResult or size equal to 1 if
 * the camera do not change its parameters during the sequence.
 * @param[in,out] vec_localizationResult The series of camera poses and point correspondences. 
 * @param[in] b_refine_pose Whether to refine the camera poses.
 * @param[in] b_refine_intrinsic Whether to refine the camera parameters.
 * @param[in] b_refine_structure Whether to refine the 3D points.
 * @return true if the bundle adjustment has success.
 */
bool refineSequence(std::vector<cameras::Pinhole_Intrinsic_Radial_K3* > vec_intrinsics,
                    std::vector<LocalizationResult> & vec_localizationResult,
                    bool b_refine_pose = true,
                    bool b_refine_intrinsic = true,
                    bool b_refine_structure = false);

} //namespace localization
} //namespace openMVG