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
 * @param[in,out] intrinsics The camera intrinsics
 * @param[in,out] localizationResult The series of camera poses and point correspondences 
 * @return true if the bundle adjustment has success
 */
bool refineSequence(cameras::Pinhole_Intrinsic_Radial_K3 *intrinsics,
                    std::vector<LocalizationResult> & localizationResult);

} //namespace localization
} //namespace openMVG

