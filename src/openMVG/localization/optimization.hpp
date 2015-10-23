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

bool refineSequence(cameras::Pinhole_Intrinsic_Radial_K3 *intrinsics,
                    std::vector<LocalizationResult> & localizationResult);

} //namespace localization
} //namespace openMVG

