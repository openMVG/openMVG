// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_UNDISTORT_IMAGE_HPP
#define OPENMVG_CAMERAS_CAMERA_UNDISTORT_IMAGE_HPP


#include "openMVG/cameras/Camera_Intrinsics.hpp"

/*
#include <limits>
#include <cmath>

#include "openMVG/cameras/Camera_Intrinsics.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/image/sample.hpp"
*/
namespace openMVG
{
namespace cameras
{

/**
* @brief  Undistort an image according a given camera & it's distortion model
* @param imageIn Input image
* @param cam Input intrinsic parameter used to undistort image
* @param[out] image_ud Output undistorted image
* @param fillcolor color used to fill pixels where no input pixel is found
*/
template <typename Image>
void UndistortImage(
  const Image& imageIn,
  const IntrinsicBase * cam,
  Image & image_ud,
  typename Image::Tpixel fillcolor = typename Image::Tpixel( 0 ) );

/**
* @brief  Undistort an image according a given camera & it's distortion model
* @param imageIn Input image
* @param cam Input intrinsic parameter used to undistort image
* @param[out] image_ud Output undistorted image
* @param fillcolor color used to fill pixels where no input pixel is found
* @param max_ud_width Maximum width of the undistorted image 
* @param max_ud_height Maximum height of the undistorted image 
* @note This function produces an image with size that try to fit the image plane
*/
template <typename Image>
void UndistortImageResized(
  const Image& imageIn,
  const IntrinsicBase * cam,
  Image & image_ud,
  typename Image::Tpixel fillcolor = typename Image::Tpixel( 0 ),
  const uint32_t max_ud_width = 10000,
  const uint32_t max_ud_height = 10000);


} // namespace cameras
} // namespace openMVG

#endif // #ifndef OPENMVG_CAMERAS_CAMERA_UNDISTORT_IMAGE_HPP

