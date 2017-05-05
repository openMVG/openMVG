// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_IO_HPP
#define OPENMVG_CAMERAS_CAMERA_IO_HPP

#include "openMVG/cameras/PinholeCamera.hpp"

namespace openMVG
{
namespace cameras
{


/**
* @brief Save a Pinhole camera to a file as a Projection matrix (12 doubles as binary values)
* @param scameraFile path to the file in which matrix will be written
* @param cam Camera to write
* @retval true If write is correct
* @retval false If there was an error during export
*/
bool save(
  const std::string & scameraFile,
  const PinholeCamera & cam );


/**
* @brief Load a Pinhole camera from a file (read 12 doubles saved in binary values)
* @param scameraFile path where camera is stored
* @param[out] cam output camera after read
* @retval true If loading is correct
* @retval false If there was an error during loading
*/
bool load(
  const std::string & scameraFile,
  PinholeCamera & cam );

} // namespace cameras
} // namespace openMVG

#endif // OPENMVG_CAMERAS_CAMERA_IO_HPP
