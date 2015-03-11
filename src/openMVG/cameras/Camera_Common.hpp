// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_COMMON_HPP
#define OPENMVG_CAMERAS_COMMON_HPP

enum EINTRINSIC
{
  PINHOLE_CAMERA_START,
  PINHOLE_CAMERA, // No distortion
  PINHOLE_CAMERA_RADIAL1, // radial distortion K1
  PINHOLE_CAMERA_RADIAL3, // radial distortion K1,K2,K3
  PINHOLE_CAMERA_BROWN, // radial + tangential
  PINHOLE_CAMERA_END
};

// Return if the camera type is a valid enum
static bool isValid(EINTRINSIC eintrinsic)
{
  return eintrinsic > PINHOLE_CAMERA_START && eintrinsic < PINHOLE_CAMERA_END;
}

#endif // OPENMVG_CAMERAS_COMMON_HPP
