// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_COMMON_HPP
#define OPENMVG_CAMERAS_COMMON_HPP

namespace openMVG {
namespace cameras {

enum EINTRINSIC
{
  PINHOLE_CAMERA_START = 0,
  PINHOLE_CAMERA = 1,         // No distortion
  PINHOLE_CAMERA_RADIAL1 = 2, // radial distortion K1
  PINHOLE_CAMERA_RADIAL3 = 3, // radial distortion K1,K2,K3
  PINHOLE_CAMERA_BROWN = 4, // radial distortion K1,K2,K3, tangential distortion T1,T2
  PINHOLE_CAMERA_END
};

// Return if the camera type is a valid enum
static bool isValid(EINTRINSIC eintrinsic)
{
  return eintrinsic > PINHOLE_CAMERA_START && eintrinsic < PINHOLE_CAMERA_END;
}

static bool isPinhole(EINTRINSIC eintrinsic)
{
  return eintrinsic > PINHOLE_CAMERA_START && eintrinsic < PINHOLE_CAMERA_END;
}

} // namespace cameras
} // namespace openMVG

#endif // OPENMVG_CAMERAS_COMMON_HPP
