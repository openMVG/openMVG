// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_COMMON_HPP
#define OPENMVG_CAMERAS_COMMON_HPP

#include <string>
#include <stdexcept>

namespace openMVG {
namespace cameras {

enum EINTRINSIC
{
  PINHOLE_CAMERA_START = 0,
  PINHOLE_CAMERA = 1,         // No distortion
  PINHOLE_CAMERA_RADIAL1 = 2, // radial distortion K1
  PINHOLE_CAMERA_RADIAL3 = 3, // radial distortion K1,K2,K3
  PINHOLE_CAMERA_BROWN = 4, // radial distortion K1,K2,K3, tangential distortion T1,T2
  PINHOLE_CAMERA_FISHEYE = 5, // a simple Fish-eye distortion model with 4 distortion coefficients
  PINHOLE_CAMERA_END
};

inline std::string EINTRINSIC_enumToString(EINTRINSIC intrinsic)
{
  switch(intrinsic)
  {
    case PINHOLE_CAMERA:
      return "pinhole";
    case PINHOLE_CAMERA_RADIAL1:
      return "radial1";
    case PINHOLE_CAMERA_RADIAL3:
      return "radial3";
    case PINHOLE_CAMERA_BROWN:
      return "brown";
    case PINHOLE_CAMERA_FISHEYE:
      return "fisheye";
    case PINHOLE_CAMERA_START:
    case PINHOLE_CAMERA_END:
      break;
  }
  throw std::out_of_range("Invalid Intrinsic Enum");
}

inline EINTRINSIC EINTRINSIC_stringToEnum(const std::string& intrinsic)
{
  if(intrinsic == "pinhole")
    return PINHOLE_CAMERA;
  if(intrinsic == "radial1")
    return PINHOLE_CAMERA_RADIAL1;
  if(intrinsic == "radial3")
    return PINHOLE_CAMERA_RADIAL3;
  if(intrinsic == "brown")
    return PINHOLE_CAMERA_BROWN;
  if(intrinsic == "fisheye")
    return PINHOLE_CAMERA_FISHEYE;
  throw std::out_of_range(intrinsic);
}

// Return if the camera type is a valid enum
static inline bool isValid(EINTRINSIC eintrinsic)
{
  return eintrinsic > PINHOLE_CAMERA_START && eintrinsic < PINHOLE_CAMERA_END;
}

static inline bool isPinhole(EINTRINSIC eintrinsic)
{
  return eintrinsic > PINHOLE_CAMERA_START && eintrinsic < PINHOLE_CAMERA_END;
}

} // namespace cameras
} // namespace openMVG

#endif // OPENMVG_CAMERAS_COMMON_HPP
