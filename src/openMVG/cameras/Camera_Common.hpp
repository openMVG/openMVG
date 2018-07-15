// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERA_COMMON_HPP
#define OPENMVG_CAMERAS_CAMERA_COMMON_HPP

#include <type_traits>

namespace openMVG
{
namespace cameras
{

/**
* @enum EINTRINSIC List of usable camera Intrinsics
* @var PINHOLE_CAMERA
*   Pinhole camera is an ideal pinhole camera with 3x3 intrinsics matrix : \n
*      \f$ K=\begin{pmatrix} f & 0 & u_0 \\ 0 & f & v_0 \\ 0 & 0 & 1 \end{pmatrix} \f$
* @var PINHOLE_CAMERA_RADIAL1
*   Same as PINHOLE_CAMERA but before projection, pixel are distorded using radial distortion using one parameter \f$k_1\f$ \n
*    Assuming input pixel is \f$X\f$, distorded pixel \f$X_d\f$ is given by the relation : \n
*      \f$ X_d = ( 1 + k_1 r^2) X \f$  \n
*    Where \f$ r^2 = (X_x - u_0)^2 + (Y_y - v_0)^2 \f$
* @var PINHOLE_CAMERA_RADIAL3
*   Same as PINHOLE_CAMERA_RADIAL1 but using 3 parameters \f$k_1, k_2, k_3\f$ : \n
*      \f$ X_d = ( 1 + k_1 r^2 + k_2 r^4 + k_3 r^6 ) X \f$
* @var PINHOLE_CAMERA_BROWN
*   Same as PINHOLE_CAMERA with radial distortion and Tangential distortion : \n
*      \f$ x_d = x_u (1 + K_1 r^2 + K_2 r^4 + K_3 r^6) + (T_2 (r^2 + 2 x_u^2) + 2 T_1 x_u y_u) \f$
*      \f$ y_d = y_u (1 + K_1 r^2 + K_2 r^4 + K_3 r^6) + (T_1 (r^2 + 2 y_u^2) + 2 T_2 x_u y_u) \f$
* @var PINHOLE_CAMERA_FISHEYE
*   Simple fisheye camera with 4 distortion coefficients
*/
enum EINTRINSIC
{
  PINHOLE_CAMERA_START = 0,
  PINHOLE_CAMERA,         // No distortion
  PINHOLE_CAMERA_RADIAL1, // radial distortion K1
  PINHOLE_CAMERA_RADIAL3, // radial distortion K1,K2,K3
  PINHOLE_CAMERA_BROWN, // radial distortion K1,K2,K3, tangential distortion T1,T2
  PINHOLE_CAMERA_FISHEYE, // a simple Fish-eye distortion model with 4 distortion coefficients
  PINHOLE_CAMERA_END,
  CAMERA_SPHERICAL = PINHOLE_CAMERA_END + 1
};

/**
* @brief test if given intrinsic value corresponds to a pinhole
* @param eintrinsic Intrinsic value to test
* @retval true if parameter is a pinhole
* @retval false if parameter is not a pinhole
*/
static inline bool isPinhole( EINTRINSIC eintrinsic )
{
  return eintrinsic > PINHOLE_CAMERA_START && eintrinsic < PINHOLE_CAMERA_END;
}

static inline bool isSpherical( EINTRINSIC eintrinsic )
{
  return eintrinsic == CAMERA_SPHERICAL;
}

/**
* @brief Test if given intrinsic value is valid
* @param eintrinsic Intrinsic value to test
* @retval true if parameter is valid
* @retval false if parameter is invalid
*/
static inline bool isValid( EINTRINSIC eintrinsic )
{
  return isPinhole(eintrinsic) || isSpherical(eintrinsic);
}

/**
* @enum Intrinsic_Parameter_Type Used to control which camera parameter must be \n
   considered as variable of held constant for non linear refinement
* @var NONE
*   Intrinsic parameters will be considered as FIXED
* @var ADJUST_FOCAL_LENGTH
*   Focal length will be considered as variable for refinement
* @var ADJUST_PRINCIPAL_POINT
*   Principal point will be considered as variable for refinement
* @var ADJUST_DISTORTION
*   Distortion parameters will be considered as variable for refinement
* @var ADJUST_ALL
*   All parameters will be considered as variable for refinement
*/
enum class Intrinsic_Parameter_Type : int
{
  // Note: Use power of two values in order to use bitwise operators.
  NONE                    = 1, // All parameters will be held constant
  ADJUST_FOCAL_LENGTH     = 2,
  ADJUST_PRINCIPAL_POINT  = 4,
  ADJUST_DISTORTION       = 8,
  ADJUST_ALL = ADJUST_FOCAL_LENGTH | ADJUST_PRINCIPAL_POINT | ADJUST_DISTORTION
};

inline constexpr Intrinsic_Parameter_Type
operator|(Intrinsic_Parameter_Type x, Intrinsic_Parameter_Type y)
{
  return static_cast<Intrinsic_Parameter_Type>
    (static_cast<std::underlying_type<Intrinsic_Parameter_Type>::type>(x) |
     static_cast<std::underlying_type<Intrinsic_Parameter_Type>::type>(y));
}

inline constexpr Intrinsic_Parameter_Type
operator&(Intrinsic_Parameter_Type x, Intrinsic_Parameter_Type y)
{
  return static_cast<Intrinsic_Parameter_Type>
    (static_cast<std::underlying_type<Intrinsic_Parameter_Type>::type>(x) &
     static_cast<std::underlying_type<Intrinsic_Parameter_Type>::type>(y));
}

} // namespace cameras
} // namespace openMVG

#endif // OPENMVG_CAMERAS_CAMERA_COMMON_HPP
