// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERAS_COMMON_COMMAND_LINE_HELPER_HPP
#define OPENMVG_CAMERAS_CAMERAS_COMMON_COMMAND_LINE_HELPER_HPP

#include <iostream>
#include <string>
#include <vector>

#include "openMVG/stl/split.hpp"

namespace openMVG
{
namespace cameras
{

// Allow to initialize an object cameras::Intrinsic_Parameter_Type BA from
// a string and delimiters('|')
//

inline
cameras::Intrinsic_Parameter_Type
StringTo_Intrinsic_Parameter_Type
(
  const std::string & rhs
)
{
  // Split the string by the '|' token.
  std::vector<std::string> items;
  stl::split(rhs, '|', items);

  cameras::Intrinsic_Parameter_Type intrinsics_opt =
    static_cast<cameras::Intrinsic_Parameter_Type>(0);

  // Look for the "STRING KEY" parameters and initialize them
  for (const std::string & item : items)
  {
    // cameras::Intrinsic_Parameter_Type
    if (item == "NONE")
    {
      return cameras::Intrinsic_Parameter_Type::NONE;
    }
    else if (item == "ADJUST_FOCAL_LENGTH")
    {
      intrinsics_opt = intrinsics_opt | cameras::Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH;
    }
    else if (item == "ADJUST_PRINCIPAL_POINT")
    {
      intrinsics_opt = intrinsics_opt | cameras::Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT;
    }
    else if (item == "ADJUST_DISTORTION")
    {
      intrinsics_opt = intrinsics_opt | cameras::Intrinsic_Parameter_Type::ADJUST_DISTORTION;
    }
    else if (item == "ADJUST_ALL")
    {
      intrinsics_opt = cameras::Intrinsic_Parameter_Type::ADJUST_FOCAL_LENGTH
        | cameras::Intrinsic_Parameter_Type::ADJUST_PRINCIPAL_POINT
        | cameras::Intrinsic_Parameter_Type::ADJUST_DISTORTION;
    }
    else
    {
      std::cerr << "WARNING: Unknow KEY: " << item << std::endl;
      intrinsics_opt = static_cast<cameras::Intrinsic_Parameter_Type>(0);
      break;
    }
  }

  return intrinsics_opt;
}

} // namespace cameras
} // namespace openMVG

#endif // OPENMVG_CAMERAS_CAMERAS_COMMON_COMMAND_LINE_HELPER_HPP
