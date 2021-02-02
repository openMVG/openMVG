// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_CAMERAS_CAMERAS_COMMON_COMMAND_LINE_HELPER_HPP
#define OPENMVG_CAMERAS_CAMERAS_COMMON_COMMAND_LINE_HELPER_HPP

#include <string>
#include <vector>

#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/stl/split.hpp"
#include "openMVG/system/logger.hpp"


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
      OPENMVG_LOG_ERROR << "WARNING: Unknow KEY: " << item;
      intrinsics_opt = static_cast<cameras::Intrinsic_Parameter_Type>(0);
      break;
    }
  }

  return intrinsics_opt;
}

} // namespace cameras

namespace sfm
{
  // Allow to initialize an object sfm::Extrinsic_Parameter_Type BA from
  // a string and delimiters('|')
  //

  inline
  sfm::Extrinsic_Parameter_Type
  StringTo_Extrinsic_Parameter_Type
  (
    const std::string & rhs
  )
  {
    // Split the string by the '|' token.
    std::vector<std::string> items;
    stl::split(rhs, '|', items);

    sfm::Extrinsic_Parameter_Type extrinsics_opt =
      static_cast<sfm::Extrinsic_Parameter_Type>(0);

    // Look for the "STRING KEY" parameters and initialize them
    for (const std::string & item : items)
    {
      // cameras::Intrinsic_Parameter_Type
      if (item == "NONE")
      {
        return sfm::Extrinsic_Parameter_Type::NONE;
      }
      else if (item == "ADJUST_TRANSLATION")
      {
        extrinsics_opt = extrinsics_opt | sfm::Extrinsic_Parameter_Type::ADJUST_TRANSLATION;
      }
      else if (item == "ADJUST_ROTATION")
      {
        extrinsics_opt = extrinsics_opt | sfm::Extrinsic_Parameter_Type::ADJUST_ROTATION;
      }
      else if (item == "ADJUST_ALL")
      {
        extrinsics_opt = sfm::Extrinsic_Parameter_Type::ADJUST_ALL;
      }
      else
      {
        OPENMVG_LOG_ERROR << "WARNING: Unknow KEY: " << item;
        extrinsics_opt = static_cast<sfm::Extrinsic_Parameter_Type>(0);
        break;
      }
    }
    return extrinsics_opt;
  }
} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_CAMERAS_CAMERAS_COMMON_COMMAND_LINE_HELPER_HPP
